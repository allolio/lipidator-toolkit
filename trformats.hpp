#ifndef HAVE_CONFIG
#include "config.hpp"
#endif
// Oracle XDR is necessary
#include <rpc/xdr.h>
#define XTC_MAGIC 1995
#define XTC_HEADER_SIZE 16

typedef float xtc_cor;

#define MAXABS INT_MAX-2
#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

//GPL ON

static const int magicints[] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0,
    8, 10, 12, 16, 20, 25, 32, 40, 50, 64,
    80, 101, 128, 161, 203, 256, 322, 406, 512, 645,
    812, 1024, 1290, 1625, 2048, 2580, 3250, 4096, 5060, 6501,
    8192, 10321, 13003, 16384, 20642, 26007, 32768, 41285, 52015, 65536,
    82570, 104031, 131072, 165140, 208063, 262144, 330280, 416127, 524287, 660561,
    832255, 1048576, 1321122, 1664510, 2097152, 2642245, 3329021, 4194304, 5284491, 6658042,
    8388607, 10568983, 13316085, 16777216 };

#define FIRSTIDX 9
/* note that magicints[FIRSTIDX-1] == 0 */
#define LASTIDX (sizeof(magicints) / sizeof(*magicints))

/*_________________________________________________________________________
 |
 | sizeofint - calculate bitsize of an integer
 |
 | return the number of bits needed to store an integer with given max size
 |
*/

static int sizeofint(const int size) {
    unsigned int num = 1;
    int num_of_bits = 0;
    
    while (size >= num && num_of_bits < 32) {
	num_of_bits++;
	num <<= 1;
    }
    return num_of_bits;
}

/*___________________________________________________________________________
 |
 | sizeofints - calculate 'bitsize' of compressed ints
 |
 | given the number of small unsigned integers and the maximum value
 | return the number of bits needed to read or write them with the
 | routines receiveints and sendints. You need this parameter when
 | calling these routines. Note that for many calls I can use
 | the variable 'smallidx' which is exactly the number of bits, and
 | So I don't need to call 'sizeofints for those calls.
*/

static int sizeofints( const int num_of_ints, unsigned int sizes[]) {
    int i, num;
    unsigned int num_of_bytes, num_of_bits, bytes[32], bytecnt, tmp;
    num_of_bytes = 1;
    bytes[0] = 1;
    num_of_bits = 0;
    for (i=0; i < num_of_ints; i++) {	
	tmp = 0;
	for (bytecnt = 0; bytecnt < num_of_bytes; bytecnt++) {
            /* BUG: overflow? */
	    tmp = bytes[bytecnt] * sizes[i] + tmp;
	    bytes[bytecnt] = tmp & 0xff;
	    tmp >>= 8;
	}
	while (tmp != 0) {
	    bytes[bytecnt++] = tmp & 0xff;
	    tmp >>= 8;
	}
	num_of_bytes = bytecnt;
    }
    num = 1;
    num_of_bytes--;
    while (bytes[num_of_bytes] >= num) {
	num_of_bits++;
	num *= 2;
    }
    return num_of_bits + num_of_bytes * 8;
}


/*___________________________________________________________________________
 |
 | receivebits - decode number from buf using specified number of bits
 | 
 | extract the number of bits from the array buf and construct an integer
 | from it. Return that value.
 |
*/

static int receivebits(int buf[], int num_of_bits) {

    int cnt, num; 
    unsigned int lastbits, lastbyte;
    unsigned char * cbuf;
    int mask = (1 << num_of_bits) -1;

    cbuf = ((unsigned char *)buf) + 3 * sizeof(*buf);
    cnt = buf[0];
    lastbits = (unsigned int) buf[1];
    lastbyte = (unsigned int) buf[2];
    
    num = 0;
    while (num_of_bits >= 8) {
	lastbyte = ( lastbyte << 8 ) | cbuf[cnt++];
	num |=  (lastbyte >> lastbits) << (num_of_bits - 8);
	num_of_bits -=8;
    }
    if (num_of_bits > 0) {
	if (lastbits < num_of_bits) {
	    lastbits += 8;
	    lastbyte = (lastbyte << 8) | cbuf[cnt++];
	}
	lastbits -= num_of_bits;
	num |= (lastbyte >> lastbits) & ((1 << num_of_bits) -1);
    }
    num &= mask;
    buf[0] = cnt;
    buf[1] = lastbits;
    buf[2] = lastbyte;
    return num; 
}

/*____________________________________________________________________________
 |
 | receiveints - decode 'small' integers from the buf array
 |
 | this routine is the inverse from sendints() and decodes the small integers
 | written to buf by calculating the remainder and doing divisions with
 | the given sizes[]. You need to specify the total number of bits to be
 | used from buf in num_of_bits.
 |
*/

static void receiveints(int buf[], const int num_of_ints, int num_of_bits,
	unsigned int sizes[], int nums[]) {
    int bytes[32];
    int i, j, num_of_bytes, p, num;
    
    bytes[1] = bytes[2] = bytes[3] = 0;
    num_of_bytes = 0;
    while (num_of_bits > 8) {
	bytes[num_of_bytes++] = receivebits(buf, 8);
	num_of_bits -= 8;
    }
    if (num_of_bits > 0) {
	bytes[num_of_bytes++] = receivebits(buf, num_of_bits);
    }
    for (i = num_of_ints-1; i > 0; i--) {
	num = 0;
	for (j = num_of_bytes-1; j >=0; j--) {
	    num = (num << 8) | bytes[j];
	    p = num / sizes[i];
	    bytes[j] = p;
	    num = num - p * sizes[i];
	}
	nums[i] = num;
    }
    nums[0] = bytes[0] | (bytes[1] << 8) | (bytes[2] << 16) | (bytes[3] << 24);
}

/*____________________________________________________________________________
 |
 | xdr3dfcoord - read or write compressed 3d coordinates to xdr file.
 |
 | this routine reads or writes (depending on how you opened the file with
 | xdropen() ) a large number of 3d coordinates (stored in *fp).
 | The number of coordinates triplets to write is given by *size. On
 | read this number may be zero, in which case it reads as many as were written
 | or it may specify the number if triplets to read (which should match the
 | number written).
 | Compression is achieved by first converting all floating numbers to integer
 | using multiplication by *precision and rounding to the nearest integer.
 | Then the minimum and maximum value are calculated to determine the range.
 | The limited range of integers so found, is used to compress the coordinates.
 | In addition the differences between succesive coordinates is calculated.
 | If the difference happens to be 'small' then only the difference is saved,
 | compressing the data even more. The notion of 'small' is changed dynamically
 | and is enlarged or reduced whenever needed or possible.
 | Extra compression is achieved in the case of GROMOS and coordinates of
 | water molecules. GROMOS first writes out the Oxygen position, followed by
 | the two hydrogens. In order to make the differences smaller (and thereby
 | compression the data better) the order is changed into first one hydrogen
 | then the oxygen, followed by the other hydrogen. This is rather special, but
 | it shouldn't harm in the general case.
 |
 */



int xdr3dfcoord(XDR *xdrs, float *fp, int *size, float *precision)
{
    int     *ip  = NULL;
    int     *buf = NULL;
    bool bRead;

    /* preallocate a small buffer and ip on the stack - if we need more
       we can always malloc(). This is faster for small values of size: */
    unsigned     prealloc_size = 3*16;
    int          prealloc_ip[3*16], prealloc_buf[3*20];
    int          we_should_free = 0;

    int          minint[3], maxint[3], mindiff, *lip, diff;
    int          lint1, lint2, lint3, oldlint1, oldlint2, oldlint3, smallidx;
    int          minidx, maxidx;
    unsigned     sizeint[3], sizesmall[3], bitsizeint[3], size3, *luip;
    int          flag, k;
    int          smallnum, smaller, larger, i, is_small, is_smaller, run, prevrun;
    float       *lfp, lf;
    int          tmp, *thiscoord,  prevcoord[3];
    unsigned int tmpcoord[30];

    int          bufsize, xdrid, lsize;
    unsigned int bitsize;
    float        inv_precision;
    int          errval = 1;
    int          rc;

    bRead         = (xdrs->x_op == XDR_DECODE);
    bitsizeint[0] = bitsizeint[1] = bitsizeint[2] = 0;
    prevcoord[0]  = prevcoord[1]  = prevcoord[2]  = 0;



        /* xdrs is open for reading */

        if (xdr_int(xdrs, &lsize) == 0)
        {
            return 0;
        }
        if (*size != 0 && lsize != *size)
        {
            fprintf(stderr, "wrong number of coordinates in xdr3dfcoord; "
                    "%d arg vs %d in file", *size, lsize);
        }
        *size = lsize;
        size3 = *size * 3;
        if (*size <= 9)
        {
            *precision = -1;
            return (xdr_vector(xdrs, (char *) fp, (unsigned int)size3,
                               (unsigned int)sizeof(*fp), (xdrproc_t)xdr_float));
        }
        if (xdr_float(xdrs, precision) == 0)
        {
            return 0;
        }

        if (size3 <= prealloc_size)
        {
            ip  = prealloc_ip;
            buf = prealloc_buf;
        }
        else
        {
            we_should_free = 1;
            bufsize        = size3 * 1.2;
            ip             = (int *)malloc((size_t)(size3 * sizeof(*ip)));
            buf            = (int *)malloc((size_t)(bufsize * sizeof(*buf)));
            if (ip == NULL || buf == NULL)
            {
                fprintf(stderr, "malloc failed\n");
                exit(1);
            }
        }

        buf[0] = buf[1] = buf[2] = 0;

        if ( (xdr_int(xdrs, &(minint[0])) == 0) ||
             (xdr_int(xdrs, &(minint[1])) == 0) ||
             (xdr_int(xdrs, &(minint[2])) == 0) ||
             (xdr_int(xdrs, &(maxint[0])) == 0) ||
             (xdr_int(xdrs, &(maxint[1])) == 0) ||
             (xdr_int(xdrs, &(maxint[2])) == 0))
        {
            if (we_should_free)
            {
                free(ip);
                free(buf);
            }
            return 0;
        }

        sizeint[0] = maxint[0] - minint[0]+1;
        sizeint[1] = maxint[1] - minint[1]+1;
        sizeint[2] = maxint[2] - minint[2]+1;

        /* check if one of the sizes is to big to be multiplied */
        if ((sizeint[0] | sizeint[1] | sizeint[2] ) > 0xffffff)
        {
            bitsizeint[0] = sizeofint(sizeint[0]);
            bitsizeint[1] = sizeofint(sizeint[1]);
            bitsizeint[2] = sizeofint(sizeint[2]);
            bitsize       = 0; /* flag the use of large sizes */
        }
        else
        {
            bitsize = sizeofints(3, sizeint);
        }

        if (xdr_int(xdrs, &smallidx) == 0)
        {
            if (we_should_free)
            {
                free(ip);
                free(buf);
            }
            return 0;
        }

        maxidx       = MIN(LASTIDX, smallidx + 8);
        minidx       = maxidx - 8; /* often this equal smallidx */
        smaller      = magicints[MAX(FIRSTIDX, smallidx-1)] / 2;
        smallnum     = magicints[smallidx] / 2;
        sizesmall[0] = sizesmall[1] = sizesmall[2] = magicints[smallidx];
        larger       = magicints[maxidx];

        /* buf[0] holds the length in bytes */

        if (xdr_int(xdrs, &(buf[0])) == 0)
        {
            if (we_should_free)
            {
                free(ip);
                free(buf);
            }
            return 0;
        }


        if (xdr_opaque(xdrs, (char *)&(buf[3]), (unsigned int)buf[0]) == 0)
        {
            if (we_should_free)
            {
                free(ip);
                free(buf);
            }
            return 0;
        }



        buf[0] = buf[1] = buf[2] = 0;

        lfp           = fp;
        inv_precision = 1.0 / *precision;
        run           = 0;
        i             = 0;
        lip           = ip;
        while (i < lsize)
        {
            thiscoord = (int *)(lip) + i * 3;

            if (bitsize == 0)
            {
                thiscoord[0] = receivebits(buf, bitsizeint[0]);
                thiscoord[1] = receivebits(buf, bitsizeint[1]);
                thiscoord[2] = receivebits(buf, bitsizeint[2]);
            }
            else
            {
                receiveints(buf, 3, bitsize, sizeint, thiscoord);
            }

            i++;
            thiscoord[0] += minint[0];
            thiscoord[1] += minint[1];
            thiscoord[2] += minint[2];

            prevcoord[0] = thiscoord[0];
            prevcoord[1] = thiscoord[1];
            prevcoord[2] = thiscoord[2];


            flag       = receivebits(buf, 1);
            is_smaller = 0;
            if (flag == 1)
            {
                run        = receivebits(buf, 5);
                is_smaller = run % 3;
                run       -= is_smaller;
                is_smaller--;
            }
            if (run > 0)
            {
                thiscoord += 3;
                for (k = 0; k < run; k += 3)
                {
                    receiveints(buf, 3, smallidx, sizesmall, thiscoord);
                    i++;
                    thiscoord[0] += prevcoord[0] - smallnum;
                    thiscoord[1] += prevcoord[1] - smallnum;
                    thiscoord[2] += prevcoord[2] - smallnum;
                    if (k == 0)
                    {
                        /* interchange first with second atom for better
                         * compression of water molecules
                         */
                        tmp          = thiscoord[0]; thiscoord[0] = prevcoord[0];
                        prevcoord[0] = tmp;
                        tmp          = thiscoord[1]; thiscoord[1] = prevcoord[1];
                        prevcoord[1] = tmp;
                        tmp          = thiscoord[2]; thiscoord[2] = prevcoord[2];
                        prevcoord[2] = tmp;
                        *lfp++       = prevcoord[0] * inv_precision;
                        *lfp++       = prevcoord[1] * inv_precision;
                        *lfp++       = prevcoord[2] * inv_precision;
                    }
                    else
                    {
                        prevcoord[0] = thiscoord[0];
                        prevcoord[1] = thiscoord[1];
                        prevcoord[2] = thiscoord[2];
                    }
                    *lfp++ = thiscoord[0] * inv_precision;
                    *lfp++ = thiscoord[1] * inv_precision;
                    *lfp++ = thiscoord[2] * inv_precision;
                }
            }
            else
            {
                *lfp++ = thiscoord[0] * inv_precision;
                *lfp++ = thiscoord[1] * inv_precision;
                *lfp++ = thiscoord[2] * inv_precision;
            }
            smallidx += is_smaller;
            if (is_smaller < 0)
            {
                smallnum = smaller;
                if (smallidx > FIRSTIDX)
                {
                    smaller = magicints[smallidx - 1] /2;
                }
                else
                {
                    smaller = 0;
                }
            }
            else if (is_smaller > 0)
            {
                smaller  = smallnum;
                smallnum = magicints[smallidx] / 2;
            }
            sizesmall[0] = sizesmall[1] = sizesmall[2] = magicints[smallidx];
        }
    
    if (we_should_free)
    {
        free(ip);
        free(buf);
    }
    return 1;
}


 
int readatoms(XDR *xdrs,  vector<Atom> &vatoms, int *size, float *precision) {
    

    static int *ip = NULL;
    static int oldsize;
    static int *buf;

    int minint[3], maxint[3], mindiff, *lip, diff;
    int lint1, lint2, lint3, oldlint1, oldlint2, oldlint3, smallidx;
    int minidx, maxidx;
    unsigned sizeint[3], sizesmall[3], bitsizeint[3], size3, *luip;
    int flag, k;
    int smallnum, smaller, larger, i, is_small, is_smaller, run, prevrun;
    float *lfp;
    int atomindex=0;
    int tmp, *thiscoord,  prevcoord[3];
    unsigned int tmpcoord[30];

    int bufsize, lsize;
    unsigned int bitsize;
    float inv_precision;
    int errval = 1;

    bitsizeint[0] = bitsizeint[1] = bitsizeint[2] = 0;
    prevcoord[0]  = prevcoord[1]  = prevcoord[2]  = 0;
    
   	
	/* xdrs is open for reading */
	
	if (xdr_int(xdrs, &lsize) == 0) 
	    return 0;
	if (*size != 0 && lsize != *size) {
	    fprintf(stderr, "wrong number of coordinates in frame; "
		    "%d arg vs %d in file", *size, lsize);
	}
	*size = lsize;
	size3 = *size * 3;
	//if (*size <= 9) {
        //    return (xdr_vector(xdrs, (char *) fp, (unsigned int)size3, 
        //            (unsigned int)sizeof(*fp), (xdrproc_t)xdr_float));/
	// }
	xdr_float(xdrs, precision);
	if (ip == NULL) {
	    ip = (int *)malloc((size_t)(size3 * sizeof(*ip)));
	    if (ip == NULL) {
		fprintf(stderr,"malloc failed\n");
		exit(1);
	    }
	    bufsize = size3 * 1.2;
	    buf = (int *)malloc((size_t)(bufsize * sizeof(*buf)));
	    if (buf == NULL) {
		fprintf(stderr,"malloc failed\n");
		exit(1);
	    }
	    oldsize = *size;
	} else if (*size > oldsize) {
	    ip = (int *)realloc(ip, (size_t)(size3 * sizeof(*ip)));
	    if (ip == NULL) {
		fprintf(stderr,"malloc failed\n");
		exit(1);
	    }
	    bufsize = size3 * 1.2;
	    buf = (int *)realloc(buf, (size_t)(bufsize * sizeof(*buf)));
	    if (buf == NULL) {
		fprintf(stderr,"malloc failed\n");
		exit(1);
	    }
	    oldsize = *size;
	}
	buf[0] = buf[1] = buf[2] = 0;
	
	xdr_int(xdrs, &(minint[0]));
	xdr_int(xdrs, &(minint[1]));
	xdr_int(xdrs, &(minint[2]));

	xdr_int(xdrs, &(maxint[0]));
	xdr_int(xdrs, &(maxint[1]));
	xdr_int(xdrs, &(maxint[2]));
		
	sizeint[0] = maxint[0] - minint[0]+1;
	sizeint[1] = maxint[1] - minint[1]+1;
	sizeint[2] = maxint[2] - minint[2]+1;
	
	/* check if one of the sizes is to big to be multiplied */
	if ((sizeint[0] | sizeint[1] | sizeint[2] ) > 0xffffff) {
	    bitsizeint[0] = sizeofint(sizeint[0]);
	    bitsizeint[1] = sizeofint(sizeint[1]);
	    bitsizeint[2] = sizeofint(sizeint[2]);
	    bitsize = 0; /* flag the use of large sizes */
	} else {
	    bitsize = sizeofints(3, sizeint);
	}
	
	if (xdr_int(xdrs, &smallidx) == 0)	
	    return 0;
	maxidx = MIN(LASTIDX, smallidx + 8) ;
	minidx = maxidx - 8; /* often this equal smallidx */
	smaller = magicints[MAX(FIRSTIDX, smallidx-1)] / 2;
	smallnum = magicints[smallidx] / 2;
	sizesmall[0] = sizesmall[1] = sizesmall[2] = magicints[smallidx] ;
	larger = magicints[maxidx];

    	/* buf[0] holds the length in bytes */

	if (xdr_int(xdrs, &(buf[0])) == 0)
	    return 0;
        if (xdr_opaque(xdrs, (char *)&(buf[3]), (unsigned int)buf[0]) == 0)
	    return 0;
	buf[0] = buf[1] = buf[2] = 0;
	
	inv_precision = 1.0 / * precision;
	run = 0;
	i = 0;
	lip = ip;
	while ( i < lsize ) {
	  atomindex=i;
	    thiscoord = (int *)(lip) + i * 3;

	    if (bitsize == 0) {
		thiscoord[0] = receivebits(buf, bitsizeint[0]);
		thiscoord[1] = receivebits(buf, bitsizeint[1]);
		thiscoord[2] = receivebits(buf, bitsizeint[2]);
	    } else {
		receiveints(buf, 3, bitsize, sizeint, thiscoord);
	    }
	    
	    i++;
	    thiscoord[0] += minint[0];
	    thiscoord[1] += minint[1];
	    thiscoord[2] += minint[2];
	    
	    prevcoord[0] = thiscoord[0];
	    prevcoord[1] = thiscoord[1];
	    prevcoord[2] = thiscoord[2];
	    
	   
	    flag = receivebits(buf, 1);
	    is_smaller = 0;
	    if (flag == 1) {
		run = receivebits(buf, 5);
		is_smaller = run % 3;
		run -= is_smaller;
		is_smaller--;
	    }
	    if (run > 0) {
		thiscoord += 3;
		for (k = 0; k < run; k+=3) {
		    receiveints(buf, 3, smallidx, sizesmall, thiscoord);
		    i++;
		    thiscoord[0] += prevcoord[0] - smallnum;
		    thiscoord[1] += prevcoord[1] - smallnum;
		    thiscoord[2] += prevcoord[2] - smallnum;
		    if (k == 0) {
			/* interchange first with second atom for better
			 * compression of water molecules
			 */
			tmp = thiscoord[0]; thiscoord[0] = prevcoord[0];
				prevcoord[0] = tmp;
			tmp = thiscoord[1]; thiscoord[1] = prevcoord[1];
				prevcoord[1] = tmp;
			tmp = thiscoord[2]; thiscoord[2] = prevcoord[2];
				prevcoord[2] = tmp;
			vatoms[atomindex].pos.x = prevcoord[0] * inv_precision*10;
			vatoms[atomindex].pos.y = prevcoord[1] * inv_precision*10;
			vatoms[atomindex].pos.z= prevcoord[2] * inv_precision*10;
		    } else {
			prevcoord[0] = thiscoord[0];
			prevcoord[1] = thiscoord[1];
			prevcoord[2] = thiscoord[2];
		    }
		    vatoms[atomindex].pos.x = thiscoord[0] * inv_precision*10;
		    vatoms[atomindex].pos.y = thiscoord[1] * inv_precision*10;
		    vatoms[atomindex].pos.z = thiscoord[2] * inv_precision*10;
		}
	    } else {
		vatoms[atomindex].pos.x = thiscoord[0] * inv_precision*10;
		vatoms[atomindex].pos.y = thiscoord[1] * inv_precision*10;
		vatoms[atomindex].pos.z = thiscoord[2] * inv_precision*10;		
	    }
	    smallidx += is_smaller;
	    if (is_smaller < 0) {
		smallnum = smaller;
		if (smallidx > FIRSTIDX) {
		    smaller = magicints[smallidx - 1] /2;
		} else {
		    smaller = 0;
		}
	    } else if (is_smaller > 0) {
		smaller = smallnum;
		smallnum = magicints[smallidx] / 2;
	    }
	    sizesmall[0] = sizesmall[1] = sizesmall[2] = magicints[smallidx] ;
	}
    
    return 1;
}
// GPL OFF

class XTCFile : public TrajecFile
{   
  public:
    XTCFile(string filename) {
    hasloaded=false;frames=0;fpos=0;lsnaps.clear(); lpos=lsnaps.begin(); 
    readcache=false;hasnames=false;
    Load(filename);
    }
     XTCFile() {
    hasloaded=false;frames=0;fpos=0;lsnaps.clear(); lpos=lsnaps.begin(); 
    readcache=false;fsize=0;
    }
    ~XTCFile()
    {
    if(hasloaded)
    {
    xdr_destroy(&xdr); 
    fclose(file); }
    }
     bool Load(string filename);
     bool SaveTo(string Filename);
     bool setAtomTypes(vector<string> types);
     void setAtomType(atomlist& indices, string name);
     bool GroAtomTypes(string grofile,int column);
  protected:
     bool ParseNextFrame();
     bool GotoFrame(int n);
     int filebybytes();
     bool hasnames;
     vector<string> anames;
     bool read_cor(xtc_cor *coor);
     bool forwardframe();
     XDR xdr; // Numerical Conversion Class
     FILE* file; // C style FIO
     off64_t fsize;
};

bool XTCFile::ParseNextFrame()
{
  if(readcache)  { GotoFrame(fpos);
#ifdef DEBUG
cout << "Last Frame from Cache!" << endl;
#endif
readcache=false;}

 int magic;   
 if (xdr_int(&xdr,&magic) == 0) return false;
  if(magic!=XTC_MAGIC) {cerr << "Magic Number is wrong" << endl;return false;}
  xdr_int(&xdr,&atoms);
  int frdump;/* number of atoms */
  xdr_int(&xdr,&frdump);    /* frame number    */
  xtc_cor time;
  read_cor(&time);
vector<Atom> vatoms;
vatoms.resize(atoms+3);
  for(int i=0;i!=3;i++)
  {
    xtc_cor a=0;
	vatoms[atoms+i].type="PBC";
	read_cor(&a);
	vatoms[atoms+i].pos.x=a*10;
		read_cor(&a);
	vatoms[atoms+i].pos.y=a*10;
		read_cor(&a);
	vatoms[atoms+i].pos.z=a*10;
//    cout << "READ PBC" <<  a << endl ;
  }
  // Check Header
  float precision=0;
  float *fp=(float*) malloc(sizeof(float)*3*atoms);
  float* lfp=fp;
  xdr3dfcoord(&xdr, fp, &atoms, &precision);
  for(int i=0;i!=atoms;i++)
  {
    vatoms[i].pos.x=*lfp*10; lfp++;
    vatoms[i].pos.y=*lfp*10; lfp++;
    vatoms[i].pos.z=*lfp*10; lfp++;
  }
  free(fp);
//  readatoms(&xdr, vatoms, &atoms, &precision);
  if(hasnames)
  {
    for(int i=0;i!=atoms;i++)
     vatoms[i].type=anames[i];
  }
  SnapShot snap;
 lsnaps.insert(lpos,make_pair( snap,fpos));
   lpos--; 
   ((*lpos).first).SnapOverwrite(vatoms,false,true);
   fpos++;
  
 return true;
}

bool XTCFile::GotoFrame(int n)
{
//  cout << "FPOS" << fpos << " N" << n << endl;
  if(n>frames) return false;
    if(fpos>n && hasloaded) { fseek(file,0,SEEK_SET); fpos=0; return true; }
    else if(!hasloaded) return false;
  while(fpos<n)
  {
    forwardframe();
    fpos++;
  } 
  return true;
}

bool XTCFile::SaveTo(string filename)
{
  cerr << "NYI - Not yet implemented" << endl;
  return false;
}

bool XTCFile::setAtomTypes(vector<string> atomtypes)
{
  anames=atomtypes;
  hasnames=true;
  return true;
}

void XTCFile::setAtomType(atomlist &index, string name)
{
  for(int i=0; i!=index.size(); i++)
   {
     anames[index[i]]=name;
   }
   hasnames=true;
   return;   
}


bool XTCFile::Load(string filename)
{
  file = 0;  
  file = fopen(filename.c_str(), "rb"); 
  xdrstdio_create(&xdr, file, XDR_DECODE);
 
  if(file == NULL) {
    return false;
  }
  int magic;
  int step;
  // Read Header
    int result;

  if (xdr_int(&xdr,&magic) == 0)
    return 0;

  if(magic!=XTC_MAGIC) {cerr << "Magic Number is wrong" << endl;return false;}
  xdr_int(&xdr,&atoms);  /* number of atoms */
  xdr_int(&xdr,&frames);    /* frame number    */
  xtc_cor time;
  read_cor(&time);
  anames.resize(atoms);
  if(!hasnames) for(int i=0;i!=atoms;i++) anames[i]="XX";
  fseek(file,0,SEEK_SET);
  hasloaded=true;
  // Box vectors
  frames=filebybytes();
  hasvelocities=false;
  fpos=0;
  //
  return true;
  
}

bool XTCFile::read_cor(xtc_cor *coor)
{
//(float *
//(*coor)=0;
  if( xdr_float(&xdr,coor) == 0) return false ;
  return true;
}


bool XTCFile::forwardframe()
{
 int magic;
 off64_t currpos= ftell(file);
 int increment=atoms*sizeof(int)/8*4+XTC_HEADER_SIZE;
  off64_t currpos2=currpos+increment;  
  fseek(file,currpos2,SEEK_SET);
  while(xdr_int(&xdr,&magic)!=0 && currpos2<fsize)
  { currpos2+=4;
    if(magic==XTC_MAGIC)
    {
      int a1;
      xdr_int(&xdr,&a1);
      if(a1==atoms)
      {
	 xdr_int(&xdr,&a1);
         if(a1>1)
	 {
	 currpos2-=4;
	 }
	  fseek(file,currpos2,SEEK_SET);
         return true;
	 
      }
    }
   }  
  return false;
}

bool XTCFile::GroAtomTypes(string Filename,int column)
{
    ifstream input(Filename.c_str());

  string line;
  getline(input,line);
  getline(input,line);
  int limit=atoi(line.c_str());
  if(atoms==0) anames.resize(limit);
  for(int i=0;i!=anames.size();i++)
  {
  getline(input,line);
     vector<string> toks;
   Tokenize(line,&toks, " \t");
   anames[i]=toks[column];
   if(i>9998 && i<100000) anames[i].resize(anames[i].size()-5);
   hasnames=true;
//   cout << anames[i] << endl;
  }
  input.close();
  return true;
}

int  XTCFile::filebybytes()
{
  frames=1;
  if(!hasloaded) return -1;
  fseek(file,0,SEEK_END);
  fsize= ftell(file);
  fseek(file,0,SEEK_SET);
  while(forwardframe()) frames++;
      fseek(file,0,SEEK_SET);
  return frames;
}
