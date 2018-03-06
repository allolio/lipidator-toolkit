#ifndef HAVE_CONFIG
#include "config.hpp"
#endif
typedef vector<int> atomlist;

#include "trformats.hpp"
//#include "lipidrecsum.hpp"
#include "histogram.hpp"
#include "kiss_fft130/tools/kiss_fftnd.h"

void buf2yvals(kiss_fft_cpx *buf, vector<yvals> &in)
{
   int ydim=in[0].size();
   for (int i=0;i!=in.size();i++)
  {
       for (int j=0;j!=ydim;j++)
       {
          in[i][j].real(buf[i*ydim+j].r);
	  in[i][j].imag(buf[i*ydim+j].i);
       }
  }    
}

void hist2yvals(Histogram<double> &h, vector<yvals> &out)
{
    Grid *g=h.getGrid();
    out.clear();
     for(int i=0;i!=g->limits[1];i++)
   {
     yvals yv;
     for(int j=0;j!=g->limits[0];j++)
     {
       complex<double> dp=h.GetRBin(j,i);
       if(dp!=dp ) dp=0.0;
       yv.push_back(dp); 
     }
     out.push_back(yv);
     //cout <<":" << yv.size() << " " << endl;
   }
}

void yval2hist(Histogram<double> &h, vector<yvals> &yvmat)
{
         double xstep=h.getGrid()->lengths[0]*((double) h.getGrid()->limits[0]/(double) yvmat[0].size());
     double ystep=h.getGrid()->lengths[1]*((double) h.getGrid()->limits[1]/(double) yvmat.size());
   for(int i=0;i!=yvmat.size();i++)
   {
     for(int j=0;j!=yvmat[0].size();j++)
     {
     h.AddDatapoint((yvmat[i])[j].real(),(xstep*(j+0.5)),(ystep*(i+0.5))) ;
     // cout  << h.getGrid()->lengths[0] << "," << endl;
  //    if (abs((yvmat[i])[j]) < 0.001 ) ((yvmat[i])[j])=0.0;
   //    cout  << abs((oldmat[i])[j]) << "/"<< abs((outmat[i])[j]) << " ";   
     }
   // cout << endl;
   }
}

void yvals2buf(vector<yvals> &in, kiss_fft_cpx *inb)
{
      int ydim=in[0].size();

   for (int i=0;i!=in.size();i++)
  {
       for (int j=0;j!=ydim;j++)
       {
	 inb[i*ydim+j].r=in[i][j].real();
	 inb[i*ydim+j].i=in[i][j].imag();
       }
  }    
}


void fft_surface(vector<yvals> &in, vector<yvals> &out, bool inv=0)
{
  out=in;
  kiss_fft_cpx inb[in[0].size()*in.size()]; 
  kiss_fft_cpx outb[in[0].size()*in.size()]; 
  yvals2buf(in,inb);
   int dims[2]; dims[1]=in[0].size(),dims[0]=in.size();
  kiss_fftnd_cfg  cfg=kiss_fftnd_alloc ( dims,2, inv,0,0);
  kiss_fftnd(cfg,(kiss_fft_cpx*) inb,  (kiss_fft_cpx*) outb);
  buf2yvals(outb,out);
  free(cfg);
}

double gfilter(double qabs)
{
  double retval=1;
  retval=exp(-(qabs)*(qabs));
  return retval;
}

void fftshift(vector<yvals> &fftarray, vector<yvals> *out)
{
  int xgrid=fftarray[0].size();
  int ygrid=fftarray.size();
  for(int i=0;i!=ygrid/2;i++)
  {
    for(int j=0;j!=xgrid/2;j++)
    {
      (out->at(i+ygrid/2)[j+xgrid/2])=fftarray[i][j];
      out->at(i)[j]=fftarray[(i+ ygrid/2)][j+xgrid/2];
      out->at(i)[xgrid/2+j]=fftarray[(i+ygrid/2)][j];
      out->at(ygrid/2+i)[j]=fftarray[i][j+xgrid/2];
      //
    }
  }
}


void zpad ( vector<yvals> &ymat,int zeros)
{
//    int zeros=200;
  vector<yvals> yout=ymat;
  double xh=ymat[0].size()/2.0;
  double yh=ymat.size()/2.0;
    fft_surface(yout,ymat);

  fftshift(ymat,&yout);

    ymat.resize(yout.size()+2*zeros);
    yvals bar; bar.resize(yout[0].size()+2*zeros,0.0);
    yvals b0;b0.resize(zeros,0.0);
    int inc=0;
    for(int i=0;i!=ymat.size();i++)
    {
        if(i<zeros || i >= (ymat.size()-zeros))
        {
            inc++;
            ymat[i]=bar;
            cout << i  << "BAM!" <<  ymat.size() << " " << inc << endl;
        }
        else {
            ymat[i]=b0;
            ymat[i].insert(ymat[i].end(),yout[i-inc].begin(),yout[i-inc].end());
           // ymat[i]=yout[i-inc];
            ymat[i].insert(ymat[i].end(),b0.begin(),b0.end());
            
            }
    }
       yout=ymat;
    cout << "done" << endl;
  fftshift(ymat,&yout);
  // ymat=yout;
//    yout.resize(ymat.size());
  
    //yout.push_back(bar);
    cout << endl;
    cout << "INVERT" << endl;
   fft_surface(yout,ymat,1);
// ymat=yout;
//   return;
    for(int x=0;x!=yout[0].size();x++)
    for(int y=0;y!=yout.size();y++)
    {     
      ymat[y][x].real()/=((yh*xh*2.0)*2.0);
//      cout << ymat[y][x].real() << " " ;
      //<< " " ;///1302.0<< " "; 
    }
    cout << endl;
   cout  << " " << xh << " "<< yh << endl;
 //  ymat=yout;
}

void lpass ( vector<yvals> &ymat)
{
  vector<yvals> yout=ymat;
  double xh=ymat[0].size()/2.0;
  double yh=ymat.size()/2.0;
    fft_surface(yout,ymat);

   fftshift(ymat,&yout);

     for(int x=0;x!=yout[0].size();x++)
    for(int y=0;y!=yout.size();y++)
    {     
      double qabs=sqrt(pow((x-xh)/xh,2)+pow((y-yh)/yh,2));
      double fac=gfilter(qabs*15.5);
      //cout << fac << endl;
//        int cutoff=yout[0].size()-40;
        int cutoff=2;
   // if(x > cutoff &&  (x <yout[0].size()-cutoff ) ) fac=0;
     //    else if (y > 2  && y < yout.size()-2) fac=0;
   // if((x+y)> 100 && (x-yout.size()-cutoff +y < 100)) fac=0;
        (yout[y])[x]=yout[y][x]*fac;

      //ymat[0][0]=0;
      //ymat[0][0]=0;
      //;
      //if(qabs > 50) fac=0;
     ///(ymat.size()*ymat[0].size());
 //     (ymat[2*yh-y][x])*=fac; 
//     (ymat[y])[2*xh-x]*=fac; 
 //    (ymat[2*yh-y])[2*xh-x]*=fac; 
//     cout << " X" << x << "Y" << y<< ymat[y][x] << "XY" << ymat[yh+y][x] << "YHX" << ymat[2*yh-1-y][x]  << "FINX" << endl;
    }
    ymat=yout;
    cout << "done" << endl;
  fftshift(ymat,&yout);
  // ymat=yout;
//    yout.resize(ymat.size());
  
    //yout.push_back(bar);
    cout << endl;
    cout << "INVERT" << endl;
   fft_surface(yout,ymat,1);
// ymat=yout;
    for(int x=0;x!=yout[0].size();x++)
    for(int y=0;y!=yout.size();y++)
    {     
      ymat[y][x].real()/=(yh*xh)*4.0;
    }
    cout << endl;
   cout  << " " << xh << " "<< yh << endl;
 //  ymat=yout;
}



//#include "diffgeom.hpp"
//#include "funcfit.hpp"



template <class T> void FromFile(string fname,Histogram<T > &h, int dim=3)
{
  FILE* f=fopen(fname.c_str(),"r");
  fseek(f,0,SEEK_SET);
  if (f==NULL) return;
  int rv=1;
  float x,y,z;
  while(rv==1)
  {
 fscanf(f,"%g",&x);
   fscanf(f,"%g",&y);
 rv= fscanf(f,"%g",&z);
 x/=10;
 y/=10;
 z/=10;
// if (x >0 && y>0 && z>4 )
 if(dim==3)
     h.AddDatapoint((T) z,(double) x, (double) y);
 else if (dim==4)
     h.AddDatapoint(1,(double) x, (double) y, (double) z);
   
  }
//cout << x << " "<< y << " "<< z << endl;  
  //   h.AddDatapoint((double) z,(double) x, (double) y);

}

void testfunc(Histogram< double >  &h,int xlim, int ylim)
{
  // 
  for(int x=0;x!=2000.0;x++)
    for(int y=0;y!=2000.0;y++)
    {
     double x1=(double)(x)*2.0*PI/(double) 2000.0; 
     double x2=(double)y *2.0*PI/(double) 2000.0;
     double f=cos(x1)*cos(2*x2)+6;
     h.AddDatapoint(f,x/2000.0*xlim,y/2000.0*ylim);
    }
}


int main(int argc, char **argv) {

//    t : 0.420722 B: -0.0240242 T: -0.0220562 DL: 0.0286894 t(old): 1.0742 p  90.6488 20.8745 69.2924

 
Histogram<double>s_in; // Surface
Histogram<double>s_fin; // Surface
Histogram<double>s_cover; // Curvature
Histogram<double>s_out; // Outsurface
Histogram<double>s_outc; // Outcover
int zeros=100;
ifstream input;
input.open(argv[1]);
/* Cutoffs
 * [x], [y], [z]
 * Bins 
 * [x],[y],[z]
 *
 */
double xmin,xmax,ymin,ymax,zmin,zmax;
int xgrid,ygrid,zgrid;
double linethickness;
/*
 * GUESS
*/
xmin=atof(argv[2]);
xmax=atof(argv[3]);
ymin=atof(argv[4]);
ymax=atof(argv[5]);
zmin=atof(argv[6]);
zmax=atof(argv[7]);

xgrid=(xmax-xmin);
ygrid=(ymax-ymin);
zgrid=(zmax-zmin)/.5;
//cout << ymin <<" "<< ymax << endl;
/** Reading grid**/
Grid g2;
Grid g3,g4;

g2.limits.push_back((int) (xgrid));
g2.limits.push_back((int) (ygrid));
g2.width.push_back(make_pair<double,double>(xmin,xmax));
g2.width.push_back(make_pair<double,double>(ymin,ymax));

g3.limits.push_back((int) (xgrid/4.0));
g3.limits.push_back((int) (ygrid/4.0));
g3.limits.push_back(zgrid);//32 12
g3.width.push_back(make_pair<double,double>(xmin,xmax));
g3.width.push_back(make_pair<double,double>(ymin,ymax));
g3.width.push_back(make_pair<double,double>(0.0,zmax-zmin));

g4.limits.push_back((int) (xgrid+2*zeros));
g4.limits.push_back((int) (ygrid+2*zeros));
g4.width.push_back(make_pair<double,double>(xmin,xmax));
g4.width.push_back(make_pair<double,double>(ymin,ymax));


/**
 * Writing Grid
**/
Grid g;
g.limits.push_back(xgrid+2*zeros);
g.limits.push_back(ygrid+2*zeros);
g.limits.push_back(zgrid);//32 12
g.width.push_back(make_pair<double,double>(0.0,xmax-xmin));
g.width.push_back(make_pair<double,double>(0.0,ymax-ymin));
g.width.push_back(make_pair<double,double>(0.0,zmax-zmin));

s_in.setGrid(g2);
s_fin.setGrid(g4);
s_cover.setGrid(g2);
s_out.setGrid(g);
s_outc.setGrid(g3);
XYZFile xyz("middle.xyz");
SnapShot s=xyz.GetCurrentSnap();
ofstream of;
of.open("atoms.dat");
of << std::fixed;
      of << std::setw( 12 );
      of << std::setprecision(6);
for(int i=0;i!=s.GetAtoms();i++)
{
   Atom &at=s.GetAtom(i);
   tuple a=at.pos;
  // if(a.x> xmin && a.y > zmin && a.z > ymin && a.z < ymax && a.x < xmax && a.y < zmax)
   {
       int type;
   if(at.type=="C") type=6;   
      if(at.type=="H") type=1;    
   if(at.type=="O") type=8;    
   if(at.type=="N") type=7;    
if(at.type=="P") type=15; 
   of << type << " "<< 0.0 << " " << (a.x-xmin)*1.88973 << " " << (a.z-ymin)*1.88973 << " " << (a.y-zmin)*1.88973 << endl;
   }
}
of.close();
string buffer;
int col=4;
while(getline(input,buffer)>0)
{
   vector<string> line;
   int htokens=Tokenize(buffer,&line, " \t");
 //  cout << line[12] << " " << line[14] << endl;
  // cout << line[col] << endl;
   tuple n,dir;
   n.x=atof(line[20].c_str());n.y=atof(line[21].c_str());n.z=atof(line[22].c_str());
   dir.x=0.0;dir.y=1.0;dir.z=0;
          n.print(); 

   if((atof(line[13].c_str()) > zmin) && (atof(line[13].c_str()) < zmax) && (dir*n > 0))
   {
   s_in.AddDatapoint( atof(line[13].c_str()), atof(line[12].c_str()), atof(line[14].c_str()));
   if(abs(atof(line[col].c_str())) < 1)
   s_cover.AddDatapoint( atof(line[col].c_str()), atof(line[12].c_str()), atof(line[14].c_str()));
   //cout << line[14] << " "  << line[13] << " "<< line[12]<< endl;

   }
//   h.AddDatapoint( atof(line[col].c_str()), atof(line[14].c_str()), atof(line[13].c_str()), atof(line[12].c_str()));
}
input.close();
//testfunc(s_in,xmax,ymax);
vector<yvals> surf;
hist2yvals(s_in,surf);
zpad(surf,zeros);
lpass(surf);
yval2hist(s_fin,surf);

vector<vector <double> > data;
vector<vector <double> > datac;
s_fin.Data(&data);
s_cover.Data(&datac);
cout << data[0].size()<< endl;
for(int i=0;i!=data[0].size();i++)
{
    vector<double> point;
    point.push_back(data[0][i]);    point.push_back(data[1][i]);
    point.push_back(data[2][i]);
    s_out.AddDatapoint(1.0,point[0]-xmin,point[1]-ymin,point[2]-zmin);
  //  cout << point[0]-xmin << " "  << point[1]-ymin<< " "<< point[2]-zmin<< endl;
    for(int j=0;j!=zgrid;j++)
    {
//   cout << datac[0][i] << " " << datac[2][i] << endl;
   if(i< datac[0].size())     s_outc.AddDatapoint(datac[2][i],datac[0][i]-xmin,datac[1][i]-ymin,((zmax-zmin)/(double) zgrid)*(0.5+j));
        s_out.AddDatapoint(1.0*((double) j/zgrid),point[0]-xmin,point[1]-ymin,(point[2]-zmin)/double(zgrid)*(j+1.0));
  
    }
    
}
ofstream of2;
of2.open("debug.dat");

s_cover.DropData(of2);
of2.close();
of2.open("sfdebug.dat");
s_fin.DropData(of2);
of2.close();
of2.open("surf.cube");
s_out.DropCube(of2);
//h.DropData(of2);
of2.close();
of2.open("mark.cube");
s_outc.DropCube(of2);
s_outc.DropData(cout);
//h.DropData(of2);
of2.close();

}
