#ifndef HAVE_CONFIG
#include "config.hpp"
#endif
typedef vector<int> atomlist;
#include <complex.h>
#include "kiss_fft130/tools/kiss_fftnd.h"
#include "trformats.hpp"
#include "lipidrecsum.hpp"
#include "histogram.hpp"
#define NFFT_PRECISION_DOUBLE

#include "funcfit.hpp"
//#include "surf.h"

int sign(double x)
{
  return x/abs(x);
}
void xydfft( vector<yvals> &yvmat, vector<yvals> &outmat,bool rev=false)
{
  
  outmat=yvmat;
  // Step1 Row Transform
  for(int i=0;i!=yvmat.size();i++)
  {
    if(!rev)
    yvmat[i]=r2fft(yvmat[i]);
    else  {  yvmat[i]=ir2fft(yvmat[i]);}
  }
  
  
  // Step 2 Column Transform // Transpose
  for(int i=0;i!=yvmat.size();i++)
  {    
    yvals oldcolumn,newcolumn;
    oldcolumn.resize(yvmat.size()); newcolumn.resize(yvmat.size());
    for(int j=0;j!=yvmat.size();j++)
    { 
      oldcolumn[j]=(yvmat[j][i]);
      
    }
    if(!rev) newcolumn=r2fft(oldcolumn);
      else {  newcolumn=ir2fft(oldcolumn);}
  
 
    outmat[i]=newcolumn;
   }
  
  
  return;
}
void bin_surface(double *grid, int gdim, double*in,int M, vector<bool> &updown, bool up=1)
{
    double *weights=(double*) malloc(gdim*gdim*sizeof(double));
  memset(weights,0,gdim*gdim*sizeof(double));
  memset(grid,0,gdim*gdim*sizeof(double));
  // bin is in center : epsilon = 0.5*width;
  // Values -0.5:0.5 --> 0:1
  // [ ] [X ] 
  // [ ] [ ]  [ . ] [X . ] weight= wa+wb = 1  wb = 1- abs(x-center)/width : wa = 1-wb 
 
  double binwidth=1.0/(double) gdim;
  
  for(int i=0; i!=updown.size();i++)
  {
    if(updown[i]==up)
    {
    double x=in[i*3]+0.5;
    double value=in[i*3+2];
    double y=in[i*3+1]+0.5;    
   // cout << x << " " << y<<  " "  << endl;
    if(x>=1) x=x-1;
    if(y>=1) y=y-1;
    if(x<0) x=1+x;
    if(y<0) y=1+y;

    // Find adjacent bins;
    int binx=(x)/binwidth;
    int biny=(y)/binwidth;     
   // cout << (x/ (double) binwidth)-binx-0.5 << endl;
    int adjx=binx+sign((x/binwidth)-binx-0.5);
    if(adjx>gdim-1) adjx=0;    
    if(adjx<0) adjx=gdim-1;    
    int adjy=biny+sign((y/binwidth)-biny-0.5);
    if(adjy>gdim-1) adjy=0;    
    if(adjy<0) adjy=gdim-1;
   // cout << x <<" " << y << " " << binx<<" " << biny <<" "<< adjx<< " " << adjy<< endl; // adjx / biny // binx / adjy // ajdx / adjy /
        double dx=abs(x-(binx+0.5)*(binwidth))/binwidth;
	double dy=abs(y-(biny+0.5)*(binwidth))/binwidth;
    double alpha=dx;
    double beta=dy;
    double gamma=sqrt(dx*dy);
    double omega=1/(1+alpha+beta+gamma);
    double wx=omega*alpha;
    double wy=omega*beta;
    double wxy=omega*gamma;
 // cout << wx+wy+wxy +omega << " a " << alpha << " b " << beta << " g " << gamma<< "  o "  << omega<< endl;
 //   cout << alpha+beta+gamma+omega << " dx " << dx << " dy " << dy << " g " << gamma<< "  o "  << omega<< endl;
    grid[biny*gdim+binx]+=value*omega;
    grid[biny*gdim+adjx]+=value*wx;
    grid[adjy*gdim+binx]+=value*wy;
    grid[adjy*gdim+adjx]+=value*wxy;
    weights[biny*gdim+binx]+=omega;
    weights[biny*gdim+adjx]+=wx;
    weights[adjy*gdim+binx]+=wy;
    weights[adjy*gdim+adjx]+=wxy;}    
  }
 // cout << "DONE FILLIN'" << endl;
 vector <int> missing;
  missing.clear();
   for(int i=0; i!=gdim*gdim;i++)
  {
    if( weights[i] !=0 ){ grid[i]/=weights[i]; }
    else {cerr << "WARNING! No data"  << endl; missing.push_back(i);} 
  }
  for(int i=0;i!=missing.size();i++)
  {
    int x=missing[i]%gdim;
    int y=missing[i]/gdim;
    double adjy=grid[((y+1)%gdim)*gdim+(x)%gdim];
    double adjx=grid[((y)%gdim)*gdim+(x+1)%gdim];
    if (x==0) x=gdim;
    if (y==0) y=gdim;
    double adjmx=grid[y*gdim+(x-1)];
    double adjmy=grid[(y-1)*gdim+(x)];
    grid[missing[i]]=0.25*(adjx+adjy+adjmx+adjmy);
  }
  // Stretch !
  free(weights);

}

void bin_surface_old(double *grid, int gdim, double*in, int M,vector<bool> &updown, bool up)
{
    double *weights=(double*) malloc(gdim*gdim*sizeof(double));
  memset(weights,0,gdim*gdim*sizeof(double));
  memset(grid,0,gdim*gdim*sizeof(double));
  // bin is in center : epsilon = 0.5*width;
  // Values -0.5:0.5 --> 0:1
  // [ ] [X ] 
  // [ ] [ ]  [ . ] [X . ] weight= wa+wb = 1  wb = 1- abs(x-center)/width : wa = 1-wb 
 
  double binwidth=1.0/(double) gdim;
  
  for(int i=0; i!=M;i++)
  {
    if(updown[i]==up){
    double x=in[i*3]+0.5;
    double value=in[i*3+2];
    double y=in[i*3+1]+0.5;    
   // cout << x << " " << y<<  " "  << endl;
    if(x>=1) x=x-1;
    if(y>=1) y=y-1;
    if(x<0) x=1+x;
    if(y<0) y=1+y;

    // Find adjacent bins;
    int binx=(x)/binwidth;
    int biny=(y)/binwidth;     
   // cout << (x/ (double) binwidth)-binx-0.5 << endl;
    int adjx=binx+sign((x/binwidth)-binx-0.5);
    if(adjx>gdim-1) adjx=0;    
    if(adjx<0) adjx=gdim-1;    
    int adjy=biny+sign((y/binwidth)-biny-0.5);
    if(adjy>gdim-1) adjy=0;    
    if(adjy<0) adjy=gdim-1;
   // cout << x <<" " << y << " " << binx<<" " << biny <<" "<< adjx<< " " << adjy<< endl; // adjx / biny // binx / adjy // ajdx / adjy /
        double dx=abs(x-(binx+0.5)*(binwidth))/binwidth;
	double dy=abs(y-(biny+0.5)*(binwidth))/binwidth;
    double alpha=dx;
    double beta=dy;
    double gamma=sqrt(dx*dy);
    double omega=1/(1+alpha+beta+gamma);
    double wx=omega*alpha;
    double wy=omega*beta;
    double wxy=omega*gamma;
 // cout << wx+wy+wxy +omega << " a " << alpha << " b " << beta << " g " << gamma<< "  o "  << omega<< endl;
 //   cout << alpha+beta+gamma+omega << " dx " << dx << " dy " << dy << " g " << gamma<< "  o "  << omega<< endl;
 
    grid[biny*gdim+binx]+=value*omega;
    grid[biny*gdim+adjx]+=value*wx;
    grid[adjy*gdim+binx]+=value*wy;
    grid[adjy*gdim+adjx]+=value*wxy;
    weights[biny*gdim+binx]+=omega;
    weights[biny*gdim+adjx]+=wx;
    weights[adjy*gdim+binx]+=wy;
    weights[adjy*gdim+adjx]+=wxy;    
  // cout << "DONE FILLIN'" << endl;
    }}
   for(int i=0; i!=gdim*gdim;i++)
  {
    if( weights[i] !=0 ){ grid[i]/=weights[i]; }
    else {cerr << "WARNING! No data"  << endl; if(i>0) grid[i]=grid[i-1];} 
  }
  
  // Stretch !
  free(weights);

}
void AssignLeaflet(vector<bool> &updown, vector <tuple> *dir)
{
  tuple plane;
  plane.x=0;
  plane.y=0;
  plane.z=1;
  updown.resize(dir->size());
  for(int i=0;i!=dir->size();i++)
  {
    double sp=dir->at(i)*plane;
    if(sp >0 ) updown[i] =true;
    else updown[i]=0;
  }
}



vector<yvals> interpol_e(int N,int M, double *in,vector<bool> &updown)
{
  // in : xyz array // out: complex double
  //N output grid
  //M number of points
  int gdim=sqrt(M/2);
  //if(!isPower2(gdim)) gdim=nextp2(gdim/2);
  double *grid=(double*) malloc(gdim*gdim*sizeof(double));
  int dims[2]; dims[0]=gdim,dims[1]=gdim;
  kiss_fft_cpx inb[gdim*gdim]; 
  kiss_fft_cpx out[gdim*gdim];
  kiss_fft_cpx out2[gdim*gdim];
  bin_surface(grid,gdim,in,M,updown,true); 
  for (int i=0;i!=gdim;i++)
  {
       for (int j=0;j!=gdim;j++)
       {
	 inb[i*gdim+j].r=grid[i*gdim+j];
	 inb[i*gdim+j].i=0;
	// cout << inc*(i+0.5)-0.5<<" " << inc*(j+0.5)-0.5 << " "   << grid[i*gdim+j] << endl;
       }
  }
  kiss_fftnd_cfg  cfg=kiss_fftnd_alloc ( dims,2, 0,0,0);
  kiss_fftnd(cfg,(kiss_fft_cpx*) inb,  (kiss_fft_cpx*) out);
  bin_surface(grid,gdim,in,M,updown,false); 
  for (int i=0;i!=gdim;i++)
  {
       for (int j=0;j!=gdim;j++)
       {
	 inb[i*gdim+j].r=grid[i*gdim+j];
	 inb[i*gdim+j].i=0;
	// cout << inc*(i+0.5)-0.5<<" " << inc*(j+0.5)-0.5 << " "   << grid[i*gdim+j] << endl;
       }
  }
  kiss_fftnd(cfg,(kiss_fft_cpx*) inb,  (kiss_fft_cpx*) out2);

  
//   kiss_fftnd_cfg cfg2=kiss_fftnd_alloc ( dims,2, 1,0,0);
//  kiss_fftnd(cfg2, (kiss_fft_cpx*) out,  (kiss_fft_cpx*) inb);
  vector<yvals> fft;
  fft.resize(gdim);
  for (int i=0;i!=gdim;i++)
  {
    fft[i].resize(gdim);
       for (int j=0;j!=gdim;j++)
       {
	 fft[i][j].real(out[i*gdim+j].r+out2[i*gdim+j].r);
	 fft[i][j].imag(out[i*gdim+j].i+out2[i*gdim+j].i);
	 fft[i][j]/=2*gdim*gdim;
//	 
	 //cout <<  out[i*gdim+j].r/1024 << " I "  << inb[i*gdim+j].i/1024 << " X "  << grid[i*gdim+j] << endl;
       }
  }
    free(grid);
free(cfg);

//  memset(out,0,N*N*2*sizeof(double));
  	 
	 /*

  for (int i=0;i!=N;i++)
  {
       for (int j=0;j!=N;j++)
       {
	 if(i<gdim/2 && j<gdim/2) 
	 {
	   
	    out[2*i*N+2*j]=ffty[i][j].real(); 
	    out[2*i*N+2*j+1]=ffty[i][j].imag();
	    
	    out[2*(N-i-1)+2*j]=ffty[gdim-i-1][j].real();
	    out[2*(N-i-1)+2*(N-j-1)]=ffty[gdim-1-i][gdim-j-1].real();
	    out[2*i+2*(N-j-1)]=ffty[i][gdim-j-1].real();

	    out[2*(N-i-1)+2*j+1]=ffty[gdim-i-1][j].imag();
	    out[2*(N-i-1)+2*(N-j-1)+1]=ffty[gdim-i-1][gdim-j-1].imag();
	    out[2*i+2*(N-j-1)+1]=ffty[i][gdim-j-i].imag();

	      
	 } 
       }
  }*/
	
  
  return fft;
}

template <class T> double SinVolumeNorm( vector<int> *coord,double freq, Histogram<T> *h)
{
double vol=1;
// Make dependent on radius!

// 2= 
//cout << " call!" << endl;
for(int i=0; i!=h->grid.limits.size(); i++)
 { int n=(*coord)[i];
  // cout << " mo " << r0 << endl;
   double alpha=h->grid.lengths[i]*1;
//  if(n==0) {vol=(r0+alpha)*(r0+alpha)*PI-(r0*r0)*PI;}
  vol*=sin(alpha*n+0.5*alpha);
  
//   cout <<" n " <<n*h->grid.lengths[i] << " " << vol << endl ;
  }

 return vol;
}


double *PrepareXYZ(vector<tuple> *c,tuple dim,double *ret=NULL)
{
  if(ret==NULL) ret=(double*) malloc(c->size()*3*sizeof(double));
  for(int i=0;i<c->size();i++)
  {
//  cout << "Y " <<  c->at(i).x << " " << c->at(i).y << " " << c->at(i).z << endl;
  ret[3*i]=c->at(i).x/dim.x-0.5;  
  ret[3*i+1]=c->at(i).y/dim.y-0.5;  
  ret[3*i+2]=c->at(i).z;  
  }
  return ret;
}



void fftshift(int xgrid,double *fftarray, vector<yvals> *out)
{

  for(int i=0;i!=xgrid/2;i++)
  {
    for(int j=0;j!=xgrid/2;j++)
    {
      (out->at(i+xgrid/2)[j+xgrid/2]).real(fftarray[2*i*xgrid+2*j]);
       out->at(i+xgrid/2)[j+xgrid/2] .imag( fftarray[2*i*xgrid+2*j+1]);
 
      out->at(i)[j].real(fftarray[2*(i+xgrid/2)*xgrid+2*(j+(int) xgrid/2)]);
      out->at(i)[j].imag(fftarray[2*(i+xgrid/2)*xgrid+2*(j+(int) xgrid/2)+1]);
      
      out->at(i)[xgrid/2+j].real(fftarray[2*(i+xgrid/2)*xgrid+  2*(j)]);
      out->at(i)[xgrid/2+j].imag(fftarray[2*(i+xgrid/2)*xgrid+  2*(j)+1]);

      out->at(xgrid/2+i)[j].real(fftarray[2*i*xgrid  + 2*(j+xgrid/2)] );
      out->at(xgrid/2+i)[j].imag(fftarray[2*i*xgrid  + 2*(j+xgrid/2)+1]);      
      
    }
  }
}

vector<yvals> nfft(double *in,double *mon,double *mon2, int N, int M)
{
 int xgrid=N;
 vector< yvals > ksurf;
  ksurf.resize(N);
  for(int i=0;i!=ksurf.size();i++) ksurf[i].resize(N);
 
//interpol(xgrid,M/2,in, (double*) mon);
//interpol(xgrid,M/2,&in[M/2*3], (double*) mon2);

for(int i=0;i!=(xgrid*xgrid);i++)
 { 
 mon[2*i]+=mon2[2*i];
 mon[2*i]/=2.0;
 mon[2*i+1]+=mon2[2*i+1];
 mon[2*i+1]/=2.0; // Stay real
 }
 fftshift(xgrid, mon, &ksurf);
return ksurf;

}

double *UpdateZ(vector<double> *c,double *ret)
{
  for(int i=0;i<c->size();i++)
  {
  *(ret+3*i+2)=c->at(i);  
  }
  return ret;
}


int main(int argc, char **argv) {

Histogram<double>h;
Grid g;
if(argc !=4) {cout << "xtcfile libraryfile indexfile " << endl; exit(-1); }
XTCFile xtc(argv[1]);
Trajectory traj(&xtc);
SnapShot s=traj.NextFrame();
tuple dim=traj.GetDimensions();
dim.print();

//traj.setCutback(true);
LipidDirectors lipids(&traj);
lipids.LoadLibrary(argv[2]);
lipids.LoadIndices(argv[3]);

lipids.Evaluate(&s);

int j=0;
vector<tuple> *dir;
vector<tuple> *c;
vector<bool> updown;

vector< complex <double> > out;
//g.limits.push_back(sq); //32 12

c=lipids.GetCenters();
dir=lipids.GetDirectors();
AssignLeaflet(updown, dir);
int ucount=0;
for(int i=0;i!=updown.size();i++) if (updown[i]==0) ucount++;
cout << "TOP " << ucount << " BOTTOM " << dir->size()-ucount << endl ; 

out.resize(c->size());
double *in=PrepareXYZ(c,dim);
//cout << "FILED" << endl;
//glacier(256, 400);
int xgrid=sqrt(c->size());
out.resize(xgrid*xgrid);
double xmax=traj.GetDimensions().x;
double qmax=2*PI*xgrid/(xmax/10);
g.limits.push_back(256);
g.width.push_back(make_pair<double,double>(0.0,qmax/2));
h.setGrid(g);

  
//h.normfunc=&(VolumeNorm<double>);
//h.normfunc=&(SinVolumeNorm<double>);
vector<double> dirx,diry;
dirx.resize(c->size());
diry.resize(c->size());
double *mon= (double*) malloc(xgrid*xgrid*sizeof(double)*2);
double *mon2= (double*) malloc(xgrid*xgrid*sizeof(double)*2);
vector<yvals> ksurf;
while(s.Valid())
{
lipids.Evaluate(&s);

//dir=lipids.GetDirectors();
c=lipids.GetCenters();
for(int i=0;i!=c->size();i++) { dirx[i]=dir->at(i).x; diry[i]=dir->at(i).y;}
//cout << "PREP" << endl;
dim=traj.GetDimensions();
xmax=dim.x;
in=PrepareXYZ(c,dim,in);
ksurf = interpol_e(xgrid,c->size(),in,updown);
double area=dim.x*dim.y/100;
// BY SURFACE
double tilesurf=dim.x/10*dim.y/10*((int) sqrt(c->size()/2));
 for(int i=0;i!=ksurf.size()/2;i++)
   for(int j=0;j!=ksurf[i].size()/2;j++)
   { 
     double qfac=2.0*PI/(xmax/10);
    h.AddDatapoint(  (ksurf[i][j]*conj(ksurf[i][j])).real()*area/100, qfac*sqrt((double) (i*i)+ (double) (j*j)));
   //  cout << ksurf[i][j].real()<< " " << ksurf[i][j].imag() << endl;  
   }

   
//exit(0);
/*
in=UpdateZ(&dirx,in);
interpol(xgrid,c->size()/2,in, (double*) &out[0]);
in=UpdateZ(&diry,in);
interpol(xgrid,c->size()/2,in, (double*) &out[0]);*/
//yvals s=xyfft(out);

//cout << "INTERPOL" << endl;
//free(mon);
//for(int i=0;i!=c->size();i++) cout << i/xgrid << " " << i%xgrid << " " << out[i].real() << " " << out[i].imag() << endl;
//exit(0);
s=traj.NextFrame();
cerr << "*" << flush;
//dim=traj.GetDimensions();
}
vector <vector < double> > data;
h.Data(&data);
ofstream of;
of.open((string(argv[1]) +".fasth").c_str());
list <double> ydat;

for(int j=0;j!=data[0].size();j++)
{
 double y=pow(data[0][j],4)*data[1][j];
 of << data[0][j] << " \t " << y << endl; 
 if(data[0][j]<0.5 && y > 0.001) ydat.push_back(y); 
}
  of.close();

ydat.sort();
double avg=0;int n=1,pos=0;
for (list< double >::iterator it = ydat.begin(); it != ydat.end(); it++)
    {
    //  if(pos > ydat.size()/2-4)
      {
      if(n>ydat.size()/2) break;
      avg+=*(it);
	      cout << (*it) << n << endl;

      n++;}pos++;
    }
    cout << "n" << n << endl;
    avg/=(n-1);
cout << "----------------HEIGHT FLUCTUATION--1.0------" << endl;
cout << "Published: C. Allolio., A. Haluts, D. Harries" << endl;
cout << "Chem. Phys. 2018                             " << endl;
cout << "https://doi.org/10.1016/j.chemphys.2018.03.004" << endl;
cout << "----------------RESULTS SUMMARY--------------" << endl;
cout << " FILE: " << argv[1] << endl;
  cout << "Monolayer Bending: kappa " << 1.0/avg/2.0 << " kt" << endl; 

// cout << mon2[2*xgrid/2*xgrid+2*xgrid/2] << mon2[2*xgrid/2*xgrid+2*xgrid/2+1] << endl;
/*
 for(int i=0;i!=c->size();i++) cout << "Xe "  << c->at(i).x << " "  << c->at(i).y << " " << c->at(i).z << endl;
 for(int i=0;i!=ksurf.size();i++)
   for(int j=0;j!=ksurf.size();j++)
   { 
    cout << ksurf[i][j].real()<< " " << ksurf[i][j].imag() << endl;  
   }
*/
free(mon);
free(mon2);

}
