#ifndef HAVE_CONFIG
#include "config.hpp"
#endif
typedef vector<int> atomlist;
#include <complex.h>
#include "kiss_fft130/tools/kiss_fftnd.h"
#include "trformats.hpp"
#include "lipidrecsum.hpp"
#include "histogram.hpp"
#include "funcfit.hpp"
#define NFFT_PRECISION_DOUBLE


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

void bin_surface(double *grid, int gdim, double*in, vector<bool> updown, bool up=1)
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
  //  if(y >=1 || x >=1 || y<=0 || x <=0)
  //  cout << x << " " << y << " " << endl; 
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
    else {/*cerr << "WARNING! No data"  << endl*/; missing.push_back(i);} 
  }
  for(int i=0;i!=missing.size();i++)
  {
    int x=missing[i]%gdim;
    int y=missing[i]/gdim;
   // cout << x << " " << y << " " << gdim << endl;
    double adjy=grid[((y+1)%gdim)*gdim+(x)%gdim];
    double adjx=grid[((y)%gdim)*gdim+(x+1)%gdim];
    if (x==0) x=gdim;
    if (y==0) y=gdim;
   // cout << x << " " << y << " " << gdim << endl;
    double adjmx=grid[y*gdim+(x-1)];
    double adjmy=grid[(y-1)*gdim+(x)];
    grid[missing[i]]=0.25*(adjx+adjy+adjmx+adjmy);
  }
  // Stretch !
  free(weights);

}

void fft_surface(double *grid, int gdim, kiss_fft_cpx* out)
{
    int dims[2]; dims[0]=gdim,dims[1]=gdim;

  kiss_fft_cpx inb[gdim*gdim]; 
   for (int i=0;i!=gdim;i++)
  {
       for (int j=0;j!=gdim;j++)
       {
	 inb[i*gdim+j].r=grid[i*gdim+j];
	 inb[i*gdim+j].i=0;
       }
  }
  kiss_fftnd_cfg  cfg=kiss_fftnd_alloc ( dims,2, 0,0,0);
  kiss_fftnd(cfg,(kiss_fft_cpx*) inb,  (kiss_fft_cpx*) out);
 // free(cfg);
}

vector<yvals> interpol_e(int N,int M, double *in, vector<bool>& updown)
{
  //N output grid
  //M number of points
  int gdim=sqrt(M/2);
 // if(!isPower2(gdim)) gdim=nextp2(gdim/2);
  double *grid=(double*) malloc(gdim*gdim*sizeof(double));
  kiss_fft_cpx out[gdim*gdim];
  kiss_fft_cpx out2[gdim*gdim];
  bin_surface(grid,gdim,in,updown,1); 
  fft_surface(grid,gdim,out);
  bin_surface(grid,gdim,in,updown,0); 
  fft_surface(grid,gdim,out2);
 // kiss_fftnd_cfg cfg2=kiss_fftnd_alloc ( dims,2, 1,0,0);
//  kiss_fftnd(cfg2, (kiss_fft_cpx*) out,  (kiss_fft_cpx*) inb);
  vector<yvals> fft;
  fft.resize(gdim);
  for (int i=0;i!=gdim;i++)
  {
    fft[i].resize(gdim);
       for (int j=0;j!=gdim;j++)
       {
	 fft[i][j].real(out[i*gdim+j].r-out2[i*gdim+j].r);
	 fft[i][j].imag(out[i*gdim+j].i-out2[i*gdim+j].i);
	 fft[i][j]/=2*gdim*gdim;
//	 
	 //cout <<  out[i*gdim+j].r/1024 << " I "  << inb[i*gdim+j].i/1024 << " X "  << grid[i*gdim+j] << endl;
       }
  }
    free(grid);
  return fft;
}
void prepare_npar(int gdim, vector<yvals> &fftx, vector<yvals> &ffty,vector<yvals> &fftxy)
{
   for (int i=0;i!=gdim;i++)
  {
    for(int j=0;j!=gdim;j++)
    {
      fftxy[i][j]=(fftx[i][j]*(double) j+ffty[i][j]*(double) i)/sqrt((double) (i*i)+ (double) (j*j));
    }
  }
   fftxy[0][0]=0;
}

void prepare_nperp(int gdim, vector<yvals> &fftx, vector<yvals> &ffty,vector<yvals> &fftxy)
{
   for (int i=0;i!=gdim;i++)
  {
    for(int j=0;j!=gdim;j++)
    {
      fftxy[i][j]=(ffty[i][j]*(double) j-fftx[i][j]*(double) i)/sqrt((double) (i*i)+ (double) (j*j));
    }
  }
   fftxy[0][0]=0;
}


/*
void prepare_nxy(int gdim, vector<yvals> &fftx, vector<yvals> &ffty,vector<yvals> &fftxy)
{
   for (int i=0;i!=gdim;i++)
  {
    for(int j=0;j!=gdim;j++)
    {
      fftxy[i][j]=(fftx[j][i]*(double) i+ffty[j][i]*(double) j)/sqrt((double) (i*i)+ (double) (j*j));
    }
  }
   fftxy[0][0]=0;
 
}*/

vector<yvals> interpol_nxy(int N,int M, double *inx,double *iny,vector<bool> &updown, vector<yvals> &nperp)
{
  //N output grid
  //M number of points
  int gdim=sqrt(M/2);
//  cout << gdim << endl ;
  //if(!isPower2(gdim)) gdim=nextp2(gdim/2);
  double *gridx=(double*) malloc(gdim*gdim*sizeof(double));
  double *gridy=(double*) malloc(gdim*gdim*sizeof(double));
  double *gridxy=(double*) malloc(gdim*gdim*sizeof(double));
    int dims[2]; dims[0]=gdim,dims[1]=gdim;
    
  kiss_fft_cpx out[gdim*gdim];
  kiss_fft_cpx out2[gdim*gdim];
  bin_surface(gridx,gdim,inx,updown,1);
 fft_surface(gridx,gdim,out);
  bin_surface(gridx,gdim,inx,updown,0); 
  fft_surface(gridx,gdim,out2);
   vector<yvals> fftx;
  fftx.resize(gdim);
  for (int i=0;i!=gdim;i++)
  {
    fftx[i].resize(gdim);
       for (int j=0;j!=gdim;j++)
       {
	 fftx[i][j].real(out[i*gdim+j].r-out2[i*gdim+j].r);
	 fftx[i][j].imag(out[i*gdim+j].i-out2[i*gdim+j].i);
	 fftx[i][j]/=gdim*gdim*2;//2
//	 
	 //cout <<  out[i*gdim+j].r/1024 << " I "  << inb[i*gdim+j].i/1024 << " X "  << grid[i*gdim+j] << endl;
       }
  }
  vector<yvals> ffty;
  ffty.resize(gdim);
  kiss_fft_cpx out3[gdim*gdim];
  kiss_fft_cpx out4[gdim*gdim];

  bin_surface(gridy,gdim,iny,updown,1);
  fft_surface(gridy,gdim,out3);
  bin_surface(gridy,gdim,iny,updown,0); 
  fft_surface(gridy,gdim,out4);
  
  for (int i=0;i!=gdim;i++)
  {
    ffty[i].resize(gdim);
       for (int j=0;j!=gdim;j++)
       {
	 ffty[i][j].real(out3[i*gdim+j].r-out4[i*gdim+j].r);
	 ffty[i][j].imag(out3[i*gdim+j].i-out4[i*gdim+j].i);
	 ffty[i][j]/=gdim*gdim*2;
//	 
	 //cout <<  out[i*gdim+j].r/1024 << " I "  << inb[i*gdim+j].i/1024 << " X "  << grid[i*gdim+j] << endl;
       }
  }
     vector<yvals> fftxy=fftx;
 if(nperp.size()!=fftxy.size()) nperp=fftxy;
  prepare_npar(gdim,fftx,ffty,fftxy);
  prepare_nperp(gdim,fftx,ffty,nperp);

 /*
  cout << " M" << endl;
    for(int j=0;j!=M;j++)
    cout << inx[3*j] << " " << inx[3*j+1] << " " <<  inx[3*j+2] << endl;
  */
  
//   kiss_fftnd_cfg cfg2=kiss_fftnd_alloc ( dims,2, 1,0,0);
//  kiss_fftnd(cfg2, (kiss_fft_cpx*) out,  (kiss_fft_cpx*) inb);
  if(gridx!=NULL)  free(gridx);
  if(gridy!=NULL)     free(gridy);

  return fftxy;
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


double *PrepareXYZ(vector<triple> *c,triple dim,double *ret=NULL)
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
/*
vector<yvals> nfft(double *in,double *mon,double *mon2, int N, int M)
{
 int xgrid=N;
 vector< yvals > ksurf;
  ksurf.resize(N);
  for(int i=0;i!=ksurf.size();i++) ksurf[i].resize(N);
 
interpol(xgrid,M/2,in, (double*) mon);
interpol(xgrid,M/2,&in[M/2*3], (double*) mon2);

for(int i=0;i!=(xgrid*xgrid);i++)
 { 
 mon[2*i]+=mon2[2*i];
 mon[2*i]/=2.0;
 mon[2*i+1]+=mon2[2*i+1];
 mon[2*i+1]/=2.0; // Stay real
 }
 fftshift(xgrid, mon, &ksurf);
return ksurf;

}*/

double *UpdateZ(vector<triple> *c ,double *ret, int axis)
{
  for(int i=0;i<c->size();i++)
  {
  *(ret+3*i+2)=c->at(i).p[axis];  
  }
  return ret;
}

void AssignLeaflet(vector<bool> &updown, vector <triple> *dir)
{
  triple plane;
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


int main(int argc, char **argv) {

vector <bool> updown;
Histogram<double>h;
Histogram<double>htilt;

Grid g;
g.limits.push_back(1024); //32 12
if(argc !=4) {cout << "xtcfile libraryfile indexfile " << endl; exit(-1); }
XTCFile xtc(argv[1]);
Trajectory traj(&xtc);
SnapShot s=traj.NextFrame();
triple dim=traj.GetDimensions();
dim.print();
cout << "-----------------PRE-RELEASE-----------------" << endl;
cout << "FOURIER_SPACE TILT FLUCTUATION               " << endl;
cout << "(c) C.Allolio " << __DATE__<<  " "<< __TIME__ << endl;
cout << "---------------------------------------------" << endl;

//traj.setCutback(true);
LipidDirectors lipids(&traj);
lipids.LoadLibrary(argv[2]);
lipids.LoadIndices(argv[3]);

lipids.Evaluate(&s);

int j=0;
vector<triple> *dir;
vector<triple> *c;
vector< complex <double> > out;
c=lipids.GetCenters();
out.resize(c->size());
dir=lipids.GetDirectors();
AssignLeaflet(updown, dir);
int ucount=0;
for(int i=0;i!=updown.size();i++) if (updown[i]==0) ucount++;
cout << "TOP " << ucount << " BOTTOM " << dir->size()-ucount << endl ;
double *inx=PrepareXYZ(c,dim);
double *iny=PrepareXYZ(c,dim);

//cout << "FILED" << endl;
//glacier(256, 400);
int xgrid=256;
out.resize(xgrid*xgrid);
double xmax=(traj.GetDimensions()).x;
double qmax=2*PI*xgrid/(xmax/10) ;///;
g.width.push_back(make_pair<double,double>(0.0,qmax/2));
h.setGrid(g);
htilt.setGrid(g);
dim=traj.GetDimensions();
  
//h.normfunc=&(VolumeNorm<double>);
//h.normfunc=&(SinVolumeNorm<double>);
double *mon= (double*) malloc(xgrid*xgrid*sizeof(double)*2);
double *mon2= (double*) malloc(xgrid*xgrid*sizeof(double)*2);
//vector<yvals> nxsurf;
vector<yvals> nxysurf;
vector<yvals> nvertsurf;

dim=traj.GetDimensions();
vector<double> vol;

while(s.Valid())
{
lipids.Evaluate(&s);

dir=lipids.GetDirectors();
c=lipids.GetCenters();
//for(int i=0;i!=c->size();i++) { dir->at(i).x=static_cast <float> (rand()) / static_cast <float> (RAND_MAX); dir->at(i).y=static_cast <float> (rand()) / static_cast <float> (RAND_MAX);}
//cout << "PREP" << endl;
xmax=dim.x;
inx=PrepareXYZ(c,dim,inx);
inx=UpdateZ(dir,inx,0);
iny=PrepareXYZ(c,dim,iny);
iny=UpdateZ(dir,iny,1);
nxysurf = interpol_nxy(xgrid,c->size(),inx,iny,updown,nvertsurf);
//nxysurf = interpol_e(xgrid,c->size(),iny);
//for(int i=0;i!=c->size();i++) { cout << dir->at(i).x - inx[i*3+2] << " " << dir->at(i).y - iny[i*3+2] << endl;}
//cout << "PREP" << endl;

double area=dim.x*dim.y/100;
vol.push_back(area);
// BY SURFACE
// nmtile =
double tilesurf=dim.x/10*dim.y/10*((int) sqrt(c->size()/2));
 for(int i=0;i!=nxysurf.size()/2;i++)
   for(int j=0;j!=nxysurf[i].size()/2;j++)
   { 
//     cout << (nxysurf[i][j]*conj(nxysurf[i][j])).real() << " "<< i << " " <<  j << endl;
     double qfac=2.0*PI/(xmax/10);
    h.AddDatapoint( (nxysurf[i][j]*conj(nxysurf[i][j])).real(), sqrt((double) (i*i)+ (double) (j*j)));
    htilt.AddDatapoint( (nvertsurf[i][j]*conj(nvertsurf[i][j])).real(), sqrt((double) (i*i)+ (double) (j*j)));
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
//cout << dim.x;
//cerr << "*" << flush;
//dim=traj.GetDimensions();
dim=traj.GetDimensions();

}
vector <vector < double> > data;
vector <vector < double> > tdata;

h.Data(&data);
htilt.Data(&tdata);

ofstream of,tilt;
double area=0;
for(int i=0;i!=vol.size();i++) area+=vol[i];
area/=vol.size();
double qfac=2.0*PI/(sqrt(area)/10);
of.open((string(argv[1]) +".brown").c_str());
tilt.open((string(argv[1]) +".btilt").c_str());

vector <double> xdat;
vector <double> ydat;
vector <double> txdat;
vector <double> tydat;
txdat.clear();tydat.clear();
for(int j=0;j!=data[0].size();j++)
{
 double x= data[0][j]*qfac/10; double y=data[0][j]*qfac*data[0][j]*qfac*data[1][j]*area/100;
 of << x << " \t " << y << endl; 
 if(x < 1.0 && y > 0 ) {xdat.push_back(x);ydat.push_back(y);} 
}
  of.close();
  for(int j=0;j!=tdata[0].size();j++)
{
 double x= tdata[0][j]*qfac/10; double y=tdata[1][j]*area;//*area/100;//*area/100;
 if(y>0) tilt << x << " \t " << y << endl; 
 if(x < 2.0 && y > 0 ) {txdat.push_back(x);tydat.push_back(1/y);} 

// if(x < 1.0 && y > 0 ) {xdat.push_back(x);ydat.push_back(y);} 
}
 tilt.close();

    MatrixXd xvals(xdat.size(),1);
    MatrixXd yvals(ydat.size(),1);
 
      FitConst cons;
        VectorXd param(2);
	param(0)= 0.02;
    cons.setparam(param);
      
  for(int i=0;i!=xdat.size();i++)
  {
    xvals(i,0)=xdat[i];
    yvals(i,0)=ydat[i];
  }
    FuncFunctord dafunc(xvals,yvals,cons);
    Eigen::LevenbergMarquardt<FuncFunctord> lm2(dafunc);
   Eigen::LevenbergMarquardtSpace::Status status = lm2.minimize(param);
cout << "-------------Director Fluctuation v.1----------" << endl;
cout << "Method: Max C. Watson, Erik G. Brandt," << endl; 
cout <<" Paul M. Welch, and Frank L. H. Brown" << endl;
cout <<" Phys. Rev. Lett. 109, 028102"<< endl;
cout <<" -----------------------------------------------" << endl;
cout << "Implemented: C. Allolio., A. Haluts, D. Harries" << endl;
cout << "Chem. Phys. 2018                             " << endl;
cout << "https://doi.org/10.1016/j.chemphys.2018.03.004" << endl;
   cout << "----------------RESULTS SUMMARY--------------" << endl;
cout << " FILE: " << argv[1] << endl;
  cout << "Monolayer Bending: kappa " << 1.0/param(0)/2.0 << " kt" << endl;
  	
	     FitParabola parb;
    
    MatrixXd txvals(txdat.size(),1);
    MatrixXd tyvals(tydat.size(),1);
      for(int i=0;i!=txdat.size();i++)
  {
    txvals(i,0)=txdat[i];
    tyvals(i,0)=tydat[i];
  }

     param(0)= 3.00;
      param(1)= tyvals(0);
     parb.setparam(param);
 //    cout << txvals << endl;
 //    cout << "TYVALS" << endl;
  //   cout << tyvals << endl;
  FuncFunctord daf2(txvals,tyvals,parb);
    Eigen::LevenbergMarquardt<FuncFunctord> lm3(daf2);
   Eigen::LevenbergMarquardtSpace::Status status2 = lm3.minimize(param);
   cout << "Bilayer Tilt: kappa . . " << param(1) << " kt/nm^2 Offset: " << 0 << endl; 
   cout << "Bilayer Twist: kappa " << param(0) << " kt/nm^2 Offset: " << 0 << endl; 
  cout << "----------------------------------------------" << endl;

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
free(inx);free(iny);
}
