#ifndef HAVE_CONFIG
#include "config.hpp"
#endif
typedef vector<int> atomlist;

#include "trformats.hpp"
#include "lipidrecsum.hpp"
#include "histogram.hpp"
#include "diffgeom.hpp"
#include "funcfit.hpp"


template <class T> double SinVolumeNorm( vector<int> *coord,double freq, Histogram<T> *h)
{
double vol=1;
for(int i=0; i!=h->grid.limits.size(); i++)
 { int n=(*coord)[i];
   double alpha=h->grid.lengths[i]*1;
  vol*=sin(alpha*n+0.5*alpha);
  }

 return vol;

}


double EstimateScaling(Surface &s)
{
  Geometry g(&s);
  Matrix2d metric;
  Vector2d pos(0,0);
  g.GetI(pos,metric);
  return 1.0/sqrt(metric.determinant());
}

double SurfaceIntegral(Surface &s,double umin, double umax, double vmin,double vmax, int points)
{
  double A=0;
  Matrix2d metric;
  Geometry g(&s);
  double ucurr=umin;
  double vcurr=vmin;
  double ustep=(umax-umin)/points;
  double vstep=(vmax-vmin)/points;
  int i=0;
  for(i=0;i!=points;i++)
   for(int j=0;j!=points;j++)
    {   ucurr+=ustep;
      vcurr=vstep*i;
       Vector2d pos(ucurr,vcurr);
      g.GetI(pos,metric);
      double sqrtg=sqrt(metric.determinant());
      A+=sqrtg*ustep*vstep;
    }
  
  return A;
}

int LipidsInArea(vector <tuple> &lipids,double xmin,double xmax, double ymin,double ymax)
{
  int counter=0;
  for(int i=0;i!=lipids.size();i++)
  {
    if(lipids[i].x > xmin && lipids[i].x < xmax && lipids[i].y > ymin && lipids[i].y < ymax) counter++;
  }
  return counter;
}

bool FitPolyLM2d(vector<tuple> &points , P2Surface &p2s )
{
  if(points.size()<6) return false;
  VectorXd c=VectorXd::Zero(6);

  MatrixXd xsvals(points.size(),2);
  MatrixXd ysvals(points.size(),3);
  for(int i=0;i!=points.size();i++)
  {
    Vector3d point(points[i].x,points[i].y,points[i].z);
    xsvals.row(i)=Vector2d(point(0),point(1)); 
    ysvals.row(i)=point; 
  }
      P2FSurface surf2;

   {  Matrix3d b; b.setIdentity(); surf2.setup(&c,&b);}
    FuncFunctord dasurf(xsvals,ysvals,surf2);
    Eigen::LevenbergMarquardt<FuncFunctord> lm2(dasurf);
    Eigen::LevenbergMarquardtSpace::Status status;
    status=lm2.minimize(c);
    if(status==1) return false;
  {  Matrix3d b; b.setIdentity(); p2s.setup(&c,&b);}
  return true;
  
}

bool FitPlaneLM2d(vector<tuple> &points , Plane &p2s )
{
  if(points.size()<4) return false;
  VectorXd c=VectorXd::Zero(3);

  MatrixXd xsvals(points.size(),2);
  MatrixXd ysvals(points.size(),3);
  for(int i=0;i!=points.size();i++)
  {
    Vector3d point(points[i].x,points[i].y,points[i].z);
    xsvals.row(i)=Vector2d(point(0),point(1)); 
    ysvals.row(i)=point; 
  }
      P2FPlane surf2;

   {  Matrix3d b; b.setIdentity(); surf2.setup(&c,&b);}
    FuncFunctord dasurf(xsvals,ysvals,surf2);
    Eigen::LevenbergMarquardt<FuncFunctord> lm2(dasurf);
    Eigen::LevenbergMarquardtSpace::Status status;
    status=lm2.minimize(c);
    if(status==1) return false;
  {  Matrix3d b; b.setIdentity(); p2s.setup(&c,&b);}
  return true;
  
}



bool FitPolyLM2dIter(vector<tuple> &points , P2Surface &p2s )
{
  VectorXd c(6);
  double coff=0.5;
  Matrix3d b; b.setIdentity();
  list<Vector3d> local;
  local.clear();
  bool converged=false;
  for(int i=0;i!=points.size();i++)
  { local.push_back(Vector3d((double*) &points[i].p));}
  while(!converged)
  {
  if(local.size()<10) return false;
  MatrixXd xsvals(local.size(),2);
  MatrixXd ysvals(local.size(),3);
  int i=0;
  for (list <Vector3d> ::iterator it = local.begin(); it != local.end(); it++)
  {
    xsvals.row(i)=Vector2d((*it)(0),(*it)(1)); 
    ysvals.row(i)=*it; i++;
  }
      P2FSurface surf2;

   {   surf2.setup(&c,&b);}
    FuncFunctord dasurf(xsvals,ysvals,surf2);
    Eigen::LevenbergMarquardt<FuncFunctord> lm2(dasurf);
    Eigen::LevenbergMarquardtSpace::Status status;
    status=lm2.minimize(c);
    int j=0;
    converged=true;
    if(status==1) {return false;coff=3;converged=false;}
   for (list <Vector3d> ::iterator it = local.begin(); it != local.end(); it++)
  {
    if(surf2.eval(xsvals.row(j))(2)-(ysvals.row(j))(2) > coff) {local.erase(it); it--;converged=false;}
    j++;
  }
   
  }
  {  p2s.setup(&c,&b);}
  return true;
  
}





bool FitGParb(vector< vector <double > > &data, vector <double> &result , bool shift=false,double nfac=1.0 )
{
  
    MatrixXd xvals(data[0].size(),1);
    
    MatrixXd yvals(data[0].size(),1);
  
  for(int i=0;i!=data[0].size();i++)
  {
    xvals(i,0)=data[0][i];
    yvals(i,0)=-log((double)(data[1][i]*nfac));
  }
      FitParabola parb;
      FitShParabola parb2;
        VectorXd param(3);
	param(0)= 5.; param(1)= -10.; param(2)=0;
  
    parb.setparam(param);
    parb2.setparam(param);
    FuncFunctord dafunc(xvals,yvals,parb);
        FuncFunctord dafunc2(xvals,yvals,parb2);

    Eigen::LevenbergMarquardt<FuncFunctord> lm1(dafunc);
    Eigen::LevenbergMarquardt<FuncFunctord> lm2(dafunc2);
   if(shift) lm2.minimize(param) ;
   else Eigen::LevenbergMarquardtSpace::Status status = lm1.minimize(param);
    result.resize(3);
    result[0]=param[0];
    result[1]=param[1];
    result[2]=param[2];

  return true;
  
}


tuple Tilt(Surface &s,tuple &n, double x=0, double y=0)
{
  Vector3d Nx = s.unormv(Vector2d(x,y));
  tuple N(Nx.data());
  if(N*n < 0) n=n*-1; 
  return n-N; ///(N*n)
}

bool Align(vector<tuple> &director, Matrix3d &bas)
{
  
  for(int i=0;i!=director.size();i++)
  { Vector3d e=bas.transpose()*Vector3d(director[i].x,director[i].y,director[i].z);
    director[i].x=e(0);director[i].y=e(1);director[i].z=e(2);
  }
  return true;
}


bool FitTilt(vector<tuple> &pos,vector<tuple> &director, Surface &s,Plane &tus, Plane &tvs)
{
 Vector2d upos(pos[0].x,pos[0].y);
 tuple N(s.unormv(upos).data());
     Vector3d sdu=s.du(upos);
 Vector3d sdv=s.dv(upos);


 tuple eu(sdu.data());
 tuple ev(sdv.data());

 
  if(director[0]*N<0 ){s.flipnorm(); N=N*-1;}
 
 vector<tuple> tu(5);
 vector<tuple> tv(5);
 for(int i=0;i!=5;i++)
 {
   upos(0)=(pos[i].x);
   upos(1)=(pos[i].y);
  N=(s.unormv(upos).data()); 
  tuple &n=director[i];
  tuple tilt=n-N;//(n/(N*n))-N
  tu[i]=pos[i];
  tu[i].z=eu*tilt;
  tv[i]=pos[i];
  tv[i].z=ev*tilt;
 }
 FitPlaneLM2d(tu,tus);
FitPlaneLM2d(tv,tvs);
}

Matrix2d TiltTensor(const Vector2d upos,Surface &tus, Surface &tvs, Geometry &g)
{
  Matrix2d tmat;
  double ta=tus.eval(upos)[2];
  double tb=tvs.eval(upos)[2];
  
  tmat(0,0)=tus.du(upos)[2]+g.Gamma(upos,0,0,0)*ta+g.Gamma(upos,0,0,1)*tb;
  tmat(0,1)=tus.dv(upos)[2]+g.Gamma(upos,0,1,0)*ta+g.Gamma(upos,0,1,1)*tb;
  tmat(1,0)=tvs.du(upos)[2]+g.Gamma(upos,1,0,0)*ta+g.Gamma(upos,1,0,1)*tb;
  tmat(1,1)=tvs.dv(upos)[2]+g.Gamma(upos,1,1,0)*ta+g.Gamma(upos,1,1,1)*tb;
  return tmat;
}

Matrix2d Btensor(const Vector2d &X, Geometry &g)
{
  Matrix2d bmat;
  g.GetShape(X,bmat);
  return bmat;
}


int main(int argc, char **argv) {

  
Histogram<double>h;
Histogram<double>ht;

Grid g,g1;
g.limits.push_back(300); //32 12
g.width.push_back(make_pair<double,double>(-0.2,0.2));
g1.limits.push_back(70); //32 12
g1.width.push_back(make_pair<double,double>(0,1.2));

h.setGrid(g);
ht.setGrid(g1);
bool witchcraft=false;
bool cutback=false;
if(argc <4) {cout << "xtcfile libraryfile indexfile " << endl; exit(-1); }
if(argc >4 && *argv[4]=='w' ) witchcraft=true;
if(argc >4 && *argv[4]=='c' ) cutback=true;

XTCFile xtc(argv[1]);
Trajectory traj(&xtc);
if(!cutback) traj.FullPbc(true);
traj.setCutback(cutback);
LipidDirectors lipids(&traj);
lipids.LoadLibrary(argv[2]);
lipids.LoadIndices(argv[3]);
SnapShot s=traj.Current();
int j=0;
tuple conc;
conc.x=0.0;
conc.y=0.0;
conc.z=1.0;
h.normfunc=&(VolumeNorm<double>);
ht.normfunc=&(SinVolumeNorm<double>);
vector<double> vol;
vector<double> avbend;
vector<double> avk;
vector<tuple> *dir;
vector<tuple> *c;
Vector2d cent;
cent=cent*0;
while(s.Valid())
{
//BUILIDING LIPIDS
lipids.Evaluate(&s);
dir=lipids.GetDirectors();
c=lipids.GetCenters();
//GETTING ENV
double AVGA=0;
double AVGPL=0;
double AVGBEND=0;
double AVGK=0;
int contrib=0;
srand48((long int) dir);
for(int i=0;i<c->size();i++)
{
 vector<tuple> lpos(20); // 25
 vector<tuple> ldir(6);
 atomlist indices;
 if(lipids.GetLipidEnv(i,lpos,indices,25,false))
 {
 P2Surface s;
 Plane tvs,tus;
 Matrix3d basis; 
 PrincipalAxes(lpos,basis);
 Align(lpos,basis);
 if(FitPolyLM2d(lpos,s))
 {
 double exc=drand48()*1.0-0.5;
 double cut=3.66; //3.7
 cut*=pow(EstimateScaling(s),2);
 double A=SurfaceIntegral(s,-cut+exc,cut+exc,-cut+exc,cut+exc,20);
 double PL=LipidsInArea(lpos,-cut+exc,cut+exc,-cut+exc,cut+exc);
 AVGA+=A;
 AVGPL+=PL;
   contrib++;

 Geometry go(&s);
 for(int j=0;j!=ldir.size();j++) {  Vector3d d((double*) &dir->at(indices[j]).p);d=basis.transpose()*d; ldir[j]=tuple(d.data());  }
 FitTilt(lpos,ldir,s,tus,tvs);
 Matrix2d tab=TiltTensor(cent,tus,tvs,go); 
 Matrix2d bab=Btensor(cent,go);
 double t=(Tilt(s,ldir[0])).abs();
 double Splay=tab.trace();
 AVGBEND+=bab.trace();
 AVGK+=bab.determinant();
 if(abs(Splay) > 5 || Splay !=Splay ) 
 {
if(argc >4 && argv[4]=="t")  for(int j=0;j!=6;j++) cerr << ldir[j]*ldir[0] << " D " << lpos[j].x << " "<< lpos[j].y << " " << lpos[j].z << " I "  << indices[j] << endl; 
  
 }
 double KG=bab.determinant(); 
 double angle=acos(dir->at(i)*conc);
  	  if(angle>PI/2) angle=PI-angle;
 //t=abs(angle);
 h.AddDatapoint(1,Splay);
 ht.AddDatapoint(1,t);//t

  

if(argc > 4) 
{   
      Vector3d Nx = s.unormv(Vector2d(0.0,0.0));
      Nx=basis*Nx;
    cout << "KG : " << bab.determinant() << " B: " << go.GetH(cent) << " T: " << tab.trace() << " DL: "<< sqrt(Splay*Splay/4-KG) << " Frame " << avk.size() << " p  " << c->at(i).x << " " << c->at(i).y << " " <<c->at(i).z << " n " << dir->at(i).x << " " << dir->at(i).y << " " << dir->at(i).z << " N " << Nx(0) << " " << Nx(1) << " " << Nx(2) << endl;}
}}}
double AVGAPL=AVGA/AVGPL;

AVGBEND/=contrib;
AVGK/=contrib;
Matrix3 box;
box=traj.GetBoxMat();
tuple dim;
dim.x=box.row[0].x;
dim.y=box.row[1].y;
dim.z=box.row[2].z;
if(!witchcraft) AVGAPL= dim.x*dim.y/c->size()*2.0;
if(AVGAPL==AVGAPL) { vol.push_back(AVGAPL); avbend.push_back(AVGBEND); avk.push_back(AVGK);}
s=traj.NextFrame();
}
double avg=0;
double stdev=0;
double ba=0;
vector<double> blocks;

for(int i=0;i<vol.size();i++)
 {
   avg+=vol[i];
   ba+=vol[i];
  stdev+=vol[i]*vol[i];
  if(i%(vol.size()/10) == 0 && i>0) {blocks.push_back((ba/(vol.size()/10))); ba=0;} 
 }
 double bend=0,gauss=0;
 for(int i=0;i<avbend.size();i++)
 {
   gauss+=avk[i];
   bend+=avbend[i];
 }
 bend/=avbend.size();
 gauss/=avk.size();
 ba=0;double ba2=0;
 for(int i=0;i!=blocks.size();i++)
 {
  ba+=blocks[i];
  ba2+=blocks[i]*blocks[i];
 }
 ba/=blocks.size();
 ba2/=blocks.size();
 ba2=sqrt(ba2-ba*ba);

 avg/=vol.size();
stdev/=vol.size();
stdev=sqrt(stdev-avg*avg);
cout << "-----------------RESIS-1.0-RELEASE------------" << endl;
cout << "(c) C.Allolio " << __DATE__<<  " "<< __TIME__ << endl;
cout << "---------------------------------------------" << endl;
cout << "Published: C. Allolio., A. Haluts, D. Harries" << endl;
cout << "Chem. Phys. 2018                             " << endl;
cout << "https://doi.org/10.1016/j.chemphys.2018.03.004" << endl;
cout << "---------------------------------------------" << endl;
cout << "----------------RESULTS SUMMARY--------------" << endl;
cout << " FILE: " << argv[1] << endl;
double lipnum=lipids.GetCenters()->size()/2.0;
if(witchcraft) cout << "MONTE CARLO APL ESTIMATE" << endl;
cout << " Area AVG [A^2]: " << avg*lipnum << " STDEV: " << stdev*lipnum << " Per Lipid [bilayer]: " << avg << " +- STDERR(10) " << ba2/sqrt(blocks.size())  << " STDEV " << stdev << endl;  
cout << " Mean Curvature /  Lipid [nm^-1]: " << bend*10 << endl;
cout << " Gausian Curvature / Lipid [nm^-2]: " << gauss*100 << endl;
ofstream of,of2;
of.open((string(argv[1]) +".spl").c_str());
of << "# Area AVG [A^2]: " << avg*lipnum << " STDEV: " << stdev*lipnum << " Per Lipid [bilayer]: " << avg << " +- STDERR(10) " << ba2/sqrt(blocks.size())  << " STDEV " << stdev << endl;    
h.DropData(of);
vector<vector <double> > data;
h.Data(&data);
vector<double> res;

FitGParb(data,res,true);
double shift=res[2];
vector<vector <double> > tmp;
tmp.resize(2);
for(int i=0;i!=data[0].size();i++)
{
if(data[0][i]>-0.05 && data[0][i]<0.05) {tmp[0].push_back(data[0][i]);tmp[1].push_back(data[1][i]);}
}
FitGParb(tmp,res);
cout << "Monolayer Bending: kappa " << res[0]*2/(avg) << " kt Offset [nm^-1]: " << shift*10 << endl; 
of.close();
of2.open((string(argv[1]) +".itilt").c_str());
of2 << "# Area AVG [A^2]: " << avg*lipnum << " STDEV: " << stdev << " Per Lipid [bilayer]: " << avg*lipnum << " +- STDERR(10) " << ba2/sqrt(blocks.size())  << " STDEV " << stdev << endl;    
ht.DropData(of2);
ht.Data(&data);
tmp=data;
double std=0,avt2=0,avg2=0,norm=0,amax=0, xmax=0;
for(int i=0;i!=data[0].size();i++)
{
tmp[1][i]*=sin(tmp[0][i]); 
norm+=tmp[1][i];
if(amax< tmp[1][i] ){ amax= tmp[1][i];xmax=tmp[0][i];}; 
}


for(int i=0;i!=data[0].size();i++)
{
tmp[1][i]/=norm;
avt2+=tmp[1][i]*tmp[0][i]*tmp[0][i];
avg2+=tmp[1][i]*tmp[0][i];
}
std=sqrt(avt2-avg2*avg2);
cout << std << " AVG " <<avg2 << endl ;
tmp[0].resize(0);
tmp[1].resize(0);
for(int i=0;i!=data[0].size();i++) if(data[0][i] > (xmax-std) &&  data[0][i] < (xmax+std) ) {tmp[0].push_back(data[0][i]);tmp[1].push_back(data[1][i]);}
FitGParb(tmp,res);
cout << "Bilayer Tilt: kappa " << res[0]*2/(avg/100) << " kt/nm^2 Offset: " << res[2] << endl; 
cout << "-------------------END SUMMARY--------------" << endl;

of2.close();

}
