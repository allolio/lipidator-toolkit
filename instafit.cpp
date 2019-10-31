#ifndef HAVE_CONFIG
#include "config.hpp"
#endif
typedef vector<int> atomlist;

#include "trformats.hpp"
#include "lipidrecsum.hpp"
#include "histogram.hpp"
#include "diffgeom.hpp"
#include "funcfit.hpp"

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
      //sqrtg;
      A+=sqrtg*ustep*vstep;
  //    cout << A << "X "<< endl;
    }
  
  return A;
}

int LipidsInArea(vector <triple> &lipids,double xmin,double xmax, double ymin,double ymax)
{
  int counter=0;
  for(int i=0;i!=lipids.size();i++)
  {
    if(lipids[i].x > xmin && lipids[i].x < xmax && lipids[i].y > ymin && lipids[i].y < ymax) counter++;
  }
  return counter;
}

bool FitPolyLM2d(vector<triple> &points , P2Surface &p2s )
{
  if(points.size()<6) return false;
  // x^2 y^2 xy x y c 
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
  //   else( surf2.setup(&c,bas));

//    Matrix3d bas;
  //  bas.setIdentity();
    FuncFunctord dasurf(xsvals,ysvals,surf2);
    Eigen::LevenbergMarquardt<FuncFunctord> lm2(dasurf);
    Eigen::LevenbergMarquardtSpace::Status status;
    status=lm2.minimize(c);
    if(status==1) return false;
  {  Matrix3d b; b.setIdentity(); p2s.setup(&c,&b);}
//else{    p2s.setup(&c,bas);}
  // coefficients=x;
  return true;
  
}

bool FitPlaneLM2d(vector<triple> &points , Plane &p2s )
{
  if(points.size()<4) return false;
  // x^2 y^2 xy x y c 
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
  //   else( surf2.setup(&c,bas));

//    Matrix3d bas;
  //  bas.setIdentity();
    FuncFunctord dasurf(xsvals,ysvals,surf2);
    Eigen::LevenbergMarquardt<FuncFunctord> lm2(dasurf);
    Eigen::LevenbergMarquardtSpace::Status status;
    status=lm2.minimize(c);
    if(status==1) return false;
  {  Matrix3d b; b.setIdentity(); p2s.setup(&c,&b);}
//else{    p2s.setup(&c,bas);}
  // coefficients=x;
  return true;
  
}



bool FitPolyLM2dIter(vector<triple> &points , P2Surface &p2s )
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
  // x^2 y^2 xy x y c 
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
  //   else( surf2.setup(&c,bas));

//    Matrix3d bas;
  //  bas.setIdentity();
    FuncFunctord dasurf(xsvals,ysvals,surf2);
    Eigen::LevenbergMarquardt<FuncFunctord> lm2(dasurf);
    Eigen::LevenbergMarquardtSpace::Status status;
    status=lm2.minimize(c);
    int j=0;
    converged=true;
    if(status==1) {return false;coff=3;converged=false;}
  //  cout << local.size() << endl;
   for (list <Vector3d> ::iterator it = local.begin(); it != local.end(); it++)
  {
    if(surf2.eval(xsvals.row(j))(2)-(ysvals.row(j))(2) > coff) {local.erase(it); it--;converged=false;}
    j++;
  }
   
  }
  {  p2s.setup(&c,&b);}
//else{    p2s.setup(&c,bas);}
  // coefficients=x;
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

  // coefficients=x;
  return true;
  
}


triple Tilt(Surface &s,triple &n, double x=0, double y=0)
{
  Vector3d Nx = s.unormv(Vector2d(x,y));
  triple N(Nx.data());
  if(N*n < 0) n=n*-1; 
  return n-N; ///(N*n)
}

bool Align(vector<triple> &director, Matrix3d &bas)
{
  
  for(int i=0;i!=director.size();i++)
  { Vector3d e=bas.transpose()*Vector3d(director[i].x,director[i].y,director[i].z);
//    cout << i<< endl;
    director[i].x=e(0);director[i].y=e(1);director[i].z=e(2);
  }
  return true;
}

double AvgStdvErr(list<double> &data,double*stderr=NULL,double *stdev=NULL, int window=10)
{
double avg=0;
double stdv=0;
double ba=0,ba2=0;
double mean=0;
list<double> blocks;

std::list<double>::iterator it;
    for (it = data.begin(); it != data.end(); ++it){
        avg+=*it;
    }
    mean=avg/data.size();
    if (stderr==NULL && stdev==NULL) return mean;
    int i=0;
   for (it = data.begin(); it != data.end(); ++it)
    {
        stdv+=(*it-mean)*(*it-mean);
        ba+=*it;
        if(i%(data.size()/window) == 0 && it!=data.begin()) {blocks.push_back((ba/(data.size()/window))); ba=0;} 
        i++;
    }
    stdv/=data.size();
    stdv=sqrt(stdv);
   for( it=blocks.begin();it!=blocks.end();it++)
   {
    ba+=*it;
    ba2+=(*it-mean)*(*it-mean);
   }
   ba/=blocks.size();
   ba2/=blocks.size();
   ba2=sqrt(ba2);
   if(stderr!=NULL) (*stderr)=ba2/sqrt(window);
   if(stdev!=NULL) (*stdev)=stdv;
   return mean;
}


bool FitTilt(vector<triple> &pos,vector<triple> &director, Surface &s,Plane &tus, Plane &tvs)
{
 Vector2d upos(pos[0].x,pos[0].y);
 triple N(s.unormv(upos).data());
     Vector3d sdu=s.du(upos);
 Vector3d sdv=s.dv(upos);


 triple eu(sdu.data());
 triple ev(sdv.data());

 
// s.VectorXYZ2UV((pos[0].p),upos);
  if(director[0]*N<0 ){s.flipnorm(); N=N*-1;}
 
/* eu.x=1;eu.y=0;eu.z=0;
 eu.x=0;eu.y=1;eu.z=0;
 N.x=0;N.y=0;N.z=1;*/
 vector<triple> tu(5);
 vector<triple> tv(5);
 for(int i=0;i!=5;i++)
 {
   upos(0)=(pos[i].x);
   upos(1)=(pos[i].y);
  N=(s.unormv(upos).data()); 
  triple &n=director[i];
  triple tilt=n-N;//(n/(N*n))-N
  tu[i]=pos[i];
  tu[i].z=eu*tilt;
  tv[i]=pos[i];
  tv[i].z=ev*tilt;
 }
 FitPlaneLM2d(tu,tus);
FitPlaneLM2d(tv,tvs);
return true;
}

Matrix2d TiltTensor(const Vector2d upos,Surface &tus, Surface &tvs, Geometry &g)
{
  // https://epje.epj.org.sci-hub.bz/articles/epje/abs/2000/11/e0025/e0025.html Eq. A6
  Matrix2d tmat;
  double ta=tus.eval(upos)[2];
  double tb=tvs.eval(upos)[2];
  // Set Basis to Diag
//  cout << ta <<" "<< tus.du(upos) << " GAMMA " <<  g.Gamma(upos,0,0,0)*ta+g.Gamma(upos,0,0,1)*tb << endl;
//  cout << tb <<" "<< tvs.dv(upos) << endl;
  
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
triple conc;
conc.x=0.0;
conc.y=0.0;
conc.z=1.0;
h.normfunc=&(VolumeNorm<double>);
ht.normfunc=&(SinVolumeNorm<double>);
list<double> vol;
list<double> avbend;
list<double> avup,avdown;
list<double> avk;
double zup=0.0, zdown=0.0;
vector<triple> *dir;
vector<triple> *c;
Vector2d cent;
cent=cent*0;
while(s.Valid())
{
//cout << "BUILIDING LIPIDS" << endl;
lipids.Evaluate(&s);
dir=lipids.GetDirectors();
c=lipids.GetCenters();
//cout << "GETTING ENV" << endl;
//#pragma omp parallel for
double AVGA=0;
double AVGPL=0;
double AVGRUP=0.0;int countup=0;
double AVGRDOWN=0.0;
double AVGBEND=0;
double AVGK=0;
int contrib=0;
srand48((long int) dir);
for(int i=0;i<c->size();i++)
{
 vector<triple> lpos(20); // 25
 vector<triple> ldir(6);
 atomlist indices;
 if(lipids.GetLipidEnv(i,lpos,indices,25,false))
 {
 //cout << lpos.size();
// cout << lpos[lpos.size()-1].abs() << endl;
 P2Surface s;
 Plane tvs,tus;
 Matrix3d basis; 
 PrincipalAxes(lpos,basis);
// basis.setIdentity();
 Align(lpos,basis);
 if(FitPolyLM2d(lpos,s))
 {

     
 if(witchcraft) // Try to estimate APL
 {
 double exc=drand48()*1.0-0.5;
 double cut=3.66; //3.7
 cut*=pow(EstimateScaling(s),2);
 double A=SurfaceIntegral(s,-cut+exc,cut+exc,-cut+exc,cut+exc,20);
 double PL=LipidsInArea(lpos,-cut+exc,cut+exc,-cut+exc,cut+exc);
 AVGA+=A;
 AVGPL+=PL;}
   contrib++;

// cout << APL << endl;
 // FitPoly2d(lpos,s);
//LM2d
// basis.setIdentity();
 Geometry go(&s);
 for(int j=0;j!=ldir.size();j++) {  Vector3d d((double*) &dir->at(indices[j]).p);d=basis.transpose()*d; ldir[j]=triple(d.data());  }
 FitTilt(lpos,ldir,s,tus,tvs);
 Matrix2d tab=TiltTensor(cent,tus,tvs,go); 
 Matrix2d bab=Btensor(cent,go);
// cout << "" << flush; //MAGIC 
 Matrix2d nab=-tab;//-tab;
 double t=(Tilt(s,ldir[0])).abs();
 double Splay=tab.trace();
 AVGBEND+=bab.trace();
 AVGK+=bab.determinant();
 if(abs(Splay) > 5 || Splay !=Splay ) 
 {
if(argc >4 && argv[4]=="t")  for(int j=0;j!=6;j++) cerr << ldir[j]*ldir[0] << " D " << lpos[j].x << " "<< lpos[j].y << " " << lpos[j].z << " I "  << indices[j] << endl; 
  
 }
 double KG=bab.determinant(); 
 double avec=dir->at(i)*conc;
 Vector3d spos=basis*s.eval(Vector2d(0.0,0.0));
 double height=spos(2)+c->at(i).z;
 if(avec>=0)  avup.push_back(height);
 else avdown.push_back(height);
 
 double angle=acos(avec);
  	  if(angle>PI/2) angle=PI-angle;
 
 //t=abs(angle);
h.AddDatapoint(1,Splay);
ht.AddDatapoint(1,t);//t



if(argc > 4) 
{ 
      Vector3d Nx = s.unormv(Vector2d(0.0,0.0));
      Nx=basis*Nx;
    cout << "KG : " << bab.determinant() << " B: " << go.GetH(cent) << " T: " << tab.trace() << " DL: "<< sqrt(Splay*Splay/4-KG) << " Frame " << avk.size() << " p  " << c->at(i).x << " " << c->at(i).y << " " <<c->at(i).z << " n " << dir->at(i).x << " " << dir->at(i).y << " " << dir->at(i).z << " N " << Nx(0) << " " << Nx(1) << " " << Nx(2) << endl;}
//else cout << flush;
}}}
double AVGAPL=AVGA/AVGPL;
//AVGAPL/=(1.0*contrib);
//cout << AVGAPL << " CONTRIB " << contrib << endl;

AVGBEND/=contrib;
AVGK/=contrib;
Matrix3 box;
box=traj.GetBoxMat();
triple dim;
dim.x=box.row[0].x;
dim.y=box.row[1].y;
dim.z=box.row[2].z;
//cout << AVGAPL  << " " << dim.x*dim.y/c->size()*2.0 << " " << contrib <<endl;
//cout << box.det() << endl;
if(!witchcraft) {
    AVGAPL= dim.x*dim.y/c->size()*2.0;
    
}
if(AVGAPL==AVGAPL) { vol.push_back(AVGAPL); avbend.push_back(AVGBEND); avk.push_back(AVGK);}
s=traj.NextFrame(); 
}


// FINAL EVALUATION
double stdev=0;
double ba=0;
double avg=AvgStdvErr(vol,&ba,&stdev,10);
double bend=AvgStdvErr(avbend);
double gauss=AvgStdvErr(avk);
double up=AvgStdvErr(avup);
double down=AvgStdvErr(avdown);

cout << "-----------------RESIS-1.01-Revise------------" << endl;
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
cout << " Area AVG [A^2]: " << avg*lipnum << " STDEV: " << stdev*lipnum << " Per Lipid [bilayer]: " << avg << " +- STDERR(10) " << ba  << " STDEV " << stdev << endl;  
cout << " Mean Curvature /  Lipid [nm^-1]: " << bend*10 << endl;
cout << " Gausian Curvature / Lipid [nm^-2]: " << gauss*100 << endl;
if(!witchcraft) cout << "Upper/Lower Pivotal Plane pos [A]: " << up << " / "<< down << endl;
ofstream of,of2;
of.open((string(argv[1]) +".spl").c_str());
of << "# Area AVG [A^2]: " << avg*lipnum << " STDEV: " << stdev*lipnum << " Per Lipid [bilayer]: " << avg << " +- STDERR(10) " << ba  << " STDEV " << stdev << endl;    
h.DropData(of);
vector<vector <double> > data;
h.Data(&data);
vector<double> res;

FitGParb(data,res,true);
//cout << "Monolayer Bending: kappa " << res[0]*2/(avg) << " kt Offset: " << res[2] << endl; 
double shift=res[2];
vector<vector <double> > tmp;
tmp.resize(2);
for(int i=0;i!=data[0].size();i++)
{
//  data[0][i]-=res[2];
if(data[0][i]>-0.05 && data[0][i]<0.05) {tmp[0].push_back(data[0][i]);tmp[1].push_back(data[1][i]);}
//  if(data[0][i]>-0.025 && data[0][i]<0.025) {tmp[0].push_back(data[0][i]);tmp[1].push_back(data[1][i]);}
}
//for(int i=0;i!=tmp[0].size();i++) cout<< tmp[0][i] << " " << tmp[1][i] << endl; 
FitGParb(tmp,res);
cout << "Monolayer Bending: kappa " << res[0]*2/(avg) << " kt Offset [nm^-1]: " << shift*10 << endl; 
/*for (int i=0;i!=data[0].size();i++)
{
  cout << data[0][i] << " "  << data[1][i] << endl;
}*/

of.close();
of2.open((string(argv[1]) +".itilt").c_str());
of2 << "# Area AVG [A^2]: " << avg*lipnum << " STDEV: " << stdev << " Per Lipid [bilayer]: " << avg*lipnum << " +- STDERR(10) " << ba  << " STDEV " << stdev << endl;    
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
//avg2/=norm;
/*cout << "XMAX " << xmax << endl;
for (int i=0;i!=data[0].size();i++)
{
  cout << tmp[0][i] << " "  << tmp[1][i] << endl;
}*/

//avt2/=(norm*norm);
//cout << std << " AVG " <<avg2 << endl ;
tmp[0].resize(0);
tmp[1].resize(0);
for(int i=0;i!=data[0].size();i++) if(data[0][i] > (xmax-std) &&  data[0][i] < (xmax+std) ) {tmp[0].push_back(data[0][i]);tmp[1].push_back(data[1][i]);}
cout << "avg" << avg2 << " " << std << " " <<  tmp[0].size() << endl;
FitGParb(tmp,res);
cout << "Bilayer Tilt: kappa " << res[0]*2/(avg/100) << " kt/nm^2 Offset: " << res[2] << endl; 
cout << "-------------------END SUMMARY--------------" << endl;

of2.close();

}
