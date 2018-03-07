#include <eigen3/Eigen/Dense>
#define PHAIL 0x0AFFE0000
using namespace std;
using namespace Eigen;


struct Surface: Function<Vector2d,Vector3d>
{
   Surface()
   {
     orient=1;
   }
   virtual Vector3d eval(const Vector2d &x) =0;
   virtual Vector3d diff(const Vector2d &x,const vector<int> &diff)=0;

   Vector3d du(const Vector2d &X,int u=1)
  {
    vector<int> da(2);
    da[0]=u;da[1]=0;
    return diff(X,da);
  }
    Vector3d duv(const Vector2d &X)
  {
    vector<int> da(2);
    da[0]=1;da[1]=1;
    return diff(X,da);
  }
 
  Vector3d dv(const Vector2d &X,int v=1)
  {
    vector<int> da(2);
    da[0]=0;da[1]=v;
    return diff(X,da);
  }
  Vector3d unormv(const Vector2d &X)
  {
    Vector3d a=this->du(X);
    Vector3d b=this->dv(X);
    a=a.cross(b);
//    a=basis*a;
    a=a/a.norm();
 //       cout << "NORMV " << a << endl;

    return a*orient;
  }
  void flipnorm()
  {orient*=-1;}
private: 
  int orient;
};

struct Sphere: public Surface
{
  
  Vector3d eval(const Vector2d &X)
  {
    Vector3d rv;
    rv[0]=radius*cos(X[0])*sin(X[1]);
    rv[1]=radius*sin(X[0])*sin(X[1]);
    rv[2]=radius*cos(X[1]);
    return rv+position;
  }
  Vector3d diff(const Vector2d &X,const vector<int> &diff)
  {
  //  cout <<"DIF" << diff[0] << " " <<diff[1] << endl;
    Vector3d rv;
    const int dx=diff[0],dy=diff[1];
  if(dx==1 && dy==0) {  rv[0]=radius*-sin(X[0])*sin(X[1]); rv[1]=radius*cos(X[0])*sin(X[1]) ; rv[2]=0;} 
  if(dx==2 && dy==0) {  rv[0]=radius*-cos(X[0])*sin(X[1]); rv[1]=radius*-sin(X[0])*sin(X[1]) ; rv[2]=0;} //+cos(X[0])*cos(X[0]))
  if(dx==0 && dy==1) {  rv[0]=radius*cos(X[0])*cos(X[1]); rv[1]=radius*sin(X[0])*cos(X[1]) ; rv[2]=radius*-sin(X[1]);}
  if(dx==0 && dy==2) {  rv[0]=radius*cos(X[0])*-sin(X[1]); rv[1]=radius*sin(X[0])*-sin(X[1]) ; rv[2]=radius*-cos(X[1]);}
  if(dx==1 && dy==1) {  rv[0]=radius*-sin(X[0])*cos(X[1]); rv[1]=radius*cos(X[0])*cos(X[1]) ; rv[2]=0;}
    return rv;
  }
  void setup(double r, Vector3d *b)
  {
    radius=r;
    position =*b;
    }
  
private:
  double radius;
  Vector3d position;
//  Matrix3d binv;
};



struct P2Surface: public Surface
{
  
  Vector3d eval(const Vector2d &X)
  {
    Vector3d rv;
    // Monge
    rv[0]=X[0];
    rv[1]=X[1];
    rv[2]=coeff[0]*X[0]*X[0]+coeff[1]*X[1]*X[1]+coeff[2]*X[0]*X[1]+coeff[3]*X[0]+coeff[4]*X[1]+coeff[5];
    rv=basis*rv;
    return rv;
  }
  Vector3d diff(const Vector2d &X,const vector<int> &diff)
  {
  //  cout <<"DIF" << diff[0] << " " <<diff[1] << endl;
    Vector3d rv;
    rv[0]=0; rv[1]=0;
    if(diff[1]>0) rv[0]=0;
    else if (diff[0]==1) rv[0]=1;
    if(diff[0]>0) rv[1]=0;
    else if (diff[1]==1) rv[1]=1;
    rv[2]=Diffz(X[0],X[1],diff[0],diff[1]);
 //   cout << basis << endl;
    rv=basis*rv;
    return rv;
  }
  bool VectorXYZ2UV(const Vector3d &X,  Vector2d &uv)
  {
  Vector3d rot;
  rot=basis*X;
  uv(0)=rot(0);
  uv(1)=rot(1);
  return true;
  }
  void setup(VectorXd *c, Matrix3d *b)
  {
    coeff=*c;
    basis=*b;
    binv=b->inverse();
  }
  
private:
  double Diffz(double x, double y, int dx=0, int dy=0)
  {
//    cout << " DIFF: DX " << dx << "DY " << dy << endl;
  if(dx==1 && dy==0) return  coeff[0]*2*x+coeff[3]+coeff[2]*y;
  if(dx==2 && dy==0) return  coeff[0]*2;
  if(dx==0 && dy==1) return  2*coeff[1]*y+coeff[2]*x+coeff[4];
  if(dx==0 && dy==2) return  2*coeff[1];
  if(dx==1 && dy==1) return  coeff[2];
  }
  VectorXd coeff;
  Matrix3d basis;
  Matrix3d binv;
};

struct Plane: public Surface
{
  // 
  Vector3d eval(const Vector2d &X)
  {
    Vector3d rv;
    // Monge
    rv[0]=X[0];
    rv[1]=X[1];
    rv[2]=coeff[0]*X[0]+coeff[1]*X[1]+coeff[2];
    rv=basis*rv;
    return rv;
  }
  Vector3d diff(const Vector2d &X,const vector<int> &diff)
  {
  //  cout <<"DIF" << diff[0] << " " <<diff[1] << endl;
    Vector3d rv;
    rv[0]=0; rv[1]=0;
    if(diff[1]>0) rv[0]=0;
    else if (diff[0]==1) rv[0]=1;
    if(diff[0]>0) rv[1]=0;
    else if (diff[1]==1) rv[1]=1;
    rv[2]=Diffz(X[0],X[1],diff[0],diff[1]);
 //   cout << basis << endl;
    rv=basis*rv;
    return rv;
  }
  bool VectorXYZ2UV(const Vector3d &X,  Vector2d &uv)
  {
  Vector3d rot;
  rot=basis*X;
  uv(0)=rot(0);
  uv(1)=rot(1);
  return true;
  }
  void setup(VectorXd *c, Matrix3d *b)
  {
    coeff=*c;
    basis=*b;
    binv=b->inverse();
  }
  
private:
  double Diffz(double x, double y, int dx=0, int dy=0)
  {
//    cout << " DIFF: DX " << dx << "DY " << dy << endl;
  if(dx==1 && dy==0) return  coeff[0];
  if(dx==2 && dy==0) return  0.0;
  if(dx==0 && dy==1) return  coeff[1];
  if(dx==0 && dy==2) return  0.0;
  if(dx==1 && dy==1) return  0.0;
  }
  VectorXd coeff;
  Matrix3d basis;
  Matrix3d binv;
};


bool PrincipalAxes(vector<tuple> &points, Matrix3d &axes, tuple *center=NULL)
{
  if(center != NULL) return false; // Always centered around zero
  // First POINT CENTER;
  Matrix3d rtemp, rd=Matrix3d::Zero(), rot;
  
  for(int i=0;i!=points.size();i++)
  {
   // points[i].print();
  /*  rtemp(0,0)=0; rtemp(0,1)=-points[i].z; rtemp(0,2)=points[i].y;
    rtemp(1,0)=points[i].z; rtemp(1,1)=0; rtemp(0,2)=-points[i].x;
    rtemp(2,0)=-points[i].y; rtemp(2,1)=points[i].x; rtemp(2,2)=0;*/
   rtemp(0,0)=points[i].y*points[i].y+points[i].z*points[i].z;
   rtemp(1,1)=points[i].x*points[i].x+points[i].z*points[i].z;
   rtemp(2,2)=points[i].y*points[i].y+points[i].x*points[i].x;
   rtemp(1,0)=-points[i].x*points[i].y;
   rtemp(2,0)=-points[i].x*points[i].z;
   rtemp(2,1)=-points[i].y*points[i].z;
    rd+=rtemp;
  }
 // cout << rd << endl;
  SelfAdjointEigenSolver<Matrix3d> es;
  es.compute(rd);
  axes=es.eigenvectors();
//  axes << 0, 0, 1 , 1,0,0 ,0,1,0 ;
  Vector3d evals=es.eigenvalues();
//  cout << evals << endl;
//  cout << axes << endl;
  return true;
}


bool FitPoly2d(vector<tuple> &points , P2Surface &p2s, Matrix3d *basis=NULL )
{
   Matrix3d bas;   
  if(basis==NULL) bas.setIdentity();
  else bas=*basis;
  if(points.size()!=6) return false;
  // x^2 y^2 xy x y c 
  MatrixXd vdm(6,6);
  VectorXd b(6),c(6);
  Vector3d point;
  for(int i=0;i!=6;i++)
   { 
      point(0)=points[i].x;point(1)=points[i].y; point(2)=points[i].z;
      //point=bas.transpose()*point;
      vdm(i,0)=point(0)*point(0);
      vdm(i,1)=point(1)*point(1);
      vdm(i,2)=point(0)*point(1);
      vdm(i,3)=point(0);
      vdm(i,4)=point(1);
      vdm(i,5)=1;
      b(i)=point(2);
   }
  FullPivHouseholderQR<MatrixXd> dec(vdm);
  dec.setThreshold(1e-99);
 // c = vdm.colPivHouseholderQr().solve(b);
  c=dec.solve(b);
 // cout << "H : " << vdm.determinant() << endl;
  p2s.setup(&c,&bas);
  // coefficients=x;
  return true;
  
}

bool FitPlane(vector<tuple> &points , Plane &ps, Matrix3d *basis=NULL )
{
   Matrix3d bas;   
  if(basis==NULL) bas.setIdentity();
  else bas=*basis;
  if(points.size()!=6) return false;
  // x^2 y^2 xy x y c 
  MatrixXd vdm(3,3);
  VectorXd b(3),c(3);
  Vector3d point;
  for(int i=0;i!=3;i++)
   { 
      point(0)=points[i].x;point(1)=points[i].y; point(2)=points[i].z;
      //point=bas.transpose()*point;
      vdm(i,0)=point(0);
      vdm(i,1)=point(1);
      vdm(i,2)=1;
      b(i)=point(2);
   }
  FullPivHouseholderQR<MatrixXd> dec(vdm);
  dec.setThreshold(1e-99);
 // c = vdm.colPivHouseholderQr().solve(b);
  c=dec.solve(b);
 // cout << "H : " << vdm.determinant() << endl;
  ps.setup(&c,&bas);
  // coefficients=x;
  return true;
  
}

struct Geometry
{
public:
  Geometry()
  {s=NULL; ev_gamma=false;}
  Geometry(Surface *surf)
  {
    setup(surf);
  }
  bool setup(Surface *surf)
  {
    s=surf;
    ev_gamma=false;
    xcur[0]=PHAIL;
    xcur[1]=PHAIL;
    return true;
  }
  bool GetI(const Vector2d& X, Matrix2d& Ia)
  {
    Ia(0,0)=s->du(X).dot(s->du(X));
    Ia(0,1)=s->du(X).dot(s->dv(X));
    Ia(1,0)=Ia(0,1);
    Ia(1,1)=s->dv(X).dot(s->dv(X));
   // Ia=&I;
    return true;
  }
  bool GetII(const Vector2d &X, Matrix2d& IIa)
  {
    Vector3d N=s->unormv(X);
    IIa(0,0)=s->du(X,2).dot(N);
    IIa(0,1)=s->duv(X).dot(N);
    IIa(1,0)=IIa(0,1);
    IIa(1,1)=s->dv(X,2).dot(N);
    
   // IIa=II;
  }
  bool GetShape(const Vector2d &X, Matrix2d &shape)
  {
    
    if(X!=xcur){
    this->GetI(X,I);
    this->GetII(X,II);}
    xcur=X;
    double denum=I(0,0)*I(1,1)-I(0,1)*I(0,1);
    // fF-eG
    S(0,0)= (I(0,1)*II(0,1)-II(0,0)*I(1,1))/denum;
    // eF-fE
    S(0,1)= (I(0,1)*II(0,0)-II(0,1)*I(0,0))/denum;
    // gF-fG
    S(1,0)= (I(0,1)*II(1,1)-II(1,0)*I(1,1))/denum;
    // fF-gE
    S(1,1)= (I(0,1)*II(0,1)-II(1,1)*I(0,0))/denum;
    shape=S;
    return true;
  }
  double GetK(const Vector2d &X)
  {
    if(X!=xcur){
    this->GetI(X,I);
    this->GetII(X,II);}
     xcur=X;
    return (II(0,0)*II(1,1)-II(0,1)*II(0,1))/(I(0,0)*I(1,1)-I(0,1)*I(0,1));
  }
  double GetH(const Vector2d &X)
  {
    if(X!=xcur){
    this->GetI(X,I);
    this->GetII(X,II);}
    xcur=X;
    // eG-2fF+gE / 2(EG-F¹)
    return (I(1,1)*II(0,0)-2.0*I(0,1)*II(0,1)+II(1,1)*I(0,0))/(I(0,0)*I(1,1)-I(0,1)*I(0,1)); // /2.0;    -
      
}
  double Gamma(const Vector2d &X,const int i,const int j,const int k)
  {
    // The Christoffel Symbols $Gamma^i_{kj} can be obtained by Do-Carmo p.176 
    if(ev_gamma==false || X!=xcur) 
      GammaMatrix(X);
    // The Matrix Calculation I*Gamma = b
 
    return (MGamma[k+j])[i];
  }
  void GammaMatrix(const Vector2d &X)
  {
    if(X!=xcur){
    this->GetI(X,I);
    this->GetII(X,II);}
    xcur=X;
  //  cout << "COMPUTE GAMMA" << endl;
   // Where b = < x_{jk} x_{u;v}> 
    Vector2d b,gamma;
    Vector3d ev;
    // 11
    for(int k=0;k!=2;k++)
    {
      for(int j=0;j!=2;j++)
      {
    if(k==1 && j==1) ev=s->duv(X);
    if(k==0 && j==0)   ev=s->du(X,2) ;
    if(j==0 && k==1) ev=s->dv(X,2);
    if(!(k==0 && j==1)) {
    b[0]=ev.dot(s->du(X));
    b[1]=ev.dot(s->dv(X));
  //  cout <<" B " << b << endl;
    gamma = I.colPivHouseholderQr().solve(b);
  //  cout << gamma;
  //  cout << "SETTING" << k+j << endl;
    MGamma[k+j]=gamma;}
      }
     }
    ev_gamma=true;
  }
private:
  Vector2d xcur;
  Surface *s;
  Matrix2d I,II,S;
  Vector2d MGamma[3]; // TOp indexed
  bool ev_gamma;
};
