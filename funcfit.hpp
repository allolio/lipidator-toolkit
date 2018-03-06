#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/NonLinearOptimization>
using namespace std;
using namespace Eigen;



template <class X,class Y> struct FitFunction: Function<X, Y>
{
// Data Are MatrixXd A (X) MatrixXd B (Y)
// X are the params
// virtual VectorXd diff(const VectorXd &x, const vector<int> &diff)=0;
// VectorXd eval(const VectorXd &x)
int setparam( const VectorXd &prams)
{
  m_params=prams;
  return 0;
}

virtual int dparam(const X &x,const vector<int> &diff, Y *y)=0; // Where X is the Xvector

VectorXd m_params;
int szparam;
};


/*template<typename _Scalar>
struct Functor
{
    typedef _Scalar Scalar;
    enum {
        InputsAtCompileTime = Eigen::Dynamic,
        ValuesAtCompileTime = Eigen::Dynamic
                          };
    typedef Eigen::Matrix<Scalar,100,1> InputType;
    typedef Eigen::Matrix<Scalar,10,1> ValueType;
    typedef Eigen::Matrix<Scalar,10,100> JacobianType;

    const int m_inputs, m_values;

    Functor() : m_inputs(InputsAtCompileTime), m_values(ValuesAtCompileTime) {}
    Functor(int inputs, int values) : m_inputs(inputs), m_values(values) {}

    int inputs() const { return m_inputs; } // number of degree of freedom (= 2*nb_vertices)
    int values() const { return m_values; } // number of energy terms (= nb_vertices + nb_edges)

    // you should define that in the subclass :
    //    void operator() (const InputType& x, ValueType* v, JacobianType* _j=0) const;
};*/


// Specialized functor warping the ikLikeCosts function
struct FuncFunctord 
{
   // typedef Eigen::AutoDiffScalar<Eigen::Matrix<Scalar,Eigen::Dynamic,1> > ADS;
   // typedef Eigen::Matrix<ADS, Eigen::Dynamic, 1> VectorXad;

    // pfirst and plast are the two extremities of the curve
    FuncFunctord(const Eigen::MatrixXd& targetx, const Eigen::MatrixXd& targety, FitFunction<VectorXd,VectorXd>& func )
        :   //Functor<double>(100, 10), 
            Xval(targetx), Yval(targety),
            m_func(func)
    {
      
    }
    int values() const { return Xval.cols()*Xval.rows(); } // number of degree of freedom (= 2*nb_vertices)
    int inputs() const { return m_func.szparam; } // number of energy terms (= nb_vertices + nb_edges)

    // input = x = {  ..., x_i, y_i, ....}
    // output = fvec = the value of each term
    int operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec)
    {
  //      cout << "SIZEY" << Yval.rows() << "FV" << fvec.size() << endl;
        fvec.resize(Yval.cols()*Yval.rows()); //LOL
	// cout << "CALLED :" << x << endl <<" " << fvec << endl;
//	fvec.fill(0.0f);
//        cout << "XVAL" << Xval  << endl; 
//	cout << "SIZEY" << Yval.rows() << "FV" << fvec.size() << endl;
	m_func.setparam(x);
	for(int i=0;i!=Xval.rows();i++)
       {
	 VectorXd out(Yval.cols());
	//fvec[i]=(Yval.row(i)[0]-
	out= m_func.eval(Xval.row(i)); // FULL POLYD
//	cout <<" ONE "  << out.size() << " "<<  Yval.row(i).size() << endl;
	out-= Yval.row(i);
	fvec.segment(i*out.size(),out.size())=out;
//	cout << Xval[i] << " " <<  Yval[i] << " "  <<  fvec[i] << endl;
       }
      
//	fvec=m_func.eval(x);
        return 0;
    }

    // Compute the jacobian into fjac for the current solution x
    int df(const Eigen::VectorXd &x, Eigen::MatrixXd &fjac)
    {
//       cout << "DF!" << x << endl;
        using namespace Eigen;
	vector<int> dvec(x.size());
	this->m_func.setparam(x);
	fjac.resize(Yval.cols()*Yval.rows(),x.size());
	//cout << "MRESI" << endl;
	int dsize=Yval.row(0).size();
	for(int i=0;i!=Xval.rows();i++)
	{
	    for(int j=0;j!=x.size();j++)
	    { for(int k=0;k!=dvec.size();k++) dvec[k]=0;
	      dvec[j]=1;
	    VectorXd out(dsize);
	    //VectorXd one=Xval[i];
	    m_func.dparam(Xval.row(i),dvec,&out);
	   // cout << out << endl;
	    fjac.block(i*dsize,j,dsize,1)=out;
	    }
	   // cout << "DONE" << i << " "<< i*dsize<< endl;
	}
//	cout << "JACOBIAN" << endl;
//	cout << fjac << endl;
        // copy the gradient of each energy term into the Jacobian
        return 0;
    }

    const MatrixXd &Xval;
    const MatrixXd &Yval;
    FitFunction<VectorXd,VectorXd>& m_func;
};

struct FitParabola : FitFunction<VectorXd, VectorXd>
{
  FitParabola()
  {szparam=2;}
  VectorXd eval(const VectorXd &X)
  {
    VectorXd f(1);
    f[0]=X[0]*X[0]*m_params[0]+m_params[1];
    return f;
  }
    VectorXd diff(const VectorXd& X, const vector<int>& df)
  {
    VectorXd f(1);
    f[0]=0;
    if(df[0]==1)    f[0]=2*X[0]*m_params[0];
    if(df[0]==2)    f[0]=2*m_params[0];
    return f;
  }
  int dparam(const VectorXd &x,const vector<int> &diff, VectorXd *rval)
  {
   double f=0;
   if(diff[0]==1) f=x[0]*x[0]; 
   if(diff[1]==1) f=1;
   (*rval)[0]=f;
  return 0;   
  }// Where X is the Xvector
};

struct FitConst : FitFunction<VectorXd, VectorXd>
{
  FitConst()
  {szparam=1;m_params.resize(1);}
  VectorXd eval(const VectorXd &X)
  {
    VectorXd f(1);
    f[0]=m_params[0];
    return f;
  }
    VectorXd diff(const VectorXd& X, const vector<int>& df)
  {
    VectorXd f(1);
    f[0]=0;
    if(df[0]==1)    f[0]=0;
    return f;
  }
  int dparam(const VectorXd &x,const vector<int> &diff, VectorXd *rval)
  {
   double f=0;
   if(diff[0]==1) f=1; 
   (*rval)[0]=f;
  return 0;   
  }// Where X is the Xvector
};


struct FitShParabola : FitFunction<VectorXd, VectorXd>
{
  FitShParabola()
  {szparam=3;m_params.resize(3);}
  VectorXd eval(const VectorXd &X)
  {
    VectorXd f(1);
    f[0]=(X[0]-m_params[2])*(X[0]-m_params[2])*m_params[0]+m_params[1];
    return f;
  }
    VectorXd diff(const VectorXd& X, const vector<int>& df)
  {
    VectorXd f(1);
    f[0]=0;
    if(df[0]==1)    f[0]=-2*(X[0]-m_params[2])*m_params[0];
    if(df[0]==2)    f[0]=2*m_params[0];
    return f;
  }
  int dparam(const VectorXd &x,const vector<int> &diff, VectorXd *rval)
  {
   double f=0;
   if(diff[0]==1) f=(x[0]-m_params[2])*(x[0]-m_params[2]); 
   if(diff[1]==1) f=1;
   if(diff[2]==1) f=-2*(x[0]-m_params[2])*m_params[0];

   (*rval)[0]=f;
  return 0;   
  }// Where X is the Xvector
};

struct P2FSurface: public FitFunction<VectorXd, VectorXd>
{
  P2FSurface()
  {
  //  m_params=VectorXd::Zero(6);
  }
  VectorXd eval(const VectorXd &X)
  {
    Vector3d rv;
    // Monge
    rv[0]=X[0];
    rv[1]=X[1];
    rv[2]=m_params[0]*X[0]*X[0]+m_params[1]*X[1]*X[1]+m_params[2]*X[0]*X[1]+m_params[3]*X[0]+m_params[4]*X[1]+m_params[5];
    rv=basis*rv;
    return rv;
  }
  VectorXd diff(const VectorXd &X,const vector<int> &diff)
  {
  //  cout <<"DIF" << diff[0] << " " <<diff[1] << endl;
    Vector3d rv;
    rv[0]=0; rv[1]=0;
    if(diff[1]>0) rv[0]=0;
    else if (diff[0]==1) rv[0]=1;
    if(diff[0]>0) rv[1]=0;
    else if (diff[1]==1) rv[1]=1;
    rv[2]=Diffz(X[0],X[1],diff[0],diff[1]);
    rv=basis*rv;
    return rv;
  }
  bool VectorXYZ2UV(const VectorXd &X,  Vector2d &uv)
  {
  VectorXd rot;
  rot=basis*X;
  uv(0)=rot(0);
  uv(1)=rot(1);
  return true;
  }
  void setup(VectorXd *c, Matrix3d *b)
  {
    m_params=*c;
    basis=*b;
    binv=b->inverse();
  }
  
  int dparam(const VectorXd &x,const vector<int> &diff, VectorXd *rval)
  {

/*    rv[0]=X[0];
    rv[1]=X[1];
    rv[2]=m_params[0]*X[0]*X[0]+m_params[1]*X[1]*X[1]+m_params[2]*X[0]*X[1]+m_params[3]*X[0]+m_params[4]*X[1]+m_params[5];
    rv=basis*rv;
 */
Vector3d rv;
   rv(0)=0;
   rv(1)=0;
    
   double f=0;
   if(diff[0]==1) f=x[0]*x[0]; 
   if(diff[1]==1) f=x[1]*x[1];
   if(diff[2]==1) f=x[0]*x[1];
   if(diff[3]==1) f=x[0];
   if(diff[4]==1) f=x[1];
   if(diff[5]==1) f=1;
   rv[2]=f;
   (*rval)=(basis*rv);
  return 0;   
  }// Where X is the Xvector
  
private:
  double Diffz(double x, double y, int dx=0, int dy=0)
  {
//    cout << " DIFF: DX " << dx << "DY " << dy << endl;
  if(dx==1 && dy==0) return  m_params[0]*2*x+m_params[3]+m_params[2]*y;
  if(dx==2 && dy==0) return  m_params[0]*2;
  if(dx==0 && dy==1) return  2*m_params[1]*y+m_params[2]*x+m_params[4];
  if(dx==0 && dy==2) return  2*m_params[1];
  if(dx==1 && dy==1) return  m_params[2];
  }
//  VectorXd m_;
  Matrix3d basis;
  Matrix3d binv;
};


struct P2FPlane: public FitFunction<VectorXd, VectorXd>
{
  P2FPlane()
  {
  //  m_params=VectorXd::Zero(6);
  }
  VectorXd eval(const VectorXd &X)
  {
    Vector3d rv;
    // Monge
    rv[0]=X[0];
    rv[1]=X[1];
    rv[2]=m_params[0]*X[0]+m_params[1]*X[1]+m_params[2];
    rv=basis*rv;
    return rv;
  }
  VectorXd diff(const VectorXd &X,const vector<int> &diff)
  {
  //  cout <<"DIF" << diff[0] << " " <<diff[1] << endl;
    Vector3d rv;
    rv[0]=0; rv[1]=0;
    if(diff[1]>0) rv[0]=0;
    else if (diff[0]==1) rv[0]=1;
    if(diff[0]>0) rv[1]=0;
    else if (diff[1]==1) rv[1]=1;
    rv[2]=Diffz(X[0],X[1],diff[0],diff[1]);
    rv=basis*rv;
    return rv;
  }
  bool VectorXYZ2UV(const VectorXd &X,  Vector2d &uv)
  {
  VectorXd rot;
  rot=basis*X;
  uv(0)=rot(0);
  uv(1)=rot(1);
  return true;
  }
  void setup(VectorXd *c, Matrix3d *b)
  {
    m_params=*c;
    basis=*b;
    binv=b->inverse();
  }
  
  int dparam(const VectorXd &x,const vector<int> &diff, VectorXd *rval)
  {

/*    rv[0]=X[0];
    rv[1]=X[1];
    rv[2]=m_params[0]*X[0]*X[0]+m_params[1]*X[1]*X[1]+m_params[2]*X[0]*X[1]+m_params[3]*X[0]+m_params[4]*X[1]+m_params[5];
    rv=basis*rv;
 */
Vector3d rv;
   rv(0)=0;
   rv(1)=0;
    
   double f=0;
   if(diff[0]==1) f=x[0]; 
   if(diff[1]==1) f=x[1];
   if(diff[2]==1) f=1.0;
   rv[2]=f;
   (*rval)=(basis*rv);
  return 0;   
  }// Where X is the Xvector
  
private:
  double Diffz(double x, double y, int dx=0, int dy=0)
  {
//    cout << " DIFF: DX " << dx << "DY " << dy << endl;
  if(dx==1 && dy==0) return  m_params[0];
  if(dx==2 && dy==0) return  0.0;
  if(dx==0 && dy==1) return  m_params[1];
  if(dx==0 && dy==2) return  0.0;
  if(dx==1 && dy==1) return  0.0;
  }
//  VectorXd m_;
  Matrix3d basis;
  Matrix3d binv;
};

//  min {1/2 || F(x(k)) + F'(x(k))(x-x(k))||^2}
// --> x(k+1) = x(k)-F'(x(k))°1 F(x(k)).
// Pseudoinverse °


