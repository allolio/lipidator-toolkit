template <class X, class Y> struct Function
{
public:
virtual Y eval(const X &x ) =0;
virtual Y diff(const X &x,const vector<int> &dff)=0;
};


/** Little Math Helper
*
**/
unsigned long int fac(unsigned long int i)
{
unsigned long int counter,result=i;
if(i==0) return 1;
for(counter=i-1;counter!=0;counter--)
{
result=result*counter;
}
return result;
}

unsigned long int dfac(unsigned long int i)
{
unsigned long int counter,result=i;
if(i==0) return 1;
for(counter=i-1;counter>0;counter-=2)
{
result=result*counter;
}
return result;
}


namespace std
{
struct triple
{ 
  triple()
  {}
  triple(double *in)
  {
    x=*in;
    y=*(in+1);
    z=*(in+2);
  }
 union{
  double p[3];
  struct {
  double x;
  double y;
  double z;
  };
  };
    triple operator-(triple b)
  { triple t;
     t.x=this->x-b.x;
     t.y=this->y-b.y;
     t.z=this->z-b.z;
     return t;}
     triple operator%(triple b)
  { triple t;
    t.x=this->y*b.z-this->z*b.y;
    t.y=this->z*b.x-this->x*b.z;
    t.z=this->x*b.y-this->y*b.x;
     return t;}
     triple operator+(triple b)
  { triple t;
     t.x=this->x+b.x;
     t.y=this->y+b.y;
     t.z=this->z+b.z;
     return t;}
       triple operator*(double b)
  { triple t;
     t.x=this->x*b;
     t.y=this->y*b;
     t.z=this->z*b;
     return t;}
          triple operator/(double b)
  { triple t;
     t.x=this->x/b;
     t.y=this->y/b;
     t.z=this->z/b;
     return t;}
  
     // Scalar Product
        double operator*(triple b)
     { double t=0;
     t+=x*b.x;
     t+=y*b.y;
     t+=z*b.z;
     return t;}
   double abs()
   {
   return sqrt(x*x+y*y+z*z);
   }
   triple positive()
   {
   triple r;
   r.x=fabs(this->x);
   r.y=fabs(this->y);
   r.z=fabs(this->z);
   return r;
   }
   void print()
   {
   cout << x << " " << y << " " << z << endl;
   }
};


struct Matrix3
{
  triple row[3];
  double p(int i,int j)
  {
    return this->row[i].p[j];
  }
    void setp(int i,int j,double val)
  {
    this->row[i].p[j]=val;
  }

  Matrix3 operator*(double l)
  { 
    Matrix3 mat;
    mat.row[0]=this->row[0]*l;
    mat.row[1]=this->row[1]*l;
    mat.row[2]=this->row[2]*l;
    return mat;    
  }
    Matrix3 operator/(double l)
  { 
    Matrix3 mat;
    mat.row[0]=row[0]/l;
    mat.row[1]=row[1]/l;
    mat.row[2]=row[2]/l;
    return mat;    
  }

  triple operator*(triple x)
  { 
    triple a;
    a.x = row[0]*x;
    a.y = row[1]*x;
    a.z = row[2]*x;
    return a;
  }
  double det()
  {
    Matrix3 *m=this;
           double dval= m->p(0, 0) * (m->p(1, 1) * m->p(2, 2) - m->p(2, 1) * m->p(1, 2)) -
             m->p(0, 1) * (m->p(1, 0) * m->p(2, 2) - m->p(1, 2) * m->p(2, 0)) +
             m->p(0, 2) * (m->p(1, 0) * m->p(2, 1) - m->p(1, 1) * m->p(2, 0));
	     return dval;
  }
  Matrix3 inv()
  {
    // computes the inverse of a matrix m
double invdet = 1 / det();
Matrix3 *m=this;
Matrix3 minv; // inverse of matrix m
minv.setp(0, 0, (m->p(1, 1) * m->p(2, 2) - m->p(2, 1) * m->p(1, 2)) * invdet);
minv.setp(0, 1, (m->p(0, 2) * m->p(2, 1) - m->p(0, 1) * m->p(2, 2)) * invdet);
minv.setp(0, 2,(m->p(0, 1) * m->p(1, 2) - m->p(0, 2) * m->p(1, 1)) * invdet);
minv.setp(1, 0, (m->p(1, 2) * m->p(2, 0) - m->p(1, 0) * m->p(2, 2)) * invdet);
minv.setp(1, 1, (m->p(0, 0) * m->p(2, 2) - m->p(0, 2) * m->p(2, 0)) * invdet);
minv.setp(1, 2, (m->p(1, 0) * m->p(0, 2) - m->p(0, 0) * m->p(1, 2)) * invdet);
minv.setp(2, 0, (m->p(1, 0) * m->p(2, 1) - m->p(2, 0) * m->p(1, 1)) * invdet);
minv.setp(2, 1, (m->p(2, 0) * m->p(0, 1) - m->p(0, 0) * m->p(2, 1)) * invdet);
minv.setp(2, 2, (m->p(0, 0) * m->p(1, 1) - m->p(1, 0) * m->p(0, 1)) * invdet);
return minv;
  }
  void print(void)
  {
      row[0].print();
            row[1].print();
      row[2].print();
  }
Matrix3 transpose(void)
{
  Matrix3 m=(*this);
  m.setp(0,1 , this->p(1,0));
  m.setp(1,0 , this->p(0,1));
  m.setp(2,0 , this->p(0,2));
  m.setp(2,1 , this->p(1,2));
  m.setp(1,2 , this->p(2,1));
  m.setp(0,2 , this->p(2,0));
  return m;  
}
/*  void T(void)
  {
    Matrix3 temp;
    for(int i=0;i!=3;i++)
    {
      row[i];
    }
  }*/
  
};
}  

double ALegendre(int m,int l, double x)
{
    // P^m_l −l ≤ m ≤ l
    if(l<0) l++;
    unsigned int labs=abs(l);
    // P^m_m=
    int mabs=abs(m);
    double sfac=1.0;
    if(mabs>labs) return 0;
    if(m<0) sfac=pow(-1,mabs)*fac(labs-mabs)/(1.0*fac(labs+mabs));
    int sign=1;
    if(labs%2!=0) sign=-1;
    double Pmm=sign*dfac(2*mabs-1)*pow(1-x*x,0.5*mabs);
    if (mabs==labs) return Pmm*sfac;
    double Pmp=x*(2*labs+1)*Pmm*sfac;
    if (labs-mabs==1) return Pmp*sfac;
    int del=labs-mabs;
    for(int i=1;i!=del;i++)
    {
        double Ptmp=((2*(mabs+i)+1)*x*Pmp-(2*mabs+i)*Pmm)/(i+1.0);
        Pmm=Pmp;
        Pmp=Ptmp;
    }
    return Pmp*sfac;
}

complex <double> Yml(double phi,double theta, int m, int l)
{
const    double PI=3.14159265;

      complex<double> im ;im.real(0); im.imag(1);
      double pfac=sqrt((2*l+1)/4.0*PI*fac(l-m)/fac(l+m));

    return exp(im*(double) m*phi)*ALegendre(m,l,cos(theta))*pfac;
}

double Legendre(double x,int order)
{
double result=0;
for(int i=0;i<=order/2;i++)
{
result+=pow(-1,i)*fac(2*order-2*i)/fac(i)/fac(order-i)/fac(order-2*i)*pow(x,order-2*i); 
}
return result*1/pow(2,order);
}




namespace std
{
typedef vector< complex<double> > yvals;
/**
*  Fourier Transforms
**/
complex<double> betak( yvals &yx, int k)
{
  complex<double> result=0;
  double PI=3.14159265;
  int sdata=yx.size();
  complex<double> im = sqrt(complex<double> (-1));
  complex<double> exponent=exp( im * complex<double>( -2*PI*k/sdata));
  complex<double> cexponent=1;
for(int n=0;n!=sdata;n++)  {result+= yx[n]*cexponent; cexponent*=exponent;}
   return result;
}

yvals dft(yvals &yx)
{
  yvals retval;
  retval.resize(yx.size());
  for (int k=0;k!=yx.size();k++)
  retval[k]=betak(yx,k);
  return retval;  
}

yvals splity( yvals &yx, bool imp=false)
{
  yvals retval;
  retval.resize(yx.size()/2);
  for(int i=0;i!=retval.size();i++)
  {
    retval[i]=yx[i*2+imp];
  }
  return retval;
}

yvals r2fft( yvals &yx)
{
    double PI=3.14159265;

  complex<double> im; im.real(0); im.imag(1);
  int n=yx.size();
  complex<double> exponent=exp((im)*complex<double>( -2*PI/n));
  complex<double> cexponent=1;  
  yvals retval,a,b;
  retval.resize(n);
  if(n==1) return yx;
  else{
  yvals u=splity(yx);
  a= r2fft(u);
  u=splity(yx,1);
  b= r2fft(u);
  }
  for(int k=0;k!=n/2;k++)
  {
    retval[k]=a[k]+b[k]*cexponent;
    retval[k+n/2]=a[k]-b[k]*cexponent;
    cexponent*=exponent; // adjudst
  }
  return retval;
}

yvals ir2fft( yvals &yx)
{
    double PI=3.14159265;

  complex<double> im ;im.real(0); im.imag(1);
  yvals retval;
  retval.resize(yx.size());
  for(int k=0;k!=retval.size();k++)
  {
    retval[k]=conj(yx[k]);
  }
  retval=r2fft(retval);
  for(int k=0;k!=retval.size();k++)
  {
    retval[k]/=retval.size();
  }
  return retval;
}

int nextp2(long int N)
{
 int tpow=log2(N);
 int rval=pow(2,tpow);
 if (rval<N) rval=pow(2,tpow+1);
 return rval;
}

/**
 *  Rotates an 
 *  @@Ar
 * **/
std::triple Rotate(std::triple &vec, double angle, std::triple &x)
{
  std::triple a;
  a.x=1;a.y=1;a.z=1;
  double center=(1-cos(angle));
  Matrix3 rot;
  rot.row[0]=a;
  rot.row[1]=a;
  rot.row[2]=a;
  rot=rot*(double) center;
  double sa=sin(angle);
  double ca=cos(angle);
  // Build Rotation Matrix
  rot.row[0]=rot.row[0]*vec.x;
  rot.row[1]=rot.row[1]*vec.y;
  rot.row[2]=rot.row[2]*vec.z;
//
 //  cout << " c " << center << " a " << angle << " t" << rot.row[0].x<< endl; 
  
  rot.row[0].x=vec.x*rot.row[0].x+ca;
  rot.row[0].y=vec.y*rot.row[0].y-vec.z*sa; 
  rot.row[0].z=vec.z*rot.row[0].z+vec.y*sa;
  
 // cout << rot.row[0].x << " " << rot.row[0].y << " " << rot.row[0].z << "\t    | " << x.x  << endl;
  
  rot.row[1].x=vec.x*rot.row[1].x+vec.z*sa;
  rot.row[1].y=vec.y*rot.row[1].y+ca; 
  rot.row[1].z=vec.z*rot.row[1].z-vec.x*sa;

 // cout << rot.row[1].x << " " << rot.row[1].y << " " << rot.row[1].z <<  "\t    | " << x.y  << endl;

  rot.row[2].x=vec.x*rot.row[2].x-vec.y*sa;
  rot.row[2].y=vec.y*rot.row[2].y+vec.x*sa; 
  rot.row[2].z=vec.z*vec.z*(1-ca)+ca; 
 
 // cout << rot.row[2].x << " " << rot.row[2].y << " " << rot.row[2].z <<  " \t   | " << x.z  << endl;

 // cout << endl;
  
  triple res=rot*x;
//  cout << res.x <<"  " << res.y << "  " << res.z << " NORM " << res.abs() << endl; 
  return rot*x;
}
};

bool isPower2(unsigned long int x)
{
    return (x != 0) && ((x & (x - 1)) == 0);
}
