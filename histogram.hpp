#include <stdarg.h>
#define NO_GRID -1
#define OUT_OF_BOUND -2
#define ADD_OK 0
#define INCR 0 
#ifndef PI
#define PI 3.1415927
#endif

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}


using namespace std;
// 1
template <class T> class Histogram;
class Grid{
public:
  Grid()
  {
    ready=false;
    limits.clear();
    width.clear();
  }
vector<int> limits; // Bin Limits
vector< pair <double, double> > width;
vector<double > lengths;
bool ready;
void setup()
{
 lengths.resize(limits.size());
 //vector<double> lengths; 
  if(limits.size()!=width.size()) 
  {width.resize(limits.size());
  for(int j=0;j!=limits.size();j++) width[j]=make_pair<double,double>(0.0, (double) limits[j]);
    
  }
   for(int i=0;i!=limits.size();i++)
  {
   lengths[i]=fabs(width[i].first-width[i].second)/(double) limits[i];
  }
 
  ready=true;
}
};

template <class T> double WeightNorm( vector<int> *coord,double freq, Histogram<T> *h)
{
return freq;  
}

template <class T> double VolumeNorm( vector<int> *coord,double freq, Histogram<T> *h)
{
double vol=1;
for(int i=0; i!=h->grid.limits.size(); i++)
 {vol*=h->grid.lengths[i];}
// cout << vol << endl ;
 return vol;
}

template <class T> double CircVolumeNorm( vector<int> *coord,double freq, Histogram<T> *h)
{
double vol=1;
// Make dependent on radius!

// 2= 
//cout << " call!" << endl;
for(int i=0; i!=h->grid.limits.size(); i++)
 { int n=(*coord)[i];
   double r0=h->grid.width[i].first;
  // cout << " mo " << r0 << endl;
   double alpha=h->grid.lengths[i]*1;
  if(n==0) {vol=(r0+alpha)*(r0+alpha)*PI-(r0*r0)*PI;}
 else { vol=(2*n+1)*PI*alpha*alpha+(2*alpha)*r0*PI;}
  
//   cout <<" n " <<n*h->grid.lengths[i] << " " << vol << endl ;
}

 return vol;
}







template <class T> int EasyBin(T &data, vector<double> *coord, Histogram<T> *h)
{
  // for(int i=0;i!=coordinate.size();i++) cout << coordinate[i] << endl;
  
  Grid *g=h->getGrid();
  
  if(g==NULL) return NO_GRID;
  if(!g->ready) g->setup();
 
  vector <int> bin;
  bin.resize(coord->size());
  for(int i=0;i!=bin.size();i++)
    
  { double d=g->width[i].first;
    double c=coord->at(i);
    float binn=((c-d)/g->lengths[i]+0.5*INCR);
    if (binn<0) return OUT_OF_BOUND;
    bin[i]=binn;
    if(bin[i]>=g->limits[i] || bin[i] < 0) return OUT_OF_BOUND;
 //  cout << bin[i]  << " @" << bin[i]*g->lengths[i] << endl;
  }
  
//  cout << endl;
  h->AddDP(data,bin);
  return ADD_OK;
} 
 
/* DimensionRecursion(vector<int> &bin,vector<int> &bin2,vector<int> Histogram<T> *h, int dimension=0)
 {
   
   // 1d case [
   // coord[i]=(1-[a_i])
   
   for(int i=0;i!=bin.size();i++)
   {
     
   }  
 }
 */
 template <class T> int EquitableBin(T &data, vector<double> *coord, Histogram<T> *h)
{
  // for(int i=0;i!=(*coord).size();i++) cout << coord->at(i) << endl;
  
  Grid *g=h->getGrid();
  
  if(g==NULL) return NO_GRID;
  if(!g->ready) g->setup();
 
  vector <int> bin,bin2;
  vector <float> part;
  bin.resize(coord->size());
    bin2.resize(coord->size());
  part.resize(coord->size());

  for(int i=0;i!=bin.size();i++)
    
  { double d=g->width[i].first;
    double c=coord->at(i);
    float binn=((c-d)/g->lengths[i]+0.5*INCR);
    part[i]=(binn-(int) binn)/g->lengths[i]-0.5;
    bin2[i]=binn+sgn(part[i]); 
    if (binn<0) return OUT_OF_BOUND;
    bin[i]=binn;
   
    if(bin2[i]<0 || bin2[i]>g->limits[i]) bin2[i]=bin[i];
     if(bin[i]>=g->limits[i] || bin[i] < 0) return OUT_OF_BOUND;
 //   cout << bin[i]  << " @" << bin[i]*g->lengths[i] << endl;
  }
  unsigned long int countdown= pow(2,bin.size()), pointer=1,cbit;
  double norm=1;
  vector<int> address=bin;
  for(int j=0;j!=countdown;j++)
    {
      for(int k=0;k!=bin.size();k++)
       {
       cbit= ((countdown) & pointer) >> (pointer-1);
//       cout << cbit << " " << cbit*(abs(part[k]))+!cbit*(1-abs(part[k])) << endl;
       norm*= cbit*(abs(part[k]))+!cbit*(1-abs(part[k]));  //      (a_i)//(1-|a_i|)
       address[k] = bin2[k]*cbit+!cbit*bin[k];
  //     cout << address[k] << " " ;
       pointer = pointer << 1;
       }
     //  cout << "<ASr" << endl;    
       T d1=data*(T) norm;
   //    cout << d1 <<" X "<< norm << endl;
       h->AddDP( d1,address);    
       norm=1;
    }
//   cout <<"ADD" << endl;
  return ADD_OK;
} 
 


 
template <class T> int PeriodicBin(T &data, vector<double> *coord, Histogram<T> *h)
{
  // for(int i=0;i!=coordinate.size();i++) cout << coordinate[i] << endl;
  
  Grid *g=h->getGrid();
  
  if(g==NULL) return NO_GRID;
  if(!g->ready) g->setup();
 
  vector <int> bin;
  bin.resize(coord->size());
  for(int i=0;i!=bin.size();i++)
    
  { double d=g->width[i].first;
    double c=coord->at(i);
    float binn=((c-d)/g->lengths[i]+0.5*INCR);
    if (binn<0) 
      
     {
         coord->at(i)=g->width[i].second-(g->width[i].first-coord->at(i)); 
    //     cout << coord->at(i) << endl;
         return PeriodicBin<T>(data, coord, h);   
     }                     
    bin[i]=binn;
    if(bin[i]>=g->limits[i] ) 
    {
      
       coord->at(i)=(coord->at(i)-g->width[i].second)+g->width[i].first; 
	//  cout << coord->at(i) << endl;
	  if(coord->at(i)<(g->width[i].first+g->lengths[i]*0.001)) bin[i]=0;
          else return PeriodicBin<T>(data, coord, h);   
    }
//   cout << bin[i]  << " @" << bin[i]*g->lengths[i] << endl;
  }
  
//  cout << endl;
  h->AddDP(data,bin);
  return ADD_OK;
} 

template <class T> class Histogram
{
public:
  Histogram(void)
  {
    histodata=NULL;
    weightdata=NULL;
    T me;
    vector<double> d;
   // EasyBin<T>(me, &d, this);g
    binfunc=&(EasyBin<T>);
    normfunc=&(WeightNorm<T>);
    multiplier=1;
  }
  Grid *getGrid()
  {
    return &grid;
  }
  bool setGrid(Grid g)
  {
    if(histodata!=NULL) {free(histodata);if(weightdata!=NULL) free(weightdata);}
    unsigned long long int gsize=sizeof(T);
    for(int i=0;i!=g.limits.size();i++) gsize*=g.limits[i];
    histodata=malloc(gsize);
    int wsize=gsize/sizeof(T)*sizeof(double);
    weightdata=malloc(wsize);
    if(histodata==NULL) return false;
    memset(histodata,0,gsize);
    memset(weightdata,0,wsize);
    size=gsize;
    grid=g;
    if(!grid.ready) grid.setup();
    
    return true;
  }
 int AddDatapoint(T data, double cp1,...)
 {
   vector <double> coordinate;
   coordinate.clear();
   coordinate.resize(grid.limits.size());
   va_list ap;
   va_start(ap, cp1); 
   coordinate[0]=cp1;
   for(int i=1;i!=coordinate.size();i++){  coordinate[i]=va_arg(ap, double);}
   va_end(ap);
   // Binning!
   return (*binfunc)(data,&coordinate,this);
//   (*selector)(&snap,this)
 //  AddDP(data,bin);
  
 }
 
 
 T GetRBin(int  cp1,...)
  {
   vector <int> coordinate;
   coordinate.clear();
   coordinate.resize(grid.limits.size());
   va_list ap;
   va_start(ap, cp1); 
   coordinate[0]=cp1;
   for(int i=1;i!=coordinate.size();i++){  coordinate[i]=va_arg(ap, int);}
   va_end(ap);
   // Binning!
   return NGridPos(coordinate);
//   (*selector)(&snap,this)
 //  AddDP(data,bin);       
  }
  
  double N()
  {
    return (*weightpos);
  }
  
 double CalcWeightsum()
  {
    double sum=0;
    double *wp=(double *)weightdata;   
    for(long long int i=0;i!=size/sizeof(double);i++) {
       sum+=*wp; wp++; }
   return sum;
  }
  
  bool setBinning(int (*pt2Func)(T &temp,vector<double> *c, Histogram<T> *h))
{
  if(pt2Func==NULL) return false;
  binfunc=pt2Func;
  return true;
}

  bool setNormalizer(double (*pt2Func)(vector<double> *c, Histogram<T> *h))
  {
  if(pt2Func==NULL) return false;
  normfunc=pt2Func;
  return true;
}

  /*
  void DropData(ostream &of, bool normalize=true)
  {
    double *wp=(double*)weightdata;
    T* point=(T*)histodata;
    double divisor=1;
    //if(normalize) divisor=CalcWeightsum();
    for(int i=0;i!=size/sizeof(T);i++)
    {
      if(i%grid.limits[0]==0)
      {

	//cout << "HERE!" << endl;
	for(int j=0;j<=grid.limits.size()-1;j++)
         if(i%grid.limits[j]==0) {cout << endl;} else break;
      }
    if(normalize && (*wp)!=0) divisor=(*wp);
    else divisor=1;
    of << (*point)/divisor << " ";  
        point++;
	wp++;
  
    }
  }
  */
  void setMultiplier(double mult)
  {
    multiplier=mult;
  }
  
   void Data(vector <vector <T> > *out, bool normalize=true)
  {
    double *wp=(double*)weightdata;
    T* point=(T*)histodata;
    double divisor=1;
    vector<int> coord;

    coord.resize(grid.limits.size());
    out->resize(grid.limits.size()+1);
    for(int k=0;k!=out->size();k++)  out->at(k).resize(size/sizeof(T));
    // Fill Wit Zeros  ?

    for(int i=0;i!=size/sizeof(T);i++)
    { 
      for(int j=0;j!=coord.size();j++) {   
out->at(j)[i]=(coord[j]*grid.lengths[j]+0.5*grid.lengths[j]) +(grid.width[j]).first;}
          if(normalize && (*wp)!=0) divisor=(*normfunc)(&coord,(*wp),this);
          else divisor=1;
      out->at(out->size()-1)[i]= (*point)/divisor*multiplier ;  
      point++;
      wp++;
            coord[0]++;
    
      if((i+1)%grid.limits[0]==0 && i!=0)
      {
//	coord[0]=0; coord[1]++;
	for(int j=0;j!=grid.limits.size()-1;j++)
        if((i+1)%grid.limits[j]==0) { if(coord[j]%grid.limits[j]==0) {coord[j]=0;coord[j+1]++;}/*cout << " :"<<j<<"*";>*/} 
        else {break;}
      }
    }
  }
  
   void DropData(ostream &of, bool normalize=true)
  {
    of << endl;
    double *wp=(double*)weightdata;
    T* point=(T*)histodata;
    double divisor=1;
    vector<int> coord;
    coord.resize(grid.limits.size());
    // Fill Wit Zeros  ?
    for(int i=0;i!=size/sizeof(T);i++)
    { 
      for(int j=0;j!=coord.size();j++) of << (coord[j]*grid.lengths[j]+0.5*grid.lengths[j]) +(grid.width[j]).first << "\t";
          if(normalize && (*wp)!=0) divisor=(*normfunc)(&coord,(*wp),this);
          else divisor=1;
      of << (*point)/divisor*multiplier << "\t" << endl;  
      point++;
      wp++;
            coord[0]++;

      if((i+1)%grid.limits[0]==0 && i!=0)
      {
//	coord[0]=0; coord[1]++;
	for(int j=0;j!=grid.limits.size()-1;j++)
        if((i+1)%grid.limits[j]==0) { if(coord[j]%grid.limits[j]==0) {coord[j]=0;coord[j+1]++;} of  <<endl;/*cout << " :"<<j<<"*";>*/} 
        else {of << endl;break;}
      }
    }
  }
  
  bool DropCube(ofstream &of, bool normalize=true,bool skiphead=false)
  {
      bool rv=true; double a2au=1.88973;
      double *wp=(double*)weightdata;
    T* point=(T*)histodata;
    double divisor=1;
    vector<int> coord;
    coord.resize(grid.limits.size());
    
      if(coord.size()!=3) {cerr<< "Warning: CUBE FILE is supposed to be 3d" << endl; rv=false;}
      of << std::fixed;
      of << std::setw( 12 );
      of << std::setprecision(6);
      of << "LIBTRAJF GAUSSIAN CUBE FILE" <<endl;
      of << "Created From Histogram" <<endl;
      double zero=0.0;
      of <<"  1 " << zero << " " << zero << " " << zero << endl;
      for(int i=grid.lengths.size()-1;i>=0;i--)
      {
          of <<" "<< grid.limits[i] << " ";
        for(int j=0;j!=i;j++) of << zero << " ";
          of << grid.lengths[i]*a2au << " " ;
        for(int j=i;j!=grid.lengths.size()-1;j++) of << zero << " ";
         of << endl;
      }
            of <<"  1 " << 1.0 << " " << zero << " " << zero << " " << zero  << endl;

     for(int i=0;i!=size/sizeof(T);i++)
    { 
      //for(int j=0;j!=coord.size();j++) of << (coord[j]*grid.lengths[j]+0.5*grid.lengths[j]) +(grid.width[j]).first << "\t";
          if(normalize && (*wp)!=0) divisor=(*normfunc)(&coord,(*wp),this);
          else divisor=1;
      of << (*point)/divisor*multiplier*a2au << " ";  
      point++;
      wp++;
      if(i%6==5) of << endl;
            coord[0]++;
    }
      return rv;
  }
  
  ~Histogram()
  {
    if(histodata!=NULL) free(histodata);
    if(weightdata!=NULL) free(weightdata);
  }
  void AddDP(T &point, vector<int> &coord, double weight=1)
  {
    T* pos=GridPos(coord);
    (*pos)+=point;
    //cout << *pos << endl;
    (*weightpos)+=weight;
  }
  T* GridPos(vector<int> &coord)
  {
  T* point=(T*) histodata;
  weightpos=(double*) weightdata;
  for(int i=0;i!=coord.size();i++)
   {
     int increment=1;
     
    for(int j=1;j<=i;j++) increment*=grid.limits[j-1]; 
//     cout << "increment_ " <<increment << " "<< coord[i] << " "<< (increment)*(coord[i]) << endl;
     point+=(increment)*(coord[i]);   
     weightpos+=(increment)*coord[i];
   }
  return point;
  }
  
  T NGridPos(vector<int> &coord)
  {
   T* result=GridPos(coord);
   T res=(*result)/(*normfunc)(&coord,(*weightpos),this); 
   return res;    
  }
int (*binfunc)(T &temp,vector<double> *c, Histogram<T> *h);
double (*normfunc)(vector<int> *c,double freq,  Histogram<T> *h);

Grid grid;

protected:
  
private:
unsigned long long int size;
void* histodata;  
void* weightdata;  
double* weightpos;
double multiplier;

};

