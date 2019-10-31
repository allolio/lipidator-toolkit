// (c) 2010-2017 CdA Christoph Allolio

#ifndef HAVE_CONFIG
#include"config.hpp"
#endif
#define HAVE_MD
#ifndef HAVE_FT
//#include "filetool.hpp"
#endif

#define XYZ_ERROR -1
#define NODATA_ERROR -1
#define BIN_WIDTH 0.02
#define PI  M_PI
//#define DEBUG



struct Atom
{
  Atom() {//pos=NULL;vel=NULL;
  extended=NULL; szextended=0;wascopied=false;
  type="N/A"; pos.x=-1e9; pos.y=pos.x; pos.z=pos.x; // Create Invalid defaultatom;
  }
  Atom(const Atom& p) {
    type = p.type;
    num = p.num;
    pos = p.pos;
    vel = p.vel;
    szextended=p.szextended;
    extended=p.extended;
    wascopied=true;
  //  cout << "Copied Atom" << endl;
   /* if(p.extended==NULL) extended=NULL;
    else
    { extended=malloc(p.szextended); 
      memcpy(extended,p.extended,szextended); }*/
                               }
string type;
int  num;
triple pos;
triple vel;
void *extended;
int szextended;
//    long int cout;
  
 ~Atom() {
      if(extended!=NULL && szextended !=0 && wascopied==false)
      {   cout << ((long int) extended) << endl;
       free(extended); // TODO Speicherzugriffsfehler
	 
      } 
  	extended=NULL;szextended=0;
    }
      private:
      bool wascopied;

};

Atom failatom;

class SnapShot
{
  friend class TrajecFile;
    friend class XYZFile;

  public:
    SnapShot()
    {
      natoms=XYZ_ERROR;
      isValid=false;
            pbc=false;

      natoms=vatoms.size();
       //cout << "Create" << endl;
    }
    SnapShot(vector<Atom> &atoms, bool hasvel=false, bool haspbc=false)
    {
     //  cout << "Create" << endl;
      natoms=atoms.size();
      isValid=true;
      vatoms=atoms;
      veloc=hasvel;
            pbc=haspbc;
       //if(haspbc) setPbc(true);
     // cout << "Constructed with " << natoms << " Atoms" << endl;
    }
     void SnapOverwrite(vector<Atom> &atoms, bool hasvel=false, bool haspbc=false)
    {
     //  cout << "Create" << endl;
      natoms=atoms.size();
      isValid=true;
      vatoms=atoms;
      veloc=hasvel;
            pbc=haspbc;

     // cout << "Constructed with " << natoms << " Atoms" << endl;
    }
    
    int GetAtoms() {      natoms=vatoms.size();
                          if(pbc) natoms-=3;
                          return natoms;}
    Atom& GetAtom(int n) {if(n<=(natoms+pbc*3) && n>=
    0 ) return vatoms[n]; else {return failatom;}}
    void AddAtom(Atom a) { if(!pbc) vatoms.push_back(a); isValid=true; natoms++;// Secure this TODO
    if(pbc) vatoms.insert(vatoms.end()-3,a);
      
    }
    bool Valid() { return isValid;}
    vector <int> GetIndicesOfType(string Type)
    {
      vector <int> indices;
      for(int i=0;i!=natoms;i++)
      {
      if(vatoms[i].type.compare(Type)==0) indices.push_back(i);
      }
      return indices;
    }
    vector<string> GetUniqueAtomTypes()
    {
    vector<string> Atomtypes;
      Atomtypes=GetAtomTypes();
         sort(Atomtypes.begin(), Atomtypes.end());
        vector<string>::iterator new_end = unique(Atomtypes.begin(), Atomtypes.end());
         int y=0; 
        for(vector<string>::iterator i=Atomtypes.begin();i!=new_end;i++) y++;
      Atomtypes.resize(y);
      return Atomtypes;
    }
    vector <string> GetAtomTypes(vector<int> *ofIndex=NULL)
    {
    vector<string> theNames;
      for(int i=0;i!=natoms;i++)
	if(ofIndex==NULL) theNames.push_back(vatoms[i].type);
       if(ofIndex!=NULL)
        for(int i=0;i!=(*ofIndex).size();i++) theNames.push_back(vatoms[i].type);
        
    return theNames;
    }
    bool setPbc(bool setpbc)
    {
      if(!pbc) {pbc=setpbc;
      natoms=-3;}
      if(pbc && setpbc==false) natoms+=3;
      return true;
    }
    bool AddVel(int index,triple &vel)
    {
      if(index>=vatoms.size()) return false;
      veloc=true;
      vatoms[index].vel=vel;
    }
    ~SnapShot()
    {//vatoms.clear();//cout << "dest" << endl;
   //  cout << "Destroyed Snapshot" << endl;
    }
    bool toStream(ostream &out)
    { 
      out.precision(8);
      out.setf(ios::fixed,ios::floatfield);  
      out << vatoms.size() << endl;
      if(headline=="") out << "XYZlib generated XYZ" << endl;
      else out<<headline << endl;
    for (int i=0;i!=vatoms.size();i++)
     {
       out << vatoms[i].type << "   \t" << vatoms[i].pos.x << "    \t" << vatoms[i].pos.y << "    \t" << vatoms[i].pos.z;
       if(veloc) out <<"   \t" << vatoms[i].vel.x << "   \t" << vatoms[i].vel.y << "    \t" << vatoms[i].vel.z;
       out << endl;
     }
     return true;
    }
      bool haspbc()
    {return pbc;}
    bool GetPbc(triple *boxdim)
    {
      if((vatoms[vatoms.size()-3].pos.x*vatoms[vatoms.size()-2].pos.y) ==0) { cerr << "Invalid Box" << endl; return false;}
      //triple pbc;
      boxdim->z=vatoms[vatoms.size()-1].pos.z;
      boxdim->y=vatoms[vatoms.size()-2].pos.y;
      boxdim->x=vatoms[vatoms.size()-3].pos.x;
      //pbc.print();
      //memcpy(&pbc,boxdim,sizeof(triple));
      return true;
        }
        
       bool GetFullPbc(Matrix3 *box)
    {
      if((vatoms[vatoms.size()-3].pos.x*vatoms[vatoms.size()-2].pos.y) ==0) { cerr << "Invalid Box" << endl; return false;} 
      box->row[0]=vatoms[vatoms.size()-3].pos;//.positive();
      box->row[1]=vatoms[vatoms.size()-2].pos;//.positive();
      box->row[2]=vatoms[vatoms.size()-1].pos;//.positive();
      return true;
        }
 
    bool hasvel()
    {return veloc; }
    bool IndextoStream(vector<int> index,ostream &out, bool XYZ=false, bool wantvel=false)
    { int cur;
   //   out << setiosflags(ios::fixed) << setprecision(6)
      out.precision(8);
      out.setf(ios::fixed,ios::floatfield);  
      out << index.size() << endl;
      if(XYZ) out << "XYZlib generated XYZ" << endl;
    for (int i=0;i!=index.size();i++)
     { cur=index[i];
       if(XYZ) out << vatoms[cur].type << "   \t";
       out << vatoms[cur].pos.x << "    \t" << vatoms[cur].pos.y << "    \t" << vatoms[cur].pos.z;
       if(veloc && wantvel) out << "    \t" << vatoms[cur].vel.x << "   \t" << vatoms[cur].vel.y << "    \t" << vatoms[cur].vel.z;
       out << endl;
     }
      
      return true;
    }
    string GetHeadline(void)
    {
      return headline;
    }
  private:
  friend class Trajectory;
  int natoms;
  bool veloc;
  vector<Atom> vatoms;
  bool isValid;
  bool pbc;
  string headline;
};

class TrajecFile
{
  public:
    TrajecFile(){hasloaded=false;frames=0;fpos=0;lsnaps.clear(); lpos=lsnaps.begin(); 
    readcache=false;}
   
    bool AddSnapshot(SnapShot &snap);
    virtual bool SaveTo(string Filename)=0;
    virtual bool Load(string filename);
    void Unload(void);
    int GetAtoms();
    SnapShot& GetNextSnap();
    SnapShot& GetCurrentSnap();
    SnapShot& GetSnap(int n);
    int GetSnapCnt(void);
    void First();
  protected:
    virtual bool ParseNextFrame() =0;
    virtual bool GotoFrame(int n) =0;
    virtual int  framesfromlines(int flength, bool hasextended);
    ifstream hfile;
    bool hasvelocities;
    bool hasloaded;
    int atoms;
    int fpos;
    string comment;
    int frames;
    typedef  pair<SnapShot,int> frame;
    list<frame> lsnaps;
    list<frame>::iterator lpos;
    bool readcache;
};


class XYZFile : public TrajecFile
{   
  public:
    XYZFile(string filename) {
    hasloaded=false;frames=0;fpos=0;lsnaps.clear(); lpos=lsnaps.begin(); 
    readcache=false;
    Load(filename);
    }
     XYZFile() {
    hasloaded=false;frames=0;fpos=0;lsnaps.clear(); lpos=lsnaps.begin(); 
    readcache=false;
    }
     bool SaveTo(string Filename);
  protected:
     bool ParseNextFrame();
     bool GotoFrame(int n);
  
};

class CPMDTrajFile : public TrajecFile
{   
  public:
    CPMDTrajFile();
    CPMDTrajFile(string fname, vector<string> atmlist);
     void setAtomTypes(vector<string> atmlist)
     {
       alist=atmlist;
     }
     bool Load(string filename); // File Must have 3 Lines!

     bool SaveTo(string Filename) { cout << "NYI My Frriend!" << endl; return false;};
  /*   ~CPMDTrajFile()
     {
       Unload();
     }*/
  protected:
     bool ParseNextFrame();
     bool GotoFrame(int n);
     int framesfromlines(int filelength, bool hasextended);
  private:
    vector<string> alist;
    // Internal conversion factors
    double bohr2a; // Conversion to Angstrom
    double au2v; // Conversion to Angstrom/fs
    bool hasforces;
};


class Trajectory
{
  public:
   Trajectory() {theTraj.clear(); spos=theTraj.begin();xyz=NULL;curate=false; addvirtual=false; fullpbc=false; triple dim; dim.x=1e9; dim.y=1e9;dim.z=1e9;setBoxDim(dim,false);
}
   Trajectory(TrajecFile *file)
   {
   xyz=file;curate=false;fullpbc=false;triple dim; addvirtual=false; dim.x=1e9; dim.y=1e9;dim.z=1e9;setBoxDim(dim,false);
   }
   void setBoxDim(triple dim, bool cur=true)
   {
    
    boxdim=dim; curate=cur; fullpbc=false; }
   void setBox(Matrix3 box, bool cur=false)
   { boxdim=boxdim*0; 
     boxmatrix=box;
     inverse=box.inv();
     fullpbc=true;
   }
   bool AddSnapShot(SnapShot s)
    {
    bool y=s.Valid();
    if(y) {theTraj.push_back(s);
    spos=theTraj.end();
    spos--;
      
    }
    return y;
    }
      inline  bool getCutback()
  {
    return curate;
  }

      SnapShot& Current()
    {
      SnapShot *s;
      if(xyz==NULL)
      {
      if(spos!=theTraj.end()) {
	s=&(*spos);
     } else s=new SnapShot;
	
      } else s=&(xyz->GetCurrentSnap());
         current = s;
       if (current->haspbc())
     {
       if(!fullpbc)
      current->GetPbc(&boxdim);
       else { current->GetFullPbc(&boxmatrix);inverse=boxmatrix.transpose().inv(); }
   //   cout << "BOXDIM" << boxdim.x << "cur " << curate << endl;
      
     }
      if(addvirtual) (vtatoms)(s,this,vtatomdata);
      if(curate) curateSnapshot();
      return *s;
    }
    int GetSnapCnt()
    {
      if(xyz!=NULL) return xyz->GetSnapCnt();
      else return theTraj.size();
    }
/*    SnapShot Current()
    {
    SnapShot s;
     if(xyz==NULL)
      {
      if(spos!=theTraj.end()) 
      { s=(*spos);
      }} else s=(xyz->GetCurrentSnap());
         current = &s;
      //if(!s.isValid) s= new *SnapShot;
      if(curate) curateSnapshot();
      return s;
    }*/
    SnapShot& NextFrame()
    {
      SnapShot *s;
      if(xyz==NULL)
      {
      if(spos!=theTraj.end()) 
      { s=&(*spos);
      spos++;
      } else {s=new SnapShot;}} else s=&(xyz->GetNextSnap());
         current = s;
	  if (current->haspbc())
       {
	       if(!fullpbc)
      current->GetPbc(&boxdim);
       else { current->GetFullPbc(&boxmatrix); inverse=boxmatrix.transpose().inv();}
        }
      if(curate) curateSnapshot();
      if(addvirtual) (vtatoms)(s,this,vtatomdata);
      return *s;
    }
    SnapShot& Frame(int n) 
    { SnapShot *s;
      int i=0; //Curation missing
      if(xyz==NULL) 
	for(spos=theTraj.begin();spos!=theTraj.end();spos++) 
	{i++;
	  if(i==n) { 
	  s=&(*spos); current=s; 	  
	  if (current->haspbc())
	    {
	       if(!fullpbc)
	       current->GetPbc(&boxdim);
	       else { current->GetFullPbc(&boxmatrix); inverse=boxmatrix.transpose().inv(); }
	    }
            if(addvirtual) (vtatoms)(current,this,vtatomdata);
      if(curate) curateSnapshot();
      return *s; } }
      else
      {
	 current=&(xyz->GetSnap(n));
	 if(current->haspbc())
	 {
	 if(!fullpbc) current->GetPbc(&boxdim);
         else { current->GetFullPbc(&boxmatrix); inverse=boxmatrix.transpose().inv();}
	      }}
         if(addvirtual) (vtatoms)(current,this,vtatomdata);
	 if(curate) curateSnapshot();
	 return *current;
    }
    void setCutback(bool cut)
    {
    curate=cut;
    }
    bool FullPbc(bool fpbc=true)
    {
      fullpbc=fpbc; curate=false;
      return fpbc;
    }
    bool BoxAngles()
    {
      return fullpbc;
    }
    void First() {// UGLY Remove
     if(xyz==NULL) spos=theTraj.begin(); 
     else xyz->First();
    }
    triple& GetDimensions()
    { return boxdim;}
    Matrix3 & GetBoxMat( Matrix3 **mat=NULL)
    {
      if(fullpbc)
      {
	if(mat!=NULL) (*mat)=&inverse;
	return boxmatrix;
      }
      else
      {cerr << "No Fullpbc Box" << endl;  	return boxmatrix;};
    }
      bool setVirtualAtoms(int (*pt2Func)(SnapShot *s,Trajectory *t, void *statdat), void *statdat)
    {
    if(pt2Func==NULL || statdat==NULL) return false;
    vtatoms=pt2Func;
    vtatomdata=statdat;
    addvirtual=true;
    return true;
    }
    ~Trajectory() { theTraj.clear();}
    bool SaveTo(string filename)
    {
      list<SnapShot>::iterator xpos;
      xpos=spos;
      First();
      ofstream out;
      out.open(filename.c_str());
      SnapShot snap;
      while(true)
     {
      snap=NextFrame();
      if(addvirtual) (vtatoms)(&snap,this,vtatomdata);
      if(snap.Valid()) snap.toStream(out);
      else {break;}
     }
     out.close();
     spos=xpos;
     return true;
     }
    private:
    void curateSnapshot(void)
    {
      if(fullpbc) {cerr << "Full PBC Activated Cutback is WRONG , FIXME" << endl;}
      // Correct Divide Off Box Vectors
    SnapShot snap=*current;
   
    
    for (int i=0;i<snap.natoms;i++)
      { 
	if(snap.vatoms[i].type!="PBC")
	{
	snap.vatoms[i].pos.x=fmod(snap.vatoms[i].pos.x,boxdim.x);
        snap.vatoms[i].pos.y=fmod(snap.vatoms[i].pos.y,boxdim.y);
	snap.vatoms[i].pos.z=fmod(snap.vatoms[i].pos.z,boxdim.z);
	if(snap.vatoms[i].pos.x<0) snap.vatoms[i].pos.x=boxdim.x+snap.vatoms[i].pos.x;
	if(snap.vatoms[i].pos.y<0) snap.vatoms[i].pos.y=boxdim.y+snap.vatoms[i].pos.y;
        if(snap.vatoms[i].pos.z<0) snap.vatoms[i].pos.z=boxdim.z+snap.vatoms[i].pos.z;
	}
      }
      *current=snap;
    }
  bool fullpbc;
  bool addvirtual;
  list<SnapShot> theTraj;
  list<SnapShot>::iterator spos;
  triple boxdim;
  Matrix3 boxmatrix;
  Matrix3 inverse;
  bool curate;
  int (*vtatoms)(SnapShot *c, Trajectory *t, void* statdat);
  void *vtatomdata;
  SnapShot *current;
  TrajecFile *xyz;
};

class Analyzer
{
  public:
    Analyzer(Trajectory *t)
    {traj=t; stride=BIN_WIDTH; boxdim=(traj->GetDimensions());  hassplit=false; // TODO: Select Smallest Vector
    }
    bool CalcTypeRdf(string typeA, string typeB)
    {
      SnapShot snap;
      traj->First();
      vector<int> indexa;
      vector<int> indexb;
      int atom1=0; int atom2=0; 
      snap=traj->NextFrame();
      indexa=snap.GetIndicesOfType(typeA);
      indexb=snap.GetIndicesOfType(typeB);
      return CalcIndexRdf(indexa,indexb);
    }
   typedef vector<int> atomlist;

    void ABAAngDist(atomlist centralatom, atomlist outeratoms, double Cutoff)
    {
     histogram.clear();
     vector <double> angles;
     // For each central Atom
     // 1.Get CN Closest Atoms;
     // 2. Get The coordination Shell atom A's closest Atom
     // 3. Avoid Double Counting!
     // 3. Calculate Angle and add to Histogram vector.
     // FIXED: Adjust Probability ---> uniform distribution
     SnapShot snap=traj->Current();
     while(snap.Valid())
     {
       
     for (int i=0;i!=centralatom.size();i++)
     {
      atomlist coord=GetAtomsTo(centralatom[i], outeratoms, Cutoff);
      list < pair < int, int > > outerpair;
      int neighbours=coord.size()-1;//-1;
      if(neighbours>=4) neighbours=3;
      if(neighbours>5) {neighbours=4; cerr << "overshell!" << endl;}
      if(neighbours>9) neighbours=5;
      Atom a=snap.GetAtom(centralatom[i]);
      for (int i2=0;i2!=coord.size();i2++)
	{
	atomlist near=GetClosestTo(coord[i2],coord, neighbours+1);// Jeder ist sich selbst der nächste
	//if(near[0]>near[1]) outerpair.push_back(make_pair< int, int> (near[0],near[1]));
	//else outerpair.push_back(make_pair< int, int> (near[1],near[0]));
	for(int i3=0;i3!=near.size();i3++)
	   if(coord[i2]!=near[i3]){
	 if(coord[i2] < near[i3] ) outerpair.push_back(make_pair(coord[i2],near[i3]));
	 else outerpair.push_back(make_pair(coord[i2],near[i3]));} 
	}
	 //}
      outerpair.sort();
      outerpair.unique();
      
      for (list< pair <int, int> >::iterator it = outerpair.begin(); it != outerpair.end(); it++)
	{
	  triple ab;
	  triple ac;
	//  cout << (*it).first << " " <<(*it).second << endl;

	  ab=boxdist(snap.GetAtom((*it).first).pos,a.pos);
	  ac=boxdist(snap.GetAtom((*it).second).pos,a.pos);
	  double abt=ab.abs();
	  double act=ac.abs();
	  ab=ab/ab.abs();
	  ac=ac/ac.abs();
	 
	  double angle=acos(ab*ac)*360/2/3.141592654;
	  if(angle<5)
	  {
	    cerr << "ERROR ! THERE IS SOMETHING unphysical going on ! check PBC!" <<endl;
	    snap.toStream(cerr);
	    cerr << "CENTER:" << centralatom[i] << " B " << (*it).first << "C " << (*it).second << "|AB|=" << abt << "|AC|=" << act << endl;
	    for(int x=0;x!=coord.size();x++)
	    {
	      cerr << snap.GetAtom(coord[x]).type << " " << coord[x] << endl; 
	    }
	    exit(-1);
	  }
	// #ifdef DEBUG
//	  cout << angle << endl;
	//  #endif
	 // if(angle < 150.0) 
	 angles.push_back(angle);
	}
//	cout << "Atom-Done" << endl;
      outerpair.clear();
	}
      
   //    cout << "Shell:" << coord.size() << endl;
     
     snap=traj->NextFrame();
    }
    normfactor=1;//(double)angles.size();
    stride=1;
    double tmax=180.0,tmin=0.0;
    // normfactor=1;
    // NORMALIZE BY 1/sin(alpha)
     ///((double) angles.size());
     nbins=181;
     histogram.resize(nbins+1);
     for(int i=0;i!=angles.size();i++)
     {
     histogram[angles[i]]+= 1;
     }
     // NORMALIZE BY 1/sin(alpha)
     double count=0;
     for(int i=0;i!=histogram.size();i++)
     {
    histogram[i]=histogram[i]/sin((i+0.5)*2*3.141592654/360); // avoid divergence
   
     count+=histogram[i];
     }
     for(int i=0;i!=histogram.size();i++)
     {
    histogram[i]=histogram[i]/count;
     }
    
    }
    
    vector <int> GetClosestTo (int index, vector<int> towhat,int howmany=1)
    {
    list< pair<double,int> > memory;
    //list<int> yours;
  //  memory.resize(howmany);
    SnapShot snap=traj->Current();
    Atom a = snap.GetAtom(index);
    // Fill Buffer Initially
    for(int i2=0;i2!=howmany;i2++)
    {
       Atom b=snap.GetAtom(towhat[i2]);
       triple dist=boxdist(a.pos,b.pos);
       double ab=dist.abs();
       memory.push_front(make_pair(ab,towhat[i2]));
    }
    memory.sort();
    
    for (int i=howmany;i!=towhat.size();i++)
     {
       Atom b=snap.GetAtom(towhat[i]);
       triple dist=boxdist(a.pos,b.pos);
       double ab=dist.abs();
       if( ab<memory.back().first )
       {
       memory.push_front(make_pair(ab,towhat[i]));
       memory.sort();
       memory.pop_back();
       //if (ab<memory.back().first) memory.push_back(make_pair<double,int>(ab,i)); memory.pop_front();
       }
     }
    vector<int> atomlist;
    atomlist.resize(howmany);
    int i=0;
    for (list< pair <double, int> >::iterator it = memory.begin(); it != memory.end(); it++)
    {atomlist[i]=(*it).second;i++;}
    return atomlist;  
    }
   
   atomlist GetAtomsTo (int index, atomlist ofwhat,double coff)
   {
   atomlist theAtoms;   
   SnapShot &snap=traj->Current();
   Atom &a=snap.GetAtom(index);
    for (int i=0;i!=ofwhat.size();i++)
     {
       Atom &b=snap.GetAtom(ofwhat[i]);
  //     double bx=dist.abs();
       triple dist=boxdist(b.pos,a.pos);
       double ab=dist.abs();
       
//        cout << ab << " " << bx << " " << endl;
       if(coff>ab){ theAtoms.push_back(ofwhat[i]); }
     }
     return theAtoms;
   }
   /**
    *  As GetAtomsTo, just with several Atoms & Returns Molecules
    **/
   
   atomlist GetMoleculesWithin(atomlist &closeto, atomlist &ofwhat, double coff)
   {
    
     if(!hassplit) SplitintoMolecules();
     
     SnapShot &snap=traj->Current(); 
     atomlist theMols;
      for(int j=0;j!=closeto.size();j++)
      {
	  Atom &a=snap.GetAtom(closeto[j]);
      for (int i=0;i!=ofwhat.size();i++)
     {
       Atom &b=snap.GetAtom(ofwhat[i]);
             triple dist=boxdist(b.pos,a.pos);
 
       double ab=dist.abs();
       if(coff>ab){ theMols.push_back(atmol[ofwhat[i]]); }
     }
      }
      
      // Collect all the Molecules
      sort(theMols.begin(),theMols.end());
      vector<int>::iterator new_end = unique(theMols.begin(), theMols.end());
      int y=0;
      for(vector<int>::iterator i=theMols.begin();i!=new_end;i++) y++;

      theMols.resize(y);
      return theMols;
   }
   
   SnapShot SnapShotfromIndices(atomlist &index)
   {
     SnapShot ret;
     SnapShot &cur=traj->Current();
     for(int i=0;i!=index.size();i++)
     {
       Atom a=cur.GetAtom(index[i]);
       ret.AddAtom(a);
     }
     return ret;
   }
   
   atomlist MoleculestoAtoms(atomlist &molecules)
   {
     atomlist molatoms;

   //cout << " atmol" << atmol.size() <<endl;
   //cout << " molecules " << molecules.size() <<endl;
   for(int i=0;i!=molecules.size();i++)
    for(int j=0;j!=atmol.size();j++)
    {
      if(molecules[i]==atmol[j]) molatoms.push_back(j);
    }
    return molatoms;
   }
   
    atomlist IndextoMolecules(atomlist &index)
   {
     if(!hassplit) SplitintoMolecules();
     atomlist molecules;
     molecules=index;
   //cout << " atmol" << atmol.size() <<endl;
   //cout << " molecules " << molecules.size() <<endl;
   for(int i=0;i!=index.size();i++) molecules[i]=atmol[index[i]];
       sort(molecules.begin(), molecules.end());
        vector<int>::iterator new_end = unique(molecules.begin(), molecules.end());
        molecules.erase(new_end,molecules.end());
    return molecules;
   }
   
   double GetCN(double distance)
    {
     
      double cn=0;
      if(histogram.empty()) return NODATA_ERROR;
      for(int i=0; i!=histogram.size()-1;i++)
      {
	if(i*stride<distance) cn+=histogram[i]/atom1/frames;
      }
      return cn;
    }
    void setIsoMasses(Isotopeinfo info)
    { masses=info;}
  
    
   double CalcIndexTemp(vector<int> &Index, Isotopeinfo *iso=NULL, int contraints=0)
   {
     const double kbol=1.3806504 * (10^-23);
     const double au2J=4.3597439 * (10^-18);
     const double cfac=315774.65;
     double Ekin,cmass, Temp;
     triple rvel; // Resulting Velocity
     double Ekincent=0; double totmasscent; // Resulting Kinetic Energy
     double cmassekin=0;double Ekininst=0;
     unsigned long int frames=0;
     if(iso==NULL) iso=&masses;
     SnapShot snap=traj->Current();
  //   snap.toStream(cout);
     Ekin=0;
          int k=0;
          cout << endl;
          cout << "step \t CMasvel \t CMassEkin[ha] \t Ekin[ha] \t <Ekin[ha]> \t T [K] \t <T[k]> \t <T_CMCorr[K]> " << endl;
     while(snap.Valid())	
     {  frames++;
     // Now
      totmasscent=0;
      rvel=rvel*0;
     for(int i=0;i!=Index.size();i++)
     {
     Atom a=snap.GetAtom(Index[i]);
     
     cmass=iso->GetIsoMass(a.type);
     // Calculate kinetic energy
     Ekininst+=0.5*cmass*(a.vel*a.vel)*1822.8884;
     //Ekin+=0.5*cmass*(a.vel*a.vel)*1822.8884;
     
     totmasscent+=cmass;
     rvel=rvel+a.vel;
   //  cout << "Vel" << a.vel.abs() << " " << rvel.abs() << endl;
     /* cout << Ekin*/;
       }
       cout << k << "\t";
     Ekin+=Ekininst;
     cout << rvel.abs()/Index.size() << "\t";
     cmassekin=0.5*totmasscent*(rvel.abs()*rvel.abs()/Index.size()/Index.size())*1822.8884 ;
     cout << cmassekin << "\t" ;
   //  cout << Index.size() << endl;
     Temp=Ekin/frames/(Index.size()-1-0.3333333*(double)contraints)*cfac*2/3;
     double Tinst=Ekininst/(Index.size()-1-0.3333333*(double)contraints)*cfac*2/3;  
     cout << Ekininst <<"\t" << Ekin/frames <<"\t";
     cout << Tinst << "\t" << Temp << "\t";
     cout << (Ekin-cmassekin)/frames/(Index.size()-1-0.3333333*(double)contraints)*cfac*2/3 << endl;
     Ekininst=0;
     k++;
     #ifdef DEBUG
     #endif
     snap=traj->NextFrame();
     }
     //Temp=Ekin/frames/(Index.size()+1);
     return Temp;
   }
   
   
   triple GetComVel(SnapShot *s=NULL)
   {
     if(s==NULL) s=&traj->Current();
     triple rvel;
     rvel=rvel*0;
     double totmass=0;
     for(int i=0;i!=s->GetAtoms();i++)
     {  
     Atom &a=s->GetAtom(i);
     double currmass=masses.GetIsoMass(a.type);
     rvel=rvel+a.vel*currmass;
     totmass+=currmass;
     }
     rvel=rvel/totmass;
     return rvel;
   }
 
  
   triple GetComPos(SnapShot *s=NULL)
   {
     if(s==NULL) s=&(traj->Current());
     triple rpos;
        rpos.x=0;rpos.y=0;rpos.z=0;
     double totmass=0;
     for(int i=0;i!=s->GetAtoms();i++)
     {  
     Atom &a=s->GetAtom(i);
     double currmass=masses.GetIsoMass(a.type);
     rpos=rpos+a.pos*currmass;
     totmass+=currmass;
     }
     rpos=rpos/totmass;
     return rpos;
   }
 
  triple GetIndComPos(atomlist &indices,SnapShot *s=NULL,bool wrap=false)
   {
     if(s==NULL) s=&(traj->Current());
     triple rpos,opos;
     rpos.x=0;rpos.y=0;rpos.z=0;
     double totmass=0;
     opos=s->GetAtom(indices[0]).pos;
     for(int i=0;i!=indices.size();i++)
     {  
     Atom &a=s->GetAtom(indices[i]);
     double currmass=masses.GetIsoMass(a.type);
     // Correct for PBC JUMP

    if(!wrap) rpos=rpos+a.pos*currmass; 
    else   rpos=rpos+(opos+boxdist(a.pos,opos))*currmass;
   //   cout <<" X "<<  a.pos.x << " " << currmass << " >" << rpos.x <<"< " ;
     totmass+=currmass;
     }
//     cout << totmass << " " << rpos.x << " BO " << endl;
     rpos=rpos/totmass;
     return rpos;
   }
      triple GetIndWAvgPos(atomlist &indices, vector<double> &masses,SnapShot *s=NULL,bool wrap=false)
   {
     if(s==NULL) s=&(traj->Current());
     triple rpos,opos;
     rpos.x=0;rpos.y=0;rpos.z=0;
     if(wrap) opos=s->GetAtom(indices[0]).pos;
     double totmass=0;
     for(int i=0;i!=indices.size();i++)
     {
     Atom &a=s->GetAtom(indices[i]);
  //   cout << "Y" << rpos.x << " " ;
     // Correct for PBC JUMP
      if(!wrap) rpos=rpos+a.pos*masses[i];
      else   rpos=rpos+(opos+boxdist(a.pos,opos))*masses[i];
  //    cout <<" X "<<  a.pos.x << " " << masses[i] << " >" << rpos.x <<"< " ;
           totmass+=masses[i];
     }
          rpos=rpos/totmass;

//     cout  << " " << rpos.x << " BO " << endl;
     return rpos;
   }

 
    triple GetIndAvgPos(atomlist &indices,SnapShot *s=NULL,bool wrap=false)
   {
     if(s==NULL) s=&(traj->Current());
     triple rpos,opos;
     rpos.x=0;rpos.y=0;rpos.z=0;
     if(wrap) opos=s->GetAtom(indices[0]).pos;
     for(int i=0;i!=indices.size();i++)
     {  
     Atom &a=s->GetAtom(indices[i]);
  //   cout << "Y" << rpos.x << " " ;
     // Correct for PBC JUMP
      if(!wrap) rpos=rpos+a.pos; 
      else   rpos=rpos+(opos+boxdist(a.pos,opos));
   //   cout <<" X "<<  a.pos.x << " " << currmass << " >" << rpos.x <<"< " ;
     }
//     cout << totmass << " " << rpos.x << " BO " << endl;
     return rpos/indices.size();
   }
 
   double GetIndMass(atomlist &indices,SnapShot *s=NULL)
   {
     if(s==NULL) s=&(traj->Current());
     double totmass=0;
     for(int i=0;i!=indices.size();i++)
     {  
     Atom &a=s->GetAtom(indices[i]);
     double currmass=masses.GetIsoMass(a.type);
     totmass+=currmass;
     }
     return totmass;
   }

   
 
 
 bool PrintForceXYZ ( ostream &hist)
    {
      SnapShot snap=traj->Current();
      
      vector <string> Atomtypes;
      vector <int> Atomindex;
        double cf=-0.529177249;
	cout << snap.GetAtoms() << endl;
	cout << "XYZlib generated Force XYZ" << endl;
	Atomtypes=snap.GetUniqueAtomTypes();
         for(int j=0;j!=Atomtypes.size();j++)
       {
       Atomindex=snap.GetIndicesOfType(Atomtypes[j]);
         for(int i=0;i!=Atomindex.size();i++)
	{
	Atom a=snap.GetAtom(Atomindex[i]); triple *tup;
	if(a.extended==NULL) {cerr << "No Forces Found" << endl; return false;}
        tup=(triple*) a.extended;
	hist << a.type  <<"  \t " << a.pos.x << "  \t " << a.pos.y << "  \t " << a.pos.z << "  \t " << tup->x*cf << "  \t " << tup->y*cf <<  "   \t " << tup->z*cf << "  \t" << endl;
	}
	
       }
	
      return true;
    }

int SplitintoMolecules(double coff=0.9, Bondinfo *bin=NULL)
{
  molatoms.clear(); // Inserted July 2012
  // 
// Get Snapshot as current snapshot from traj
  SnapShot snap=traj->Current();
  Bondinfo binfo;
  if(bin!=NULL) binfo= *bin;
 int num=0; 
#ifdef HAVE_BOOST
using namespace boost;
typedef property<vertex_name_t, std::string> VertexProperty;
typedef property<edge_weight_t, float> EdgeProperty;
typedef adjacency_list <vecS, vecS, undirectedS, VertexProperty, EdgeProperty> Graph;
dynamic_properties dp;

typedef  boost::graph_traits<Graph>::vertex_iterator   v;
    Graph G;
    property_map<Graph, vertex_name_t>::type name = get(vertex_name_t(), G);
    dp.property("label", get(edge_weight, G));
    dp.property("id", get(vertex_name, G));
  string pfileline;
  // get length of file:
  vector<string> filenames;
  vector<string>::iterator uniqueparams;

  // read filenames into string vector

  // Build Overlap
	list<string> fparam1,fparam2;
	int size1,size2;
	float ovrlp;

     for(unsigned int i=0;i!=snap.GetAtoms();i++){
		float shortcut=coff;int conn=-1;
	for(unsigned int i2=i+1;i2<snap.GetAtoms();i2++)
	{
	//fparam1=tio.ReadFileParams(filenames[i] + ".out");
	//fparam2=tio.ReadFileParams(filenames[i2] + ".out");
	
	Atom a=snap.GetAtom(i);
	Atom b=snap.GetAtom(i2);
        std::triple ab;
         ab=boxdist(a.pos,b.pos);
        double dist=ab.abs();
	// In den Graphen Einfügen:
	//	if(a.type!="H") {
  
        if( dist<=binfo.GetBondMax(a.type,b.type)) { add_edge(i, i2,EdgeProperty(ovrlp), G);conn=1;}
		std::vector<int>::size_type i;
	  
	//}
		
      //	else if(dist<shortcut) { conn=i2; shortcut=dist;}
		
	//if (conn!=-1)
	//add_edge(i, conn,EdgeProperty(ovrlp), G);
       }
       if (conn==-1)
	   add_edge(i, i,EdgeProperty(ovrlp), G);
       //  if(conn==-1) add_vertex(G);
    }
	std::vector<int> component(num_vertices(G));
	 num = connected_components(G, &component[0]);
/*	for(int i=0;i!=component.size();i++)
	{
	  Atom a=snap.GetAtom(i);
	cout << " Atom " << i << "( " << a.type << " ) " << " is in Molecule :" <<  component[i] << endl;
	}*/
//	cout << " A total of " << num << " Molecules found" << endl;
	
       atmol=component;
       vector<int> submol;
       // Split onto Map
       for(int j=0;j!=num;j++){
       for(int i=0;i!=component.size();i++)
       {
	 if(component[i]==j) submol.push_back(i);
       }
       molatoms.insert(make_pair(j,submol));    
       submol.clear();
       }
       hassplit=true;

#endif
       nmol=num;
       //Component Structure
       return num;

  }
  /**
  Returns the Atom indices of Atoms in one Molecule --> Tool to convert atmol information
  **/
    vector<int> GetMoleculeAtoms(int inMol)
    {
      vector<int> indices;
      indices.clear();
    if(!hassplit) SplitintoMolecules();
    int molof=atmol[inMol];
    for(int i=0;i!=atmol.size();i++)
      if(atmol[i]==molof) indices.push_back(i);
    return indices;
    }
    
    int GetMolOf(int index)
    {
      if(!hassplit) SplitintoMolecules();
      return atmol[index];
    }
    
    vector<int> GetMoleculeNo(int Mol)
    {
      vector<int> indices;
      indices.clear();
    if(!hassplit) SplitintoMolecules();
    for(int i=0;i!=atmol.size();i++)
      if(atmol[i]==Mol) indices.push_back(i);
    return indices;
    }
    double GetDist(int atoma, int atomb)
    {
      SnapShot snap=traj->Current();
      Atom a=snap.GetAtom(atoma);
      Atom b=snap.GetAtom(atomb);
      return boxdist(a.pos,b.pos).abs();
    }
    bool  GridVolumeNorm(vector<int> indexa, vector<double> &normvec, double strde=0.5)
    {
                SnapShot s=traj->Current();
        //stride=0.5;
        triple box=traj->GetDimensions();
        triple boxdots=box/strde;
        triple bins=boxdots;
               int counter=0;
        bins.x=round(bins.x);bins.y=round(bins.y);bins.z=round(bins.z);
        if(boxdots.x>1e4) {cerr << "Too Many Gridpoints" << endl; return false;}
        vector <double> weight(box.x/strde);
        cout << "stride " << strde << endl;
        bins.print();
    //    box.print();
     //   exit(0);
        while(s.Valid())
        {
        counter++;
        box=traj->GetDimensions();
        triple stride; stride.x=box.x/bins.x;stride.y=box.x/bins.y;stride.z=box.z/bins.z;
        for(int x=0;x<bins.x;x++)
        {
            for(int y=0;y<bins.y;y++)
            {
                for(int z=0;z<bins.z;z++)
                {
                    double mindist=1e9;
                    #pragma omp parallel for
                    for(int i=0;i<indexa.size();i++)
                    {   
                        triple t;
                        t.x=(x+0.5)*stride.x;t.y=(y+0.5)*stride.y;t.z=(z+0.5)*stride.z;
                        double dist=boxdist(s.GetAtom(indexa[i]).pos,t).abs();
#pragma omp critical
                        if(dist<mindist) mindist=dist;
                    }
//     cout<< (int)(mindist/stride) << endl;
                     weight[(int)(mindist/strde)]+=stride.x*stride.y*stride.z;
               //     cout << mindist << endl;
                }
            }
        }
        s=traj->NextFrame();
        cout << "*" << flush;
        }
        for(int i=0;i!=weight.size();i++)
        {
         weight[i]/=(double) counter;
        // weight[i]*=2.0;
       //  cout << (i+0.5)*strde << " \t " << weight[i] << endl;
        }
        normvec=weight;
        return true;
     }
    
    bool CalcSurfaceRdf(vector<int> &indexa, vector<int> &indexb,double strde=0.5)
    {
     atom1=indexa.size();
     atom2=indexb.size();
     stride=strde;
     double cutoff;
     boxdim=traj->GetDimensions();
     if(boxdim.x < boxdim.y && boxdim.x< boxdim.z) cutoff=boxdim.x/2;
      else if(boxdim.y < boxdim.x && boxdim.y< boxdim.z) cutoff=boxdim.y/2;
      else cutoff=boxdim.z/2;
      // Endinsert
     nbins=int(cutoff/stride)+1;
     frames=0;
        histogram.clear();
      histogram.resize(nbins);
            for(int i=0;i!=nbins;i++) histogram[i]=0;
       boxdim=boxdim*0;
      SnapShot snap=traj->Current(); 
      while(snap.Valid())	
      { //cout << "Frame Entered  Atoms:" << snap.GetAtoms() << endl;
        frames++;
        boxdim=boxdim+traj->GetDimensions();
	for (int j=0;j!=atom2;j++)
	{
	  Atom b=snap.GetAtom(indexb[j]);
	  double minim=1e9;
#pragma omp parallel for
	  for (int i=0;i<atom1;i++)
	  { Atom a=snap.GetAtom(indexa[i]);
	    triple diff=boxdist(a.pos,b.pos);
	    double r=diff.abs();
#pragma omp critical
	    if(minim>r) minim=r;
	  }
	    if(minim<cutoff) 
	      {
	      if(int(minim/stride)>(nbins)) {cout << "Histogram Fuckup! --- Crashing!" << endl; exit(1);}
	      else { int bin=int(minim/stride);   histogram[bin]++;}
	      }          
        } 
        snap=traj->NextFrame();
      }
       boxdim=boxdim/frames;
       double volume=boxdim.x*boxdim.y*boxdim.z;
       double region= cutoff*2;
    //   if(!same)
          normfactor = (volume)/
          (atom2*frames);
         //(region*region*region stride*stride*stride*
  //     cout <<" Pi " << PI << "Frames " << frames << " atom1 " <<  atom1 << " atom2 "<< atom2 <<  " Region " << region<< endl ; 

      // Test:
    }
    bool CalcIndexRdf(vector<int> &indexa, vector<int> &indexb, double strde=0.2, bool same=false)
    {
     atom1=indexa.size();
     atom2=indexb.size();
     stride=strde;
     // Inserted here May2012
     double cutoff;
     triple boxdim=traj->GetDimensions();
     if(boxdim.x < boxdim.y && boxdim.x< boxdim.z) cutoff=boxdim.x/2;
      else if(boxdim.y < boxdim.x && boxdim.y< boxdim.z) cutoff=boxdim.y/2;
      else cutoff=boxdim.z/2;
      // Endinsert
   //   cout << "CUtoff" << cutoff<< " " << nbins << endl;
     nbins=int(cutoff/stride)+1;
     frames=0;
        histogram.clear();
      histogram.resize(nbins);
            for(int i=0;i!=nbins;i++) histogram[i]=0;
     SnapShot snap; 
      boxdim=boxdim*0;
      while((snap=traj->NextFrame()).Valid())	
      { //cout << "Frame Entered  Atoms:" << snap.GetAtoms() << endl;
        boxdim=boxdim+traj->GetDimensions();
        frames++;
	for (int j=0;j!=atom1;j++)
	{
	  Atom a=snap.GetAtom(indexa[j]);
	  
	  for (int i=0;i!=atom2;i++)
	  { Atom b=snap.GetAtom(indexb[i]);
	    triple diff=boxdist(a.pos,b.pos);
	    double r=diff.abs();
	    if(r<cutoff) 
	      {
	      if(int(r/stride)>(nbins)) {cout << "Histogram Fuckup! --- Crashing!" << endl; exit(1);}
	      else { int bin=int(r/stride);   histogram[bin]++;}
	      }          
	  }
        }           
      }
       boxdim=boxdim/frames;
       double volume=boxdim.x*boxdim.y*boxdim.z;
    //   double region= cutoff*2; volume=region*region*region;
     if(!same) normfactor = (volume)/(4.0*PI*stride*stride*stride*atom2*atom1*frames);
     else
     { normfactor = (volume)/
      (2.0*PI*stride*stride*stride*atom2*atom1*frames);
      cout << "ATOMS ARE SAME " << endl;
     }
  //     cout <<" Pi " << PI << "Frames " << frames << " atom1 " <<  atom1 << " atom2 "<< atom2 <<  " Region " << region<< endl ; 
     return true;
      // Test:
    }
    /**
     * Use this for RDFs - it corrects for the Spherical Coordinates 1/r^2
     **/
    bool PrintHistogram (ostream &hist,bool cn=false, int dimens=3)
    {
        double kz=0;
      for(int i=0; i<histogram.size()-1;i++)
      {   double binvol;
          if(dimens==3) binvol=((i+1)*(i+1)*(i+1)-(i*i*i))/3.0;
          else if(dimens==2) binvol=((i+1)*(i+1)-(i*i));
    	hist << (i+0.5)*stride << "     " << histogram[i]*normfactor/binvol;
        if(cn==true) { kz+=histogram[i]/atom1/frames ;  hist << "     "  << kz;  }
    	hist << endl; //normfactor/((i+0.5)*(i+0.5))
      }
      return true;
    }
    /**
     * Do not Use this for RDFs - it uses rectangular coordinates - Everything but RDF use this
     **/
     bool PrintHistogramR (ostream &hist,double offset=0,vector<double> *dynorm=NULL)
    {
      for(int i=0; i<histogram.size()-1;i++)
      {if(dynorm==NULL)
    	hist << (i+0.5)*stride+offset << "     " << histogram[i]* normfactor << endl;
        else hist << (i+0.5)*stride+offset << "     " << histogram[i]* normfactor /dynorm->at(i) << endl;
      }
      return true;
    }
    
    bool PrintForces (ostream &hist)
    {
      traj->First();
      SnapShot snap; 
      while((snap=traj->NextFrame()).Valid())	
      {
	cout << "  >>>NewFrame<<<" << endl;
	for(int i=0;i!=snap.GetAtoms();i++)
	{
	Atom a=snap.GetAtom(i); triple *tup;
	if(a.extended==NULL) {cerr << "No Forces Found" << endl; return false;}
        tup=(triple*) a.extended;
	hist << a.type <<"  \t " << tup->x << "  \t " << tup->y <<  "   \t " << tup->z << "  \t" << endl;
	}
	
      }
    }
    /** Gives the Center of the Molecule, starting with the atom closest to the box center cpos
     **/
    void MolCenter(int mol, SnapShot &s, triple *cpos=NULL, double cutoff=0)
    {
      if(cutoff==0) cutoff=90E90;
      triple pos,closer;
      atomlist molatoms=GetMoleculeNo(mol);
      if(cpos==NULL) 
      {
	pos=boxdim/2;
      }
      else {pos=*cpos;}
      int theatom=0;
      double closest=5E90;
      SnapShot &s1=traj->Current();
      // Center on Atom closest to middle;
      for(int i=0;i!=molatoms.size();i++)
      { Atom a= s1.GetAtom(molatoms[i]);
        double dist=boxdist(a.pos,pos).abs();
	if(dist<closest)
	{closest=dist; theatom=molatoms[i];}
      }
      if(closest<=cutoff)
      {
      closer= s1.GetAtom(theatom).pos; 
      triple relbas=boxdist(closer,pos);
      for(int i=0;i!=molatoms.size();i++)
      { Atom &a= s1.GetAtom(molatoms[i]);
	pos=a.pos;
      triple fpos=boxdist(pos,closer)+relbas;
      a.pos=fpos;
      s.AddAtom(a);
      }}
    }
    
    SnapShot MolWrap(triple centerpos, double cutoff=0)
    {
      SnapShot result;
      if(!hassplit)SplitintoMolecules();
      for (int i=0;i!=nmol;i++)
      {
	MolCenter(i,result, &centerpos, cutoff);	
      }
      return result;
    }
   inline triple boxdist(triple a, triple b)
   {
     if(!traj->BoxAngles()) { return boxdiff(a-b);}
     Matrix3 *inv;
  /*   triple add100;
     add100.x=100;
     add100.y=add100.x;
     add100.z=add100.y;*/
     Matrix3 &box=traj->GetBoxMat(&inv);
   //  box=box.transpose();
  //   Matrix3 inv2=box.transpose().inv();
   //  inv=&inv2;
   //  box=box.transpose();
  
    // cout << "BOX" << endl;
   //  box.print();
  //   inv->print();
   //  inv->transpose().print();
    // Matrix3 check=box*inv;
    // check.print();
   //  cout << "INVERSE" << endl;
 //   inv->print();
/*     triple ax=a; ax.y=0;ax.z=0;
     triple ay=a; ay.x=0;ay.z=0;
     triple az=a; az.y=0;az.x=0;
     triple bx=b; bx.y=0;bx.z=0;
     triple by=b; by.x=0;by.z=0;
     triple bz=b; bz.y=0;bz.x=0;
*/
     triple as=(*inv)*a;
    // as=as+add100;
/*   as.x=as.x-round(as.x);
     as.y=as.y-round(as.y);
     as.z=as.z-round(as.z);*/

  //   triple as= ((*inv)*ax)+((*inv)*ay)+((*inv)*az);
 //   cout << "AS" << endl;
  //   as.print();
   //  cout << "BS" << endl;
  //  triple bs= ((*inv)*bx)+((*inv)*by)+((*inv)*bz);
    triple bs=(*inv)*b;
      //  bs=bs+add100;

 /*    bs.x=bs.x-round(bs.x);
     bs.y=bs.y-round(bs.y);
     bs.z=bs.z-round(bs.z);*/

   // (bs).print();
     triple dist=as-bs;
   //       cout << "DISTX" << endl;
   //  dist.print();

  
     dist.x=dist.x-round(dist.x);
     dist.y=dist.y-round(dist.y);
     dist.z=dist.z-round(dist.z);
    // cout << "DIST" << endl;
   //  dist.print();
     return box*dist;
     
   }
   inline triple boxdiff(triple a)
    {  
      triple &boxdim=traj->GetDimensions();
     if(!traj->getCutback())
      {
//	cout << "fmod" << endl;
      a.x=fmod(a.x,boxdim.x);
      a.y=fmod(a.y,boxdim.y);
      a.z=fmod(a.z,boxdim.z);
      }
   if(a.x >= 0.5*boxdim.x)
    a.x -= boxdim.x;
  else if(a.x < -0.5*boxdim.x)
    a.x += boxdim.x;
  if(a.y >= 0.5*boxdim.y)
    a.y -= boxdim.y; 
  else if(a.y < -0.5*boxdim.y)
    a.y += boxdim.y;

  if(a.z >= 0.5*boxdim.z)
    a.z -= boxdim.z; 
  else if(a.z < -0.5*boxdim.z)
    a.z += boxdim.z;
	return a;
    }
    
    int GetMolnum(void)
    {
      if(!hassplit) SplitintoMolecules();
      return nmol;
    }
     protected:
   // typedef pair<const int, vector<int> > atomap;
    Isotopeinfo masses;
    Trajectory* traj;
    Topology top;
    double stride;
    vector<double> histogram;
    vector<int> atmol;
    map <const int, vector<int> > molatoms;
    bool hassplit;
    int nbins;
    int atom1;
    int atom2;
    int frames;
    triple boxdim;
    double normfactor;
    int nmol;
};

bool TrajecFile::Load(string filename) // File Must have 3 Lines!
{
  int filelength=0;
  string testvel;
  vector<string> ntokens;
  hfile.open(filename.c_str());
  if(hfile.fail()) return false;
  hfile >> atoms; // Only Valid for XYZ beware
  getline(hfile,comment);
  bool hasextended=false;
  getline(hfile,testvel);
  getline(hfile,testvel);
  
//  if(strlen(testvel.c_str())==0) getline(hfile,testvel);
//  cout << testvel << " " << comment << endl;
  if(Tokenize(testvel, NULL, " \t")>6) hasvelocities=true;
    else hasvelocities=false;  
  if(Tokenize(testvel, NULL," \t")>9) hasextended=true;
       #ifdef DEBUG
      if(hasvelocities) cout << "Found Velocities" << endl;
      else cout << "Cant Find Velocities " <<  testvel << endl;
      #endif
   hfile.seekg(0,ios::beg); fpos=0;
   while ( !hfile.eof() )
   {
   getline(hfile, testvel);
   filelength++;
   }
   frames=framesfromlines(filelength-1,hasextended);
  hfile.clear();
  hfile.seekg(0,ios::beg); fpos=0;
  hasloaded=true;
  return true;    
}
int TrajecFile::framesfromlines(int filelength, bool hasextended)
{
 frames=filelength/(atoms+2);
  if((filelength%(atoms+2))!=0) {
      #ifdef DEBUG
      cout << "Broken Trajectory" << endl;
      cout << filelength << " and " << atoms << endl;
      #endif
    //     frames--;
      }
            return frames;

}
void TrajecFile::Unload(void )
{
hasloaded=false;
hfile.close();
}

int TrajecFile::GetSnapCnt(void )
{
  if(hasloaded) return frames;
    else return XYZ_ERROR;
}

SnapShot& TrajecFile::GetNextSnap()
{
  SnapShot *snap;    list<frame>::iterator elbeg=lsnaps.begin() ;
  if (fpos>=frames) return *(new SnapShot); // Outside Frame Limit
       for (lpos=lsnaps.begin(); lpos != lsnaps.end(); lpos++)
             if((*lpos).second==(fpos))
	     {readcache=true;
	      #ifdef DEBUG
		cout << "using cached frame" << endl;
	        cout << (*lpos).second <<": " << fpos << endl;
	      #endif
	       snap= &(*lpos).first;
	       fpos= (*lpos).second+1; 
	       return *snap;
	     }
   //if (lpos==lsnaps.end() || (*lpos).second<fpos) // Cached Frame N/A Forwarded Frame
   lpos=lsnaps.end();
   if(lsnaps.size()>CACHEFRAMES) { lsnaps.pop_front();} 
   if(!ParseNextFrame()) return *(new SnapShot); //  Cache Next Frame
   snap = &(*lpos).first;
   #ifdef DEBUG
   cout << "Provided next Snapshot - Validity: " << snap->Valid() << endl;
   #endif
   lpos++; 
   return *snap;
}

SnapShot& TrajecFile::GetCurrentSnap()
{
  SnapShot *snap;    list<frame>::iterator elbeg=lsnaps.begin() ;
    if (fpos>frames) return *(new SnapShot); // Outside Frame Limit

       for (lpos=lsnaps.begin(); lpos != lsnaps.end(); lpos++)
             if((*lpos).second==(fpos-1))
	     {readcache=true;
	      #ifdef DEBUG
		cout << "using cached frame" << endl;
		cout << (*lpos).second <<": " << fpos << endl;
	      #endif
	       snap= &(*lpos).first;
	       return *snap;
	     } else
	       #ifdef DEBUG
	       cout << (*lpos).second <<": " << fpos << endl
	       #endif
	       ;
	     //Current frame MUST_BE CACHED
	     // THIS IS ONLY FOR FIRST FRAME
	      #ifdef DEBUG
	     cout << "Warning: Current Frame " << fpos << " not cached - using first Frame!" << endl;
	      #endif 
	     snap=&GetSnap(0);
	     lsnaps.insert(lpos,make_pair(*snap,-1));

	     return *snap;
   //if (lpos==lsnaps.end() || (*lpos).second<fpos) // Cached Frame N/A Forwarded Frame
   lpos=lsnaps.end();
   if(lsnaps.size()>CACHEFRAMES) { lsnaps.pop_front();} 
   if(!ParseNextFrame()) return *(new SnapShot); //  Cache Next Frame
   snap = &(*lpos).first;
   fpos--;
   #ifdef DEBUG
   cout << "Provided next Snapshot - Validity: " << snap->Valid() << endl;
   #endif
   //lpos++; 
   return *snap;

}

int TrajecFile::GetAtoms()
{
if(hasloaded) return atoms;
else return XYZ_ERROR;
}


bool TrajecFile::AddSnapshot(SnapShot &snap)
{
if(frames==0) atoms=snap.GetAtoms();
else if(atoms!=snap.GetAtoms()) return false;
frames++;
lsnaps.push_back(make_pair(snap,frames));
return true;
}


SnapShot& TrajecFile::GetSnap(int n)
{
  SnapShot *snap;
  if (n>frames) return *(new SnapShot); // Outside Frame Limit
     for (lpos=lsnaps.begin(); lpos != lsnaps.end(); lpos++)
      if((*lpos).second==n) { fpos=(*lpos).second+1; snap=&(*lpos).first;  cout << "Frame was cached! (of)" << lsnaps.size() << endl; return *snap; }
	     if(GotoFrame(n)) ParseNextFrame(); //  Cache Next Frame
	       else return *(new SnapShot);
   snap = &(*lpos).first;
   return *snap;
}

void TrajecFile::First()
{
  GotoFrame(0);
}

bool XYZFile::GotoFrame(int n)
{ string buffer;
  if (n>frames ) return false;
  if(fpos>n && hasloaded) {  hfile.clear();
    hfile.seekg(0,ios::beg); fpos=0;  }
    else if(!hasloaded) return false;
  while(fpos<n)
  {
    for(int i=0; i!=(atoms+2);i++) getline(hfile,buffer);
    fpos++;
  } 
  return true;
}

bool XYZFile::SaveTo(string filename)
{
  ofstream out;
  out.open(filename.c_str());
  GotoFrame(0); SnapShot snap;
  while(true)
  {
    snap=GetNextSnap();
    if(snap.Valid()) snap.toStream(out);
    else {break;}
  }
 out.close();
return true;
}

bool XYZFile::ParseNextFrame()
{
string buffer;
if(readcache)  { GotoFrame(fpos);
#ifdef DEBUG
cout << "Last Frame from Cache!" << endl;
#endif
readcache=false;}
// Ignore First two lines
getline(hfile,buffer);
getline(hfile,buffer);
int htokens;

if(fpos>frames) return false;
if(!hfile.is_open()) return false;
#ifdef DEBUG
 cout << "parsing new frame" << endl;
#endif
SnapShot snap;
 lsnaps.insert(lpos,make_pair(snap,fpos));
  lpos--;

 SnapShot &sn=(*lpos).first;
sn.headline=buffer;
sn.vatoms.resize(atoms);
sn.veloc=hasvelocities;
sn.natoms=atoms;

vector<string> filelines; filelines.resize(atoms);
for(int i=0;i!=atoms;i++)   getline(hfile,filelines[i]); 

	#pragma omp parallel for
for(int i=0;i<atoms;i++)
  {
   vector<string> line;
   
   htokens=Tokenize(filelines[i],&line, " \t");
 //  if(htokens<3) { cout << "Broken Frame ! Cannot Parse" << endl; return false;}
   sn.vatoms[i].type=line[0]; 
   sn.vatoms[i].pos.x=atof(line[1].c_str());
   sn.vatoms[i].pos.y=atof(line[2].c_str());
   sn.vatoms[i].pos.z=atof(line[3].c_str());
   if(hasvelocities && htokens > 3)
   {
    // triple *vel=new(triple);
   sn.vatoms[i].vel.x=atof(line[4].c_str());
   sn.vatoms[i].vel.y=atof(line[5].c_str());
   sn.vatoms[i].vel.z=atof(line[6].c_str());
   }
   line.clear();
  }
  
  //  cout <<" Cache: " << lsnaps.size() << "parsed frame " << fpos << endl;
// SnapShot snap(satoms,hasvelocities);
sn.isValid=true;
  
 fpos++;
 return true;
}

int CPMDTrajFile::framesfromlines(int filelength, bool hasextended)
{
 frames=filelength/atoms;
  if((filelength%atoms)!=0) {
      #ifdef DEBUG
      cout << "Broken Trajectory" << endl;
      #endif
    //     frames--;
      }
      hasforces=hasextended;
      return frames;
}

bool CPMDTrajFile::ParseNextFrame()
{
string buffer;
if(readcache)  { GotoFrame(fpos);
#ifdef DEBUG
cout << "Last Frame from Cache!" << endl;
#endif
readcache=false;}
// Ignore First two lines
int htokens;

if(fpos>frames) return false;
if(!hfile.is_open()) return false;
#ifdef DEBUG
 cout << "parsing new frame" << endl;
#endif
vector<Atom> satoms;

  
satoms.resize(atoms);
for(int i=0;i!=atoms;i++)
  {
   vector<string> line;
   
   getline(hfile,buffer);
   if(buffer.find("<<<<<<  NEW DATA  >>>>>>")!=string::npos) getline(hfile,buffer);
   htokens=Tokenize(buffer,&line);
   if(htokens<3) return false;
   satoms[i].type=alist[i]; 
   satoms[i].pos.x=atof(line[1].c_str())*bohr2a;
   satoms[i].pos.y=atof(line[2].c_str())*bohr2a;
   satoms[i].pos.z=atof(line[3].c_str())*bohr2a;
//   cout<< buffer<< endl
   if(hasvelocities && htokens > 3)
   {
   //    if(hasvelocities) cout << "Writing velocities" << endl;

    // triple *vel=new(triple);
   satoms[i].vel.x=atof(line[4].c_str())*au2v;
   satoms[i].vel.y=atof(line[5].c_str())*au2v;
   satoms[i].vel.z=atof(line[6].c_str())*au2v;
   }
   if(hasforces && htokens > 9)
   { //cout << "Extended Data = ?" << endl;
     triple t; // Force is in Atomic units
     t.x=atof(line[7].c_str());
     t.y=atof(line[8].c_str());
     t.z=atof(line[9].c_str());
     satoms[i].szextended=sizeof(triple);
     void *force=malloc(sizeof(triple)); // Extended Quantity in CPMD Trajectories is Force!
     memcpy(force,(void*) &t,sizeof(t));
     satoms[i].extended=force;
   }
 
   line.clear();
  }
  #ifdef DEBUG
  cout <<" Cache: " << lsnaps.size() << "parsed frame " << fpos << endl;
  #endif
 //SnapShot snap;
 
 lsnaps.insert(lpos,make_pair(SnapShot(satoms,hasvelocities),fpos));
  
  lpos--;
 fpos++;
 return true;
}

CPMDTrajFile::CPMDTrajFile(string filename, vector< string > atomlist)
{
    hasloaded=false;frames=0;fpos=0;lsnaps.clear(); lpos=lsnaps.begin(); 
    readcache=false;    
    bohr2a = 0.529177249; alist=atomlist; au2v=1; // Keep AU//au2v=0.045710289;
    Load(filename);
    if(atoms!=alist.size()) cout << " Provided Atoms do not Match!" << endl;
}
CPMDTrajFile::CPMDTrajFile()

 {
    hasloaded=false;frames=0;fpos=0;lsnaps.clear(); lpos=lsnaps.begin(); 
    readcache=false;    
    bohr2a = 0.529177249; alist.clear();  au2v=1; // Lets keep the a.u.//au2v= 0.045710289;
    }
    
bool CPMDTrajFile::GotoFrame(int n)
{
 string buffer;
  if (n>frames ) return false;
  if(fpos>n && hasloaded) {  hfile.clear();
    hfile.seekg(0,ios::beg); fpos=0;  }
    else if(!hasloaded) { cout << "TRAJECTORY not loaded!" << endl;
    return false; }
  while(fpos<n)
  {
    for(int i=0; i!=(atoms);i++) 
              {
		getline(hfile,buffer); 
		if(buffer.find("<<<<<<  NEW DATA  >>>>>>")!=string::npos) 
		                    {
				      cerr << "Warning: Fractured Traj" << endl;i--;
				    }
	      }
    fpos++;
  } 
  return true;
}

bool CPMDTrajFile::Load(string filename) 
{
  int filelength=0;
  string testvel;
  vector<string> ntokens;
  hfile.open(filename.c_str());
  if(!hfile) return false;
   atoms=0;
  bool hasextended=false;
  getline(hfile,testvel);
  if(Tokenize(testvel,&ntokens)>6) hasvelocities=true;   
    else hasvelocities=false;
  if(Tokenize(testvel)>9) hasextended=true;
   else hasextended=false;
  // if(hasextended) cout << "Found extended Data!" << endl;
   hfile.seekg(0,ios::beg); fpos=0;
  string frame=ntokens[0];
  ntokens.clear();
   while(frame.compare(ntokens[0])==0)
   { getline(hfile,testvel);
       ntokens.clear();
     atoms++; 
     Tokenize(testvel,&ntokens);
   } 
   atoms--;
   hfile.clear();
   hfile.seekg(0,ios::beg); fpos=0;
 
   while ( !hfile.eof() )
   {
   getline(hfile, testvel);
   filelength++;
   }
   frames=framesfromlines(filelength,hasextended);
  hfile.clear();
  hfile.seekg(0,ios::beg); fpos=0;
  hasloaded=true;
  return true;    
}
