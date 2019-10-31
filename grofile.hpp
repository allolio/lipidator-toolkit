#ifndef HAVE_CONFIG
#include "config.hpp"
#endif
#include <stdarg.h>
//#include "topology.hpp"
//#define DEBUG

std::string string_format(const char *fmt, ...) {
std::vector<char> str(100);
va_list ap;
while (1) {
va_start(ap, fmt);
int n = vsnprintf(&str[0], str.size(), fmt, ap);
va_end(ap);
if (n > -1 && n < str.size()) {
str.resize(n);
return &str[0];
}
str.resize(str.size() * 2);
}
}

bool WriteGro(ostream &out,SnapShot &Snap, Topology &Top);

class GroFile : public TrajecFile
{   
  public:
    GroFile(string filename) {
    hasloaded=false;frames=0;fpos=0;lsnaps.clear(); lpos=lsnaps.begin(); 
    readcache=false;
    Load(filename);
    }
     GroFile() {
    hasloaded=false;frames=0;fpos=0;lsnaps.clear(); lpos=lsnaps.begin(); 
    readcache=false;fsize=0;
    }
    ~GroFile()
    {
    if(hasloaded)
    {
    hfile.close(); }
    }
     bool Load(string filename);
     bool SaveTo(string Filename);
     Topology &GetTopology();
     bool SetTopology (Topology &);
  protected:
     bool ParseNextFrame();
     bool GotoFrame(int n);
     int getmolecules(void);
     int framesfromlines(int filelength, bool hasextended);

     off64_t fsize;
     Topology top;
};


bool GroFile::ParseNextFrame()
{
string buffer;
if(readcache)  { GotoFrame(fpos);
#ifdef DEBUG
cout << "Last Frame from Cache!" << endl;
#endif
readcache=false;}
// Ignore First three lines
getline(hfile,buffer);
getline(hfile,buffer);
int htokens;

if(fpos>frames) return false;
if(!hfile.is_open()) return false;
#ifdef DEBUG
 cout << "parsing new frame" << endl;
#endif

vector<Atom> vatoms;
vatoms.resize(atoms+3); // 3 PBC Vectors

vector<string> filelines; filelines.resize(atoms);
for(int i=0;i!=atoms;i++)   getline(hfile,filelines[i]); 
#ifdef DEBUG
 cout << "read frame" << endl;
#endif
	#pragma omp parallel for
for(int i=0;i<atoms;i++)
  {
    double cfac =0.00045710285;
    bool flipflag=false;
   vector<string> line;
   htokens=Tokenize(filelines[i],&line, " \t");
  //    cout << htokens << i <<  " "<< line[4] <<   " "<< line[3] << " "<< line[2] << endl;
   vatoms[i].type=top.molecules[top.hastype[i]].atoms[top.moltypes[i]];   
   if(i>=9999) flipflag=true;
      if(htokens>=6) {  
   vatoms[i].pos.x=atof(line[3-flipflag].c_str())*10;
   vatoms[i].pos.y=atof(line[4-flipflag].c_str())*10;
   vatoms[i].pos.z=atof(line[5-flipflag].c_str())*10;
   }
   if(hasvelocities && htokens > 6)
   {
    // a/
   vatoms[i].vel.x=atof(line[6-flipflag].c_str())*cfac;
   vatoms[i].vel.y=atof(line[7-flipflag].c_str())*cfac;
   vatoms[i].vel.z=atof(line[8-flipflag].c_str())*cfac;
   }
   line.clear();
  }
     vector<string> line;
  #ifdef DEBUG
 cout << "ok" << endl;
#endif  
  // Get PBC
  getline(hfile,filelines[0]);
  htokens=Tokenize(filelines[0],&line, " \t");
   for(int i=0;i!=3;i++)
  {
	vatoms[atoms+i].type="PBC";
	vatoms[atoms+i].pos=vatoms[atoms+i].pos*0;  

	if(i==0) vatoms[atoms+i].pos.x=atof(line[0].c_str());
        if(i==1) vatoms[atoms+i].pos.y=atof(line[1].c_str());
	if(i==2) vatoms[atoms+i].pos.z=atof(line[2].c_str());
  }
  #ifdef DEBUG
 cout << "ok" << endl;
#endif  
 SnapShot snap;
 lsnaps.insert(lpos,make_pair<SnapShot,int>( snap,fpos));
   lpos--; 
   ((*lpos).first).SnapOverwrite(vatoms,hasvelocities,true)   
;
   fpos++;

 return true;
}

bool GroFile::GotoFrame(int n)
{ string buffer;
  if (n>frames ) return false;
  if(fpos>n && hasloaded) {  hfile.clear();
    hfile.seekg(0,ios::beg); fpos=0;  }
    else if(!hasloaded) return false;
  while(fpos<n)
  {
    for(int i=0; i!=(atoms+3);i++) getline(hfile,buffer);  

    fpos++;
  } 
  return true;
}

bool GroFile::SaveTo(string filename)
{
   ofstream out;
  out.open(filename.c_str());
  GotoFrame(0); SnapShot snap;
  while(true)   

  {
    snap=GetNextSnap();
    if(snap.Valid()) WriteGro(out,snap,top);
    else {break;}
  }
 out.close();
return true;
}


bool GroFile::Load(string filename) // File Must have 3 Lines!
{
  int filelength=0;
  string testvel;
  vector<string> ntokens;
  hfile.open(filename.c_str());
  if(hfile.fail()) return false;   

  getline(hfile,comment);

  hfile >> atoms; // Only Valid for XYZ beware
  bool hasextended=false;
  getline(hfile,testvel);
  getline(hfile,testvel);
  
//  if(strlen(testvel.c_str())==0) getline(hfile,testvel);
//  cout << testvel << " " << comment << endl;
  if(Tokenize(testvel, NULL, " \t")>9) hasvelocities=true;
    else hasvelocities=false;  
      #ifdef DEBUG
      if(hasvelocities) cout << "Found Velocities" << endl;
      else cout << "Cant Find Velocities " <<  testvel << endl;
      #endif
         hfile.seekg(0,ios::beg); fpos=0;
 getline(hfile,testvel);
  getline(hfile,testvel);
   getmolecules();
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

int GroFile::getmolecules(void)
{
 vector<string> filelines; filelines.resize(atoms);
for(int i=0;i!=atoms;i++)   getline(hfile,filelines[i]); 
vector<string> molname; molname.resize(atoms); 
vector<string> atomname; atomname.resize(atoms);
// Parse File
     #ifdef DEBUG
 cout << "Reading Topology " << endl;
    #endif
 top.atmol.resize(atoms);
 top.hastype.resize(atoms);
 top.moltypes.resize(atoms);
 
	#pragma omp parallel for
for(int i=0;i<atoms;i++)
  {
   vector<string> line;   
   Tokenize(filelines[i],&line, " \t");
 //  if(htokens<3) { cout << "Broken Frame ! Cannot Parse" << endl; return false;}
   top.atmol[i]=atoi(line[0].c_str());
   int digits=log10(top.atmol[i])+1;
   int slen=(line[0].length())-digits;
   molname[i].resize(slen);
   line[0].copy((char*) molname[i].c_str(),slen,digits);
   atomname[i]=line[1];
   if(i>=9999) atomname[i].resize(atomname[i].length()-5);
   line.clear();
  }
      #ifdef DEBUG
 cout << "Read Topology " << endl;
    #endif
  // Create Topology
     int typenr=0;
     int lastmol=top.atmol[0];
  for(int i=0;i<atoms;i++)
  {
    top.hastype[i]=-1;
    for(int j=0;j!=top.molecules.size();j++)
     if(molname[i]==top.molecules[j].name) {top.hastype[i]=j;}
     if(top.atmol[i]!=lastmol) {typenr=0;lastmol=top.atmol[i];}
    top.moltypes[i]=typenr;

    if(top.hastype[i]==-1)   
    {
      #ifdef DEBUG
 cout << "New Molecule Type " << molname[i] << endl;
    #endif
      top.hastype[i]=top.molecules.size();
      if(top.molecules.size()==0) top.hastype[i]=0;
      Molecule mol;
      top.moltypes[i]=0;
      mol.name=molname[i];
      mol.atoms.push_back(atomname[i]);
      typenr=1;
      while(top.atmol[i+1]==top.atmol[i])
      {
       i++; 
       top.hastype[i]=top.molecules.size();
       mol.atoms.push_back(atomname[i]);
       top.moltypes[i]=typenr;
       typenr++;
      }     
      top.molecules.push_back(mol); 
            #ifdef DEBUG
     cout << " Has " << typenr+1 << "Atoms" << mol.atoms.size() << endl;
    #endif

      typenr=-1;
    }
    typenr++;
  }
 top.nmol=top.atmol[top.atmol.size()-1];
    #ifdef DEBUG
 cout << "Built Topology " << endl;
    #endif
 return top.nmol;  
 
 
}

Topology &GroFile::GetTopology()
{
  return top;
}

bool GroFile::SetTopology(Topology &topol)
{
  if(topol.atmol.size()!=top.atmol.size()) {cerr << "Topology Invalid" << endl; return false;}
  top=topol;
  return true;
}

int GroFile::framesfromlines(int filelength, bool hasextended)
{
 frames=filelength/(atoms+3);
  if((filelength%(atoms+3))!=0) {
      #ifdef DEBUG
      cout << "Broken Trajectory" << endl;
      cout << filelength << " and " << atoms << endl;
      #endif
    //     frames--;
      }
            return frames;

}

bool WriteGro(ostream &out,SnapShot &Snap, Topology &Top)
{
      double cfac =1/0.00045710285;

  int natoms=Snap.GetAtoms();
  out << "Generated by XYZlib  " << endl;
  out << natoms << endl;
  int nmol=Top.nmol;
      out.precision(3);
      out.setf(ios::fixed,ios::floatfield);  
  int running=1;
          int counter=1;

  for(int j=0;j!=Top.molecules.size();j++)
  {
    vector<int> mols=Top.GetMolsOfType(j);
    cout << "Mols" << mols.size() << endl;
    for(int k=0;k!=mols.size();k++)
    {
    int atomnum=0;
    while(atomnum<Top.molecules[j].atoms.size())
    {
    for(int i=0;i!=natoms;i++)
     {
      if (mols[k]==Top.atmol[i] && Top.moltypes[atomnum]==Top.moltypes[i])
      {

      Atom &a=Snap.GetAtom(i);
      //
//	out << " "<< counter << Top.molecules[j].name <<%8.4f%8.4f%8.4f " " << Top.molecules[j].atoms[atomnum] << " " << running  << " "<< a.pos.x << "  "  << a.pos.y << " " << a.pos.z;
      out << string_format("%5d%-5s%5s%5d%8.3f%8.3f%8.3f", counter, Top.molecules[j].name.c_str() , Top.molecules[j].atoms[atomnum].c_str(), running , a.pos.x/10 , a.pos.y/10 , a.pos.z/10);
	if(Snap.hasvel()) out <<  string_format("%8.4f%8.4f%8.4f", a.vel.x*cfac , a.vel.y*cfac , a.vel.z*cfac);
	out << endl;
	running++;
     
	atomnum++;
   //   cout << i << " @MOL" << Top.atmol[i] << " Atomtype: " << Top.moltypes[i] << "Moltype"<<  Top.hastype[i] <<  endl;
      }
     }
    }
       counter++;  
    
    }
  }
  // Add PBC
        out.precision(5);

  out << Snap.GetAtom(natoms).pos.x << " " <<   Snap.GetAtom(natoms+1).pos.y << " " << Snap.GetAtom(natoms+2).pos.z << endl;
return true;
  
}


