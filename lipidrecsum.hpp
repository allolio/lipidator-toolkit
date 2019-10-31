#ifndef HAVE_CONFIG
#include "config.hpp"
#endif
#define READ_FAIL -1
#define WRITE_FAIL -1
#include "ndxtool.hpp"
using namespace std;
typedef vector<int> atomlist;

//1 2  ---  3 3 3
struct node
{
atomlist indices;
vector<double> masses;
int center;
string name;
};

typedef vector<node> nodes;


class LipidTopology
{
public:
  int ReadTopology(ifstream &ndx);
  int WriteTopology(ostream &ndx);
  nodes* GetMap(void);
  bool SetMap(nodes &nod);
  string GetName();
  void SetName(string name);
private:
  nodes thenodes;
  string lipidname;
};

class LipidDirectors : public Analyzer
{
public:
    LipidDirectors(Trajectory* t): Analyzer(t)
    {
      
    }
  bool LoadLibrary(string filename);
  bool LoadIndices(string ndxfile);
  bool Evaluate(SnapShot *current=NULL);
  bool GetLipidEnv(int index, vector<triple> &lpos, atomlist &indices,double cutoff,bool strict);
  vector<triple> *GetDirectors(void);
  vector<triple> *GetCenters(void);
  void PrintNodes();
  
private:
  vector<LipidTopology> thelipids;
  vector<triple> thedirectors;
  vector< node>  indexnodes;
  vector<triple> centers;
};

void LipidDirectors::PrintNodes()
{
  cout << "N NODES: " << indexnodes.size() << endl;
  for(int i=0;i!=indexnodes.size();i++)
  {
    node n=indexnodes[i];
    cout << n.name << endl;
    cout << n.center << endl;
    for(int j=0;j!=n.indices.size();j++)
    cout << n.indices[j] << " " ;
    cout << endl;
  }
}

bool LipidDirectors::LoadLibrary(string filename)
{
  ifstream lib;
  lib.open(filename.c_str());
  bool running=true;
  while(running)
  {
  LipidTopology top;
  int result=top.ReadTopology(lib);
  //cout << "READ "<< result << "GROUPS" <<endl;
  if(result < 0) break; 
  thelipids.push_back(top);
  }
  //cout << "READ" << thelipids.size() << "LIPIDS" << endl;
  lib.close();
  return true;
}

bool LipidDirectors::LoadIndices(string ndxfile)
{

  if(thelipids.size()==0) {cerr<< "Load Library First" << endl; return false;}
  int globalcount=0;
//  cout << thelipids.size() << "LIPIDTYPES" << endl;;
  for(int i=0;i!=thelipids.size();i++)
    {
  //   cout << "READINDEX" << endl;
      atomlist offsets=NDX2Atomlist(ndxfile,thelipids[i].GetName());
   //  cout << "FOUND " << offsets.size() << endl;
	nodes *t=thelipids[i].GetMap();
//	cout << "LIPID_COL" << t->size() << " " <<endl;;      
      //if(offsets.size()==0) return false;
      for(int j=0;j!=offsets.size();j++)
      {
	for(int k=0;k!=t->size();k++)
	{  
	 node &cur=t->at(k);
	 node newnode;
//	 cout <<"IND" << cur.indices.size() << endl;;
	 for(int l=0;l!=cur.indices.size();l++)
	  {
	   newnode.indices.push_back(cur.indices[l]+offsets[j]);
	  }
	 newnode.masses=cur.masses;
	 newnode.center=cur.center;
	 newnode.name=SSTR(globalcount);
//	 newnode.name="MUPPETX";
//	  newnode.name=SSTR(globalcount);
//	 cout << globalcount << "LOCAL " << offsets.size() << endl;
	 indexnodes.push_back(newnode);
	}
	       globalcount++;

       }

    }
    thedirectors.resize(globalcount);
    centers.resize(globalcount);
//   cout << "ASSIGNED " << globalcount << " LIPIDS" <<endl;
  return true;
}

bool LipidDirectors::Evaluate(SnapShot* current)
{
  if(current==NULL)  current=&traj->Current(); 
  int globalcount=0;
  for( int i=0;i!=indexnodes.size();i++)
  {
    string id=indexnodes[i].name;
    vector<triple> pos1;
    vector<triple> pos2;
    triple center;
    while(id==indexnodes[i].name && i!=indexnodes.size())
    {
     // for(int k=0;k!=indexnodes[i].indices.size();k++) cout << indexnodes[i].indices[k] << endl;
      triple pos=GetIndWAvgPos(indexnodes[i].indices, indexnodes[i].masses,current,true);
  //    cout << "CENTER" << indexnodes[i].center << endl;
      switch (indexnodes[i].center)
      {
	case 1:
	pos1.push_back(pos);
	break;
	case 2:
	pos2.push_back(pos);
	case 3: 
	center=pos;
      }
      i++;
    }
    i--;
//    cout << "POS1" << pos1.size() << " POS2 " << pos2.size() << endl;
   triple o1=pos1[0],o2=pos2[0];

    for(int j=1;j!=pos1.size();j++)
    {
	pos1[0]=pos1[0]+o1+boxdist(pos1[j],o1);
    }
    for(int j=1;j!=pos2.size();j++)
	pos2[0]=pos2[0]+o2+boxdist(pos2[j],o2);
    pos1[0]=pos1[0]/pos1.size();
    pos2[0]=pos2[0]/pos2.size();
  
//    thedirectors[globalcount]=boxdist(pos2[0],pos1[0]);
        thedirectors[globalcount]=boxdist(pos1[0],pos2[0]);

  //  thedirectors[globalcount].print();
    thedirectors[globalcount]=thedirectors[globalcount]/thedirectors[globalcount].abs();
    centers[globalcount]=center;
    pos1.clear();pos2.clear();
    globalcount++;
  }
  return true;
}

vector< triple >* LipidDirectors::GetDirectors(void)
{
 if(thedirectors.size()==0) return NULL;
 return &thedirectors;
}

vector< triple >* LipidDirectors::GetCenters(void)
{
 if(centers.size()==0) return NULL;
 return &centers;
}

bool LipidDirectors::GetLipidEnv(int index, vector< triple >& lpos, atomlist& indices, double cutoff=25, bool strict=false)
{// CUTOFF 15
  list<pair<double,int> > memory(0);
  double worst=0;
  int sz=lpos.size();
  for(int i=0;i!=centers.size();i++)
  {
    double diff=boxdist(centers[index],centers[i]).abs();
    if(diff < cutoff)
    {
      if(memory.size() < 2*sz || diff < worst)
      {
      if(diff> worst) worst=diff;
      memory.push_back(make_pair(diff,i));
      }
    }
  }

  if(memory.size()<sz && strict) return false;
  memory.sort();
  if(memory.size()>3*sz) memory.resize(sz*3);
  int oldind=-1000;
  for (list< pair <double, int> >::iterator it = memory.begin(); it != memory.end(); it++)
    {
          int ind=(*it).second;
      double sprod=(thedirectors[index]*thedirectors[ind]);
      if(sprod < 0.1) {
	      list< pair <double, int> >::iterator kill=it; it--;
	      memory.erase(kill);}
    }
 //   cout << "memory " << memory.size();

  if(memory.size()<sz && strict || memory.size() < 6) {//cout << "FAIL" << endl ;
    return false;}
if(strict)  memory.resize(lpos.size());
  else lpos.resize(memory.size());
  int i=0;
  indices.resize(lpos.size());
  for (list< pair <double, int> >::iterator it = memory.begin(); it != memory.end(); it++)
    {
    int ind=(*it).second;  
//      if(ind==oldind) cerr << "DOUBLETTE WTF!" << endl;
      
      indices[i]=ind;
      lpos[i]=boxdist(centers[index],centers[(*it).second]); 
    //  lpos[i].z=lpos[i].z*0.0001;
      i++;
      //;i++;
//      centers[(*it).second].print();
 //     cout << lpos[i].abs() << " " << index << " " << (*it).second << endl ; i++;
          oldind=ind;
  
    }

  return true;
}


int LipidTopology::ReadTopology(ifstream &ndx)
{
  Isotopeinfo iso;
  string buffer;
  int result=READ_FAIL;
  vector<string> ltokens;
  while(ltokens.size()==0)
  {  if(!getline(ndx,buffer)) return -1;
     Tokenize(buffer,&ltokens,"#");
//     cout << buffer << endl;
  }
  //if(ltokens.size()==0) return -1;
  lipidname=ltokens[0];
   while (getline(ndx,buffer) && ltokens.size() >0)
  { ltokens.clear(); 
    Tokenize(buffer,&ltokens,"[]");
    while(ltokens.size()>0)     
    {
//    cout << "NODE" << endl;
      node current;
      current.name=ltokens[0];
      if(!isdigit( *((ltokens[1]).c_str()))) return READ_FAIL;
      current.center=atoi(ltokens[1].c_str());
      while (getline(ndx,buffer))
      {
      vector<string> tokens;
      Tokenize(buffer,&tokens);
      if(tokens.size()==0) break;
      if(buffer[0]=='#' || buffer[0]=='[') break;
      if(!isalnum( *((tokens[0]).c_str()))) break;
	for(int i=0;i!=tokens.size()/2;i++)
	  {
//	    cout << buffer[0] << " " << tokens[i] << endl;
	    current.indices.push_back(atoi(tokens[2*i].c_str())-1);
	    double mass=1;
	    if(isalpha(*(tokens[2*i+1].c_str()))) mass=iso.GetIsoMass(tokens[2*i+1]);
	    else mass=atof(tokens[2*i+1].c_str());
	 //   cout << "INDEX" << tokens[2*i] << endl;
	    current.masses.push_back(mass);
	  }
      }
      thenodes.push_back(current);
      ltokens.clear();
      Tokenize(buffer,&ltokens,"[]");
      if(buffer[0]=='#') break;
     
    }
   }
   if (thenodes.size()==0) return READ_FAIL;
   return  thenodes.size();
}

int LipidTopology::WriteTopology(ostream &ndx)
{
  int result=WRITE_FAIL;
  ndx << "#"<< lipidname << endl;
  for (int i=0;i!=thenodes.size();i++)
  {
    node cur=(thenodes[i]);
    ndx << "[" << cur.name <<"][" << cur.center << "]" << endl;
    for(int j=0;j!=cur.indices.size();j++)
      ndx << cur.indices[j]+1 << " " ;
    ndx << endl;
  }
  return thenodes.size();
}

bool LipidTopology::SetMap(nodes& nod)
{
thenodes=nod;
return true;
}

nodes* LipidTopology::GetMap(void)
{
return &thenodes;
}

string LipidTopology::GetName()
{
return lipidname;
}

void LipidTopology::SetName(string name)
{
lipidname=name;
}
