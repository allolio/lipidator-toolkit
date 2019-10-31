bool MkCPMDAtoms(SnapShot &s, ostream &out, string ecpdir="../ecp")
{
  typedef   map<string,string> inputmap;
  typedef   map<string,string>::iterator inputkey;
  
  inputmap cpmddict;
  inputkey pspkey;
  cpmddict.insert(pair<string,string>("C","C_MT_BLYP KLEINMAN BYLANDER\nLMAX=P"));
  cpmddict.insert(pair<string,string>("O","O_MT_BLYP KLEINMAN BYLANDER\nLMAX=P"));
  cpmddict.insert(pair<string,string>("N","N_MT_BLYP KLEINMAN BYLANDER\nLMAX=P"));
  cpmddict.insert(pair<string,string>("H","H_MT_BLYP \nLMAX=S"));
  vector<string> Atomtypes;
  vector<int> Atomindex;
  Atomtypes=s.GetUniqueAtomTypes();
  for(int i=0;i!=Atomtypes.size();i++)
  {
  pspkey=cpmddict.find(Atomtypes[i]);
  out << "*" << ecpdir << "/" <<(*pspkey).second << endl;
  Atomindex=s.GetIndicesOfType(Atomtypes[i]);
  s.IndextoStream(Atomindex,out);
  }
  
  return true;
}

bool ForceMap( SnapShot s1, SnapShot s2, double cf=1, string fname="vmd")
{   double sumsq=0,hival=0,loval=0;
    ofstream theforce;
    theforce.open((fname+ ".state").c_str());
    cout << fname << endl;
        cout.precision(6);
        //double cf=0.529177249;
          
    for(int i=0;i!=s1.GetAtoms();i++)
	{
	 
	Atom a=s1.GetAtom(i); triple tup;
	Atom b=s2.GetAtom(i); 
        tup=a.vel+b.vel*cf;
        //cout << tup.x << " " << tup.y << " " << tup.z << "  " << a.type << "  " << a.vel.x << " " << b.vel.x*cf << endl;
        if(hival<tup.abs()) hival=tup.abs();
	if(loval>tup.abs()) loval=tup.abs();
	}
	
	for(int i=0;i!=s1.GetAtoms();i++)
	{
	 
	Atom a=s1.GetAtom(i); triple tup;
	Atom b=s2.GetAtom(i); 
        tup=a.vel+b.vel*cf;
	double thev= 15-((tup.abs()+loval)/(hival+loval)*15);
	theforce << lround(thev) << " " << i  << endl;
       }
	//sumsq/=s1.GetAtoms();
        sumsq=sqrt(sumsq);
	//suma=sqrt(suma);
	//sumb=sqrt(sumb)*sqrt(cf*cf);
	theforce.close();
    //    cout << " ||delta|| : " << sumsq << " Max Dev Abs " << hival << " ||A|| " << suma << " || B || " << sumb << " [% of ||A|| :" << sumsq/suma*100 << "]  [% of ||B|| :" << sumsq/sumb*100  << "] fac " << cf << endl;
	return true;
}

bool CalcRMS( SnapShot s1, SnapShot s2, double cf=1)
{
  double sumsq=0,hival=0,suma=0, sumb=0,carbopr;
        cout.precision(6);
        //double cf=0.529177249;
       
    for(int i=0;i!=s1.GetAtoms();i++)
	{
	triple zero;
        zero=zero*0;
     
	Atom a=s1.GetAtom(i); triple tup;
	Atom b=s2.GetAtom(i); 
        tup=a.vel+b.vel*cf;
        suma+=a.vel.abs()*a.vel.abs();
        sumb+=b.vel.abs()*b.vel.abs();
   //   cout << tup.x << " " << tup.y << " " << tup.z << "  " << a.type << "  " << a.vel.x << " " << b.vel.x*cf << endl;
        sumsq+=tup.abs()*tup.abs();
	if(i==10)  carbopr=sqrt(sumsq);
        if(hival<tup.abs()) hival=tup.abs();
	}
	//sumsq/=s1.GetAtoms();
        sumsq=sqrt(sumsq);
	suma=sqrt(suma);
	sumb=sqrt(sumb)*sqrt(cf*cf);
	
        cout << " ||delta|| : " << sumsq << " Max Dev Abs " << hival << " ||A|| " << suma << " || B || " << sumb << " [% of ||A|| :" << sumsq/suma*100 << "]  [% of ||B|| :" << sumsq/sumb*100  << "] fac " << cf << " ";
	cout << "%C :" << carbopr << " /  % of ||delta||:" << carbopr/sumsq*100 << " ";
	return true;
}

bool CalcInterRMS( SnapShot s1, SnapShot s2, double cf, vector<int> MoleculeAtoms)
{
  triple zero;
  zero=zero*0;  
  double sumsq=0,hival=0,sum=0,carbopr;
        cout.precision(6);
        //double cf=0.529177249;
	int frag =0,maxfrag=0,mol=0;
	 vector<int> molecule;
    for(int f=0;f!=s1.GetAtoms();f++)
       if(MoleculeAtoms[f]>=maxfrag) maxfrag=MoleculeAtoms[f]; 
    
    for(frag=1;frag!=maxfrag;frag++) // Do not take MQ==frag 0!
    {
    for(int f=0;f!=s1.GetAtoms();f++)
    {
     if(MoleculeAtoms[f]==frag) molecule.push_back(f); 
    }
    
    triple forcea=zero;
    triple forceb=zero;
  //  cout << "mol";
    for(int i=0;i!=molecule.size();i++)
	{
//	 cout << molecule[i] << " " ;
	 Atom a=s1.GetAtom(molecule[i]); 
         Atom b=s2.GetAtom(molecule[i]);
	 forcea=forcea+a.vel;
	 forceb=forceb+b.vel*cf;
	}
//	cout << endl; 
	triple diff=(forcea+forceb);
        sum=diff.abs();
	if(sum>hival) {hival=sum;mol=frag;}
	sumsq+=sum;
	molecule.clear();
	forcea=forcea*0;
	forceb=forceb*0;
    }
      cout << "Avg. Diff. Inter :" << sumsq/maxfrag << " Max. Diff. " << hival << " @ Mol " << mol << " of" << maxfrag << " ";
     #ifdef DEBUG
      cout << "Atom Indices of Molecule: "<< mol << endl;
      for(int f=0;f!=s1.GetAtoms();f++)
       if(MoleculeAtoms[f]==mol) cout << f << " " ;
       #endif
          
 	return true;
}
 triple boxdiff(triple a,double region)
    {   
        int indicator;
   if(a.x >= 0.5*region)
    {a.x -= region; indicator =1;}
  else if(a.x < -0.5*region)
    {a.x += region;indicator=1;}

  if(a.y >= 0.5*region)
    {a.y -= region; indicator =1;}
  else if(a.y < -0.5*region)
    {a.y += region;indicator=1;}

  if(a.z >= 0.5*region)
    {a.z -= region; indicator =1;}
  else if(a.z < -0.5*region)
    {a.z += region;indicator=1;}

	return a;
    }
   
bool CalcInterRMSCutoff( SnapShot s1, SnapShot s2, double cf, vector<int> MoleculeAtoms, int around, double cutoff, triple box)
{ // We calculate the Average Delta of Molecules inside cutoff except the one we are in!
  // box is for PBC
   triple zero;
  zero=zero*0;  
 
  double sumsq=0,hival=0,sum=0,carbopr;
        cout.precision(6);
        //double cf=0.529177249;
	int frag =0,maxfrag=0,mol=0,usedmol=0;
	 vector<int> molecule;
    for(int f=0;f!=s1.GetAtoms();f++)
    if(MoleculeAtoms[f]>=maxfrag) maxfrag=MoleculeAtoms[f]; 
    int fragofa = MoleculeAtoms[around]; // Get fragment that is being developed around
//    cout << "Atom in Fragment" << fragofa << endl;
    for(frag=0;frag!=maxfrag;frag++)
    {
//     cout << "Frag of A " << fragofa << endl;
    if(frag!=fragofa) // No interest in my own molecule
    {// Get Atoms in Fragment
    for(int f=0;f!=s1.GetAtoms();f++)
    {
     if(MoleculeAtoms[f]==frag) molecule.push_back(f); 
    }
    // Check if inside cutoff
     bool evmol=false; 
    Atom center=s1.GetAtom(around);
    triple dcenter=zero;
    triple forcea=zero;
    triple forceb=zero;
  
  //  cout << "mol";
    for(int i=0;i!=molecule.size();i++)
	{ Atom a=s1.GetAtom(molecule[i]); 
	 // PBC
	  a.pos.x=fmod(a.pos.x,box.x);
	 a.pos.y=fmod(a.pos.y,box.y);
	 a.pos.z=fmod(a.pos.z,box.z);
	
	   dcenter=center.pos-a.pos;
           dcenter=boxdiff(dcenter,box.x);
          
	// cout << molecule[i] << " " << dcenter.abs();
	  if(dcenter.abs()<cutoff){
	    
	    evmol=true;
	  }
	 /* ivec=boxdiff(a.pos,box.x);
	  triple jvec=boxdiff(center.pos,box.x);
	  dcenter=ivec-jvec;*/
	  
        }
   // Calculate RMS
   	if(evmol){
	  usedmol++; cout << endl;
    for(int i=0;i!=molecule.size();i++)
	{
	  cout << frag << ":" <<molecule[i] << endl;
	 Atom a=s1.GetAtom(molecule[i]); 
         Atom b=s2.GetAtom(molecule[i]);
	 forcea=forcea+a.vel;
	 forceb=forceb+b.vel*cf;
	 // trivial
  //        if(evmol) cout << molecule[i] << " "; 

	}

//	cout << endl; 
	triple diff=(forcea+forceb);
        sum=diff.abs();
	if(sum>hival) {hival=sum;mol=frag;}
	sumsq+=sum;
	forcea=forcea*0;
	forceb=forceb*0;}
        molecule.clear();
    }
    }
      cout << "Av Diff @: " << around << " ctf " << cutoff << "  " << "mols " << usedmol << "  "<< sumsq/usedmol << " Maximum Difference " << hival << " @ Molecule " << mol << " ";
    // #ifdef DEBUG
      cout << "Atom Indices of Molecule: "<< mol << endl;
      for(int f=0;f!=s1.GetAtoms();f++)
       if(MoleculeAtoms[f]==mol) cout << f << " " ;
    //   #endif
          
 	return true;
}