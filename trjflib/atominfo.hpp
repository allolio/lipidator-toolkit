//#include "boost/assign.hpp"
class Isotopeinfo
{
  public:
   Isotopeinfo(){AtomInit();}
   double GetIsoMass(string Type)
   {
       typedef   map<string,double>::iterator elementkey;
       elementkey elem=weighttable.find(Type);
       if(elem==weighttable.end()) { Type.resize(1); 
	 elem=weighttable.find(Type);}
       if(elem==weighttable.end()) 
	return 0.00; 
       return (*elem).second;

   }
   bool SetIsoMass(string Type, double mass)
   {
    typedef   map<string,double>::iterator elementkey;
    elementkey elem=weighttable.find(Type);
    if (elem!=weighttable.end()) (*elem).second=mass;
    else return false;
    return true;
   }
  private:
    typedef   map<string,double> typemass;
    typemass weighttable;
void AtomInit(void)
{
    weighttable = { 
                           make_pair("H",1.01),
                           make_pair("He",4.00),
                           make_pair("Li",6.94),
                           make_pair("Be",9.01),
                           make_pair("B",10.81),
                           make_pair("C",12.01),
                           make_pair("N",14.01),
                           make_pair("O",16.00),
                           make_pair("F",19.00),
                           make_pair("Ne",20.18),
                           make_pair("Na",22.99),
                           make_pair("Mg",24.31),
                           make_pair("Al",26.98),
                           make_pair("Si",28.09),
                           make_pair("P",30.98),
                           make_pair("S",32.06),
                           make_pair("Cl",35.45),
                           make_pair("Ar",39.95),
                           make_pair("K",39.10),
                           make_pair("Ca",40.08),
                           make_pair("Sc",44.96),
                           make_pair("Ti",47.90),
                           make_pair("V",50.94),
                           make_pair("Cr",52.00),
                           make_pair("Mn",54.94),
                           make_pair("Fe",55.85),
                           make_pair("Co",58.93),
                           make_pair("Ni",58.71),
                           make_pair("Cu",63.54),
                           make_pair("Zn",65.37),
                           make_pair("Ga",69.72),
                           make_pair("Ge",72.59),
                           make_pair("As",74.99),
                           make_pair("Se",78.96),
                           make_pair("Br",79.91),
                           make_pair("Kr",83.80),
                           make_pair("Rb",85.47),
                           make_pair("Sr",87.62),
                           make_pair("Y",88.91),
                           make_pair("Zr",91.22),
                           make_pair("Nb",92.91),
                           make_pair("Mo",95.94),
                           make_pair("Tc",96.91),
                           make_pair("Ru",101.07),
                           make_pair("Rh",102.90),
                           make_pair("Pd",106.40),
                           make_pair("Ag",107.87),
                           make_pair("Cd",112.40),
                           make_pair("In",114.82),
                           make_pair("Sn",118.69),
                           make_pair("Sb",121.75),
                           make_pair("Te",127.60),
                           make_pair("I",126.90),
                           make_pair("Xe",131.30),
                           make_pair("Cs",132.90),
                           make_pair("Ba",137.34),
                           make_pair("La",138.91),
                           make_pair("Ce",140.12),
                           make_pair("Pr",140.91),
                           make_pair("Nd",144.24),
                           make_pair("Pm",144.91),
                           make_pair("Sm",150.35),
                           make_pair("Eu",151.96),
                           make_pair("Gd",157.25),
                           make_pair("Tb",158.92),
                           make_pair("Dy",162.50),
                           make_pair("Ho",164.93),
                           make_pair("Er",167.26),
                           make_pair("Tm",168.93),
                           make_pair("Yb",173.04),
                           make_pair("Lu",174.97),
                           make_pair("Hf",178.49),
                           make_pair("Ta",180.95),
                           make_pair("W",183.85),
                           make_pair("Re",186.20),
                           make_pair("Os",190.20),
                           make_pair("Ir",192.22),
                           make_pair("Pt",195.09),
                           make_pair("Au",196.97),
                           make_pair("Hg",200.59),
                           make_pair("Tl",204.37),
                           make_pair("Pb",207.19),
                           make_pair("Bi",208.98),
                           make_pair("Po",208.98),
                           make_pair("At",209.99),
                           make_pair("Rn",222.02),
                           make_pair("Fr",223.02),
                           make_pair("Ra",226.00),
                           make_pair("Ac",227.03),
                           make_pair("Th",232.04),
                           make_pair("Pa",231.04),
                           make_pair("U",238.03),
                           make_pair("Np",237.00),
                           make_pair("Pu",242.00),
                           make_pair("Am",243.06),
                           make_pair("Cm",247.07),
                           make_pair("Bk",247.07),
                           make_pair("Cf",251.08),
                           make_pair("Es",254.09),
                           make_pair("Fm",257.10),
                           make_pair("Md",257.10),
                           make_pair("No",255.09),
                           make_pair("Lr",256.10)};
			   // Forcefield Type;
}

};

class Chargeinfo
{
  public:
   Chargeinfo(){AtomInit();}
   double GetAtomCharge(string Type)
   {
       typedef   map<string,double>::iterator elementkey;
       elementkey elem=chargetable.find(Type);
       return (*elem).second;

   }
   bool SetAtomCharge(string Type, double charge)
   {
    typedef   map<string,double>::iterator elementkey;
    elementkey elem=chargetable.find(Type);
    if (elem!=chargetable.end()) (*elem).second=charge;
    else return false;
    return true;
   }
  private:
    typedef   map<string,double> typecharge;
    typecharge chargetable;
void AtomInit(void)
{
    chargetable ={ 
                           make_pair("H",1),
			    make_pair("X",-2), // Wannier Center
                           make_pair("He",2.00),
                           make_pair("Li",3),
                           make_pair("Be",4),
                           make_pair("B",5),
                           make_pair("C",6),
                           make_pair("N",7),
                           make_pair("O",8),
                           make_pair("F",9),
                           make_pair("Ne",10),
                           make_pair("Na",11),
                           make_pair("Mg",12),
                           make_pair("Al",13),
                           make_pair("Si",14),
                           make_pair("P",15),
                           make_pair("S",16),
                           make_pair("Cl",17),
                           make_pair("Ar",18),
                           make_pair("K",19),
                           make_pair("Ca",20)};
  /*                         make_pair("Sc",44.96),
                           make_pair("Ti",47.90),
                           make_pair("V",50.94),
                           make_pair("Cr",52.00),
                           make_pair("Mn",54.94),
                           make_pair("Fe",55.85),
                           make_pair("Co",58.93),
                           make_pair("Ni",58.71),
                           make_pair("Cu",63.54),
                           make_pair("Zn",65.37),
                           make_pair("Ga",69.72),
                           make_pair("Ge",72.59),
                           make_pair("As",74.99),
                           make_pair("Se",78.96),
                           make_pair("Br",79.91),
                           make_pair("Kr",83.80),
                           make_pair("Rb",85.47),
                           make_pair("Sr",87.62),
                           make_pair("Y",88.91),
                           make_pair("Zr",91.22),
                           make_pair("Nb",92.91),
                           make_pair("Mo",95.94),
                           make_pair("Tc",96.91),
                           make_pair("Ru",101.07),
                           make_pair("Rh",102.90),
                           make_pair("Pd",106.40),
                           make_pair("Ag",107.87),
                           make_pair("Cd",112.40),
                           make_pair("In",114.82),
                           make_pair("Sn",118.69),
                           make_pair("Sb",121.75),
                           make_pair("Te",127.60),
                           make_pair("I",126.90),
                           make_pair("Xe",131.30),
                           make_pair("Cs",132.90),
                           make_pair("Ba",137.34),
                           make_pair("La",138.91),
                           make_pair("Ce",140.12),
                           make_pair("Pr",140.91),
                           make_pair("Nd",144.24),
                           make_pair("Pm",144.91),
                           make_pair("Sm",150.35),
                           make_pair("Eu",151.96),
                           make_pair("Gd",157.25),
                           make_pair("Tb",158.92),
                           make_pair("Dy",162.50),
                           make_pair("Ho",164.93),
                           make_pair("Er",167.26),
                           make_pair("Tm",168.93),
                           make_pair("Yb",173.04),
                           make_pair("Lu",174.97),
                           make_pair("Hf",178.49),
                           make_pair("Ta",180.95),
                           make_pair("W",183.85),
                           make_pair("Re",186.20),
                           make_pair("Os",190.20),
                           make_pair("Ir",192.22),
                           make_pair("Pt",195.09),
                           make_pair("Au",196.97),
                           make_pair("Hg",200.59),
                           make_pair("Tl",204.37),
                           make_pair("Pb",207.19),
                           make_pair("Bi",208.98),
                           make_pair("Po",208.98),
                           make_pair("At",209.99),
                           make_pair("Rn",222.02),
                           make_pair("Fr",223.02),
                           make_pair("Ra",226.00),
                           make_pair("Ac",227.03),
                           make_pair("Th",232.04),
                           make_pair("Pa",231.04),
                           make_pair("U",238.03),
                           make_pair("Np",237.00),
                           make_pair("Pu",242.00),
                           make_pair("Am",243.06),
                           make_pair("Cm",247.07),
                           make_pair("Bk",247.07),
                           make_pair("Cf",251.08),
                           make_pair("Es",254.09),
                           make_pair("Fm",257.10),
                           make_pair("Md",257.10),
                           make_pair("No",255.09),
                           make_pair("Lr",256.10)}*/
}

};


class Bondinfo
{
public:
 Bondinfo(void)
  {
    BondInit();
  }
   double GetBondMax(string TypeA,string TypeB)
   {  
       const double bdefault=1.6;
       typedef   map<string,double>::iterator bondkey;
       bondkey blen=bondtable.find(TypeA+"."+TypeB);
       if(blen==bondtable.end()) return bdefault;
	 else return (*blen).second;
   }
   void SetBondMax(string TypeA, string TypeB, double coff)
   {
       typedef map<string,double>::iterator bondkey;
       bondkey blen=bondtable.find(TypeA+"."+TypeB);
       if(blen==bondtable.end()) { 
	  bondtable.insert(make_pair(TypeA+"."+TypeB,coff));
	  bondtable.insert(make_pair(TypeB+"."+TypeA,coff));
          }
          else
	  {
	    (*blen).second=coff;
	    blen=bondtable.find(TypeB+"."+TypeA);
            (*blen).second=coff;  	    
	  }
   }
private:
  typedef   map<string,double> typebon;
    typebon bondtable;
    
  void BondInit(void){
    bondtable = { 
                           make_pair("H.O",1.2),
			   make_pair("O.H",1.2),
			   make_pair("H.N",1.3),
			   make_pair("N.H",1.3),
			   make_pair("1.8",1.2),
			   make_pair("8.1",1.2),
                           make_pair("Si.O",1.9),
			   make_pair("H.H",1.2),
			   make_pair("O.Si",1.9),
			   make_pair("X.H",1.1),
			   make_pair("X.O",1.2),
		           make_pair("X.X",1.1),
			   make_pair("XX.XX",0.001),
			   make_pair("XX.H",0.001),
			   make_pair("XX.C",0.001),
		           make_pair("XX.N",0.001),
			   make_pair("XX.S",0.001),
                           make_pair("XX.O",0.001),
                          make_pair("C.Cl",1.98),
                           make_pair("C.O",1.65),
                            make_pair("O.C",1.65),
                             make_pair("C.N",1.7),
                            make_pair("N.C",1.7),
		           make_pair("Cl.C",1.98),
                           make_pair("H.Cl",1.4),
		           make_pair("Cl.H",1.4),
                           make_pair("C.C",1.7),
			   make_pair("C.H",1.6),
                           make_pair("H.C",1.6),
                           make_pair("Au.Au",3.0),
                           make_pair("S.S",3.5),
                           make_pair("Au.S", 2.4),
                           make_pair("S.Au",2.4)};
	//		   ("H.",1.9)
  //   	                  ("X.N",0.001)
	//		   ("X.C",0.001)
		       //    ("H.X",1.1);
	//		   ("O.X",0.001)
	//	           ("C.X",0.001)
	//		   ("N.X",0.001);
			   // Add relevant bond lengths here
  }

};
