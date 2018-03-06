#include "boost/assign.hpp"
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
    using namespace boost::assign;
    weighttable = map_list_of
                           ("H",1.01)
                           ("He",4.00)
                           ("Li",6.94)
                           ("Be",9.01)
                           ("B",10.81)
                           ("C",12.01)
                           ("N",14.01)
                           ("O",16.00)
                           ("F",19.00)
                           ("Ne",20.18)
                           ("Na",22.99)
                           ("Mg",24.31)
                           ("Al",26.98)
                           ("Si",28.09)
                           ("P",30.98)
                           ("S",32.06)
                           ("Cl",35.45)
                           ("Ar",39.95)
                           ("K",39.10)
                           ("Ca",40.08)
                           ("Sc",44.96)
                           ("Ti",47.90)
                           ("V",50.94)
                           ("Cr",52.00)
                           ("Mn",54.94)
                           ("Fe",55.85)
                           ("Co",58.93)
                           ("Ni",58.71)
                           ("Cu",63.54)
                           ("Zn",65.37)
                           ("Ga",69.72)
                           ("Ge",72.59)
                           ("As",74.99)
                           ("Se",78.96)
                           ("Br",79.91)
                           ("Kr",83.80)
                           ("Rb",85.47)
                           ("Sr",87.62)
                           ("Y",88.91)
                           ("Zr",91.22)
                           ("Nb",92.91)
                           ("Mo",95.94)
                           ("Tc",96.91)
                           ("Ru",101.07)
                           ("Rh",102.90)
                           ("Pd",106.40)
                           ("Ag",107.87)
                           ("Cd",112.40)
                           ("In",114.82)
                           ("Sn",118.69)
                           ("Sb",121.75)
                           ("Te",127.60)
                           ("I",126.90)
                           ("Xe",131.30)
                           ("Cs",132.90)
                           ("Ba",137.34)
                           ("La",138.91)
                           ("Ce",140.12)
                           ("Pr",140.91)
                           ("Nd",144.24)
                           ("Pm",144.91)
                           ("Sm",150.35)
                           ("Eu",151.96)
                           ("Gd",157.25)
                           ("Tb",158.92)
                           ("Dy",162.50)
                           ("Ho",164.93)
                           ("Er",167.26)
                           ("Tm",168.93)
                           ("Yb",173.04)
                           ("Lu",174.97)
                           ("Hf",178.49)
                           ("Ta",180.95)
                           ("W",183.85)
                           ("Re",186.20)
                           ("Os",190.20)
                           ("Ir",192.22)
                           ("Pt",195.09)
                           ("Au",196.97)
                           ("Hg",200.59)
                           ("Tl",204.37)
                           ("Pb",207.19)
                           ("Bi",208.98)
                           ("Po",208.98)
                           ("At",209.99)
                           ("Rn",222.02)
                           ("Fr",223.02)
                           ("Ra",226.00)
                           ("Ac",227.03)
                           ("Th",232.04)
                           ("Pa",231.04)
                           ("U",238.03)
                           ("Np",237.00)
                           ("Pu",242.00)
                           ("Am",243.06)
                           ("Cm",247.07)
                           ("Bk",247.07)
                           ("Cf",251.08)
                           ("Es",254.09)
                           ("Fm",257.10)
                           ("Md",257.10)
                           ("No",255.09)
                           ("Lr",256.10);
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
    using namespace boost::assign;
    chargetable = map_list_of
                           ("H",1)
			    ("X",-2) // Wannier Center
                           ("He",2.00)
                           ("Li",3)
                           ("Be",4)
                           ("B",5)
                           ("C",6)
                           ("N",7)
                           ("O",8)
                           ("F",9)
                           ("Ne",10)
                           ("Na",11)
                           ("Mg",12)
                           ("Al",13)
                           ("Si",14)
                           ("P",15)
                           ("S",16)
                           ("Cl",17)
                           ("Ar",18)
                           ("K",19)
                           ("Ca",20)
  /*                         ("Sc",44.96)
                           ("Ti",47.90)
                           ("V",50.94)
                           ("Cr",52.00)
                           ("Mn",54.94)
                           ("Fe",55.85)
                           ("Co",58.93)
                           ("Ni",58.71)
                           ("Cu",63.54)
                           ("Zn",65.37)
                           ("Ga",69.72)
                           ("Ge",72.59)
                           ("As",74.99)
                           ("Se",78.96)
                           ("Br",79.91)
                           ("Kr",83.80)
                           ("Rb",85.47)
                           ("Sr",87.62)
                           ("Y",88.91)
                           ("Zr",91.22)
                           ("Nb",92.91)
                           ("Mo",95.94)
                           ("Tc",96.91)
                           ("Ru",101.07)
                           ("Rh",102.90)
                           ("Pd",106.40)
                           ("Ag",107.87)
                           ("Cd",112.40)
                           ("In",114.82)
                           ("Sn",118.69)
                           ("Sb",121.75)
                           ("Te",127.60)
                           ("I",126.90)
                           ("Xe",131.30)
                           ("Cs",132.90)
                           ("Ba",137.34)
                           ("La",138.91)
                           ("Ce",140.12)
                           ("Pr",140.91)
                           ("Nd",144.24)
                           ("Pm",144.91)
                           ("Sm",150.35)
                           ("Eu",151.96)
                           ("Gd",157.25)
                           ("Tb",158.92)
                           ("Dy",162.50)
                           ("Ho",164.93)
                           ("Er",167.26)
                           ("Tm",168.93)
                           ("Yb",173.04)
                           ("Lu",174.97)
                           ("Hf",178.49)
                           ("Ta",180.95)
                           ("W",183.85)
                           ("Re",186.20)
                           ("Os",190.20)
                           ("Ir",192.22)
                           ("Pt",195.09)
                           ("Au",196.97)
                           ("Hg",200.59)
                           ("Tl",204.37)
                           ("Pb",207.19)
                           ("Bi",208.98)
                           ("Po",208.98)
                           ("At",209.99)
                           ("Rn",222.02)
                           ("Fr",223.02)
                           ("Ra",226.00)
                           ("Ac",227.03)
                           ("Th",232.04)
                           ("Pa",231.04)
                           ("U",238.03)
                           ("Np",237.00)
                           ("Pu",242.00)
                           ("Am",243.06)
                           ("Cm",247.07)
                           ("Bk",247.07)
                           ("Cf",251.08)
                           ("Es",254.09)
                           ("Fm",257.10)
                           ("Md",257.10)
                           ("No",255.09)
                           ("Lr",256.10)*/;
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
	  bondtable.insert(make_pair<string,double>(TypeA+"."+TypeB,coff));
	  bondtable.insert(make_pair<string,double>(TypeB+"."+TypeA,coff));
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
    using namespace boost::assign;
    bondtable = map_list_of
                           ("H.O",1.2)
			   ("O.H",1.2)
			   ("H.N",1.3)
			   ("N.H",1.3)
			   ("1.8",1.2)
			   ("8.1",1.2)
                           ("Si.O",1.9)
			   ("H.H",1.2)
			   ("O.Si",1.9)
			   ("X.H",1.1)
			   ("X.O",1.2)
		           ("X.X",1.1)
			   ("XX.XX",0.001)
			   ("XX.H",0.001)
			   ("XX.C",0.001)
		           ("XX.N",0.001)
			   ("XX.S",0.001)
                           ("XX.O",0.001)
                          ("C.Cl",1.98)
		           ("Cl.C",1.98)
                           ("H.Cl",1.4)
		           ("Cl.H",1.4)
                           ("C.C",1.7)
			   ("C.H",1.6)
                           ("H.C",1.6)

	//		   ("H.",1.9)
  //   	                  ("X.N",0.001)
	//		   ("X.C",0.001)
		           ("H.X",1.1);
	//		   ("O.X",0.001)
	//	           ("C.X",0.001)
	//		   ("N.X",0.001);
		   
			   // Add relevant bond lengths here
  }

};
