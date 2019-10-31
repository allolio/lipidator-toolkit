#define HAVE_BOOST
#include <string>
#include <fstream>
#include <vector>
#include <list>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include <algorithm>
#include <complex>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/random.hpp>
#include <boost/nondet_random.hpp>
//#define DEBUG
using namespace std;
#define SSTR( x ) static_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()
#define CACHEFRAMES 5
#define HAVE_CONFIG
#define STR_MAXCHARS 2000
#define INVALID_INDEX -1
//class XYZFile;
int RandomInt(int min, int max)
{
 static bool seeded=false;
 using namespace boost;
 static mt19937 gen;
  if(!seeded) {gen.seed(time(NULL));seeded=true;}  

 boost::uniform_int<> range(min, max);
 boost::variate_generator<boost::mt19937&, boost::uniform_int<> > next(gen, range);
 // this is optional, but wise to do
 return  next();  
}
//class XYZFile;
int Tokenize(const string& str,
                      vector<string>* tokens = NULL,
                      const string& delimiters = " ")
{
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);
    int ntokens=1;
    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        if (tokens!=NULL) tokens->push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
	ntokens++;
    }
    return ntokens;
}


#include "minimath.hpp"
/*#ifndef HAVE_FT
//#include "filetool.hpp"
#endif*/
#ifndef HAVE_MD

#include "topology.hpp"
#include "atominfo.hpp"
#include "md.hpp"
//#include "hbonding.hpp"
#endif
