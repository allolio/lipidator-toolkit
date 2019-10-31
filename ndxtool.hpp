#ifndef HAVE_CONFIG
#include "config.hpp"
#endif

typedef vector<int> atomlist;

atomlist NDX2Atomlist(string filename, string key)
{
  string buffer;
  string teststring="[ "+key+" ]" ;
  atomlist result;
  result.clear();
  ifstream ndx;
  ndx.open(filename.c_str());
  while (getline(ndx,buffer))
     if(buffer==teststring) break;
  if(buffer!=teststring) return result;
   while (getline(ndx,buffer))
   { vector<string> tokens;
     Tokenize(buffer,&tokens);
   if(tokens.size()==0) break;
   if(!isdigit( *((tokens[0]).c_str()))) break;
      for(int i=0;i!=tokens.size();i++)
	{
	  result.push_back(atoi(tokens[i].c_str())-1);
	}
   }
   ndx.close();
   return  result;
}
