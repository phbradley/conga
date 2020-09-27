// some "useful" typedefs and simple classes
//

#ifndef INCLUDED_types_HH
#define INCLUDED_types_HH


// #include <boost/foreach.hpp>
// #include <boost/algorithm/string.hpp>

#include <iostream>
#include <fstream>
#include <sstream>
// #include <iterator>
#include <algorithm>
#include <string>
#include <map>
#include <vector>
#include <cassert>

#include <math.h>

using namespace std;

// #define foreach_ BOOST_FOREACH

typedef std::size_t  Size;
typedef double  Real;
typedef std::pair< Size, Size > SizePair;

typedef vector< SizePair > SizePairs;
typedef vector< string > strings;
typedef vector< Size > Sizes;
typedef vector< bool > bools;
typedef vector< Real > Reals;

// struct Feature {
// 	string name;
// 	Sizes poslist;
// 	Sizes neglist;
// };

// typedef vector< Feature > Features;



#endif
