// Miscellaneous generic helper functions
//

#ifndef INCLUDED_misc_HH
#define INCLUDED_misc_HH

#include "types.hh"

#include <random>
#include <iomanip>

#include <tclap/CmdLine.h>


void
my_exit(
	string const & file,
	int const line,
	string const & message
)
{

	ostringstream oss;
	if ( ! message.empty() ) oss << "\n" << "ERROR: " << message << "\n";
	oss << "ERROR:: Exit from: " << file << " line: " << line << "\n";
	string failure_message = oss.str();
	cerr << failure_message << flush;
	assert( false ); // for gdb?
	exit(1);
}

// this macro is borrowed from Rosetta
#define runtime_assert(_Expression) if ( !(_Expression) ) my_exit(__FILE__, __LINE__, #_Expression)

/// @details  Split a string to a vector
vector< string >
split_to_vector( string const & s )
{
	vector< string > v;

	istringstream l( s );
	string tag;
	l >> tag;
	while ( !l.fail() ) {
		v.push_back( tag );
		l >> tag;
	}
	return v;
}

/// @details  Split a string to a vector using sep as a separator
vector< string >
split_to_vector(
	string s,
	string const & sep
)
{
	vector< string > v;

	Size pos( s.find( sep ) );
	while ( pos != string::npos ) {
		v.push_back( s.substr( 0, pos ) );
		s.erase( 0, pos + sep.size() );
		pos = s.find( sep );
	}
	assert( s.find( sep ) == string::npos );
	v.push_back( s );
	return v;
}

// more hacky stuff
inline
bool
is_whitespace( string const & s )
{
	if ( s.empty() ) {
		return true;
	} else {
		return ( s.find_last_not_of( " \t\000" ) == string::npos );
	}
}

template< typename T >
inline
bool
is_type( string const & s )
{
	if ( is_whitespace( s ) ) {
		return false;
	} else {
		istringstream t_stream(s);
		T t;
		t_stream >> t;
		return ( ( t_stream ) && ( t_stream.eof() ) );
	}
}

inline
bool
is_int( string const & s )
{
	return is_type< int >( s );
}

inline
int
int_of( string const & s ) {
	istringstream l(s);
	int val;
	l >> val;
	runtime_assert( !l.fail() );
	return val;
}

inline
Real // actually a double
float_of( string const & s ) {
	istringstream l(s);
	Real val;
	l >> val;
	runtime_assert( !l.fail() );
	return val;
}


//
// eg, "V04"
//
string
get_v_family_from_v_gene( string const & g )
{
	runtime_assert( g[3] == 'V' );
	// runtime_assert( g.substr(0,2) == "TR" && g[3] == 'V' );
	// runtime_assert( g.substr(0,4) == "TRBV" );
	Size numlen(0);
	while ( is_int( g.substr(4,numlen+1) ) )++numlen;
	if (!numlen) { // for example, TRGVA (which is a pseudogene)
		return "V00";
	}
	runtime_assert( numlen );
	Size const vno( int_of( g.substr(4,numlen) ) );
	string const zeropad( vno < 10 ? "0" : "" );
	string const v_family( "V" + zeropad + to_string( vno ) );
	return v_family;
}


map<string,strings>
setup_v_family2v_genes( strings const & v_genes )
{
	map<string,strings> v_family2v_genes;
	for ( string g : v_genes ) {
		v_family2v_genes[ get_v_family_from_v_gene( g ) ].push_back( g );
	}
	return v_family2v_genes;
}

template < typename T1, typename T2 >
vector< T1 >
get_keys( map< T1, T2 > const & m )
{
	vector< T1 > ks;
	for ( typename map< T1,T2 >::const_iterator it= m.begin(); it != m.end(); ++it ) ks.push_back( it->first );
	return ks;
}


template < class T >
bool
has_element( vector< T > const & v, T const & t )
{
	return ( find( v.begin(), v.end(), t ) != v.end() );
}

/// never can remember the order...
template < class T >
bool
has_element( T const & t, vector< T > const & v )
{
	return ( find( v.begin(), v.end(), t ) != v.end() );
}


template < class T >
Size
vector_index( T const & t, vector< T > const & v )
{
	return find( v.begin(), v.end(), t ) - v.begin();
}



#endif
