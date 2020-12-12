// Miscellaneous generic helper functions
//

#ifndef INCLUDED_io_HH
#define INCLUDED_io_HH

#include "types.hh"
#include "tcrdist.hh"
#include "misc.hh"

// #include <random>
// #include <iomanip>

// silly helper function
Size
get_tsv_index(
	string const & headerline,
	strings const & fields
)
{
	strings const header(split_to_vector(headerline, "\t"));
	for ( string const & field : fields ) {
		if ( has_element(field, header) ) {
			return vector_index(field, header);
		}
	}
	cerr << "tcrs .tsv file is missing column. Possible fields ";
	for ( string const & field : fields ) cerr << ' '<< field;
	cerr << endl;
	exit(1);
	return 0;
}


//////////////////////////////////// READ THE TCRS
void
read_paired_tcrs_from_tsv_file(
	string const filename,
	TCRdistCalculator const & atcrdist,
	TCRdistCalculator const & btcrdist,
	vector< PairedTCR > & tcrs
)
{
	ifstream data(filename.c_str());
	if ( !data.good() ) {
		cerr<< "unable to open " << filename << endl;
		exit(1);
	}

	string line;
	getline(data, line);
	strings const header(split_to_vector(line, "\t") ); // for debugging only
	Size const va_index(get_tsv_index(line, split_to_vector("va va_gene")));
	Size const vb_index(get_tsv_index(line, split_to_vector("vb vb_gene")));
	Size const cdr3a_index(get_tsv_index(line, split_to_vector("cdr3a")));
	Size const cdr3b_index(get_tsv_index(line, split_to_vector("cdr3b")));

	while ( getline( data, line ) ) {
		strings const l(split_to_vector(line, "\t"));
		if ( l.size() != header.size() ) {
			cerr << "bad line length: " << line << endl;
			exit(1);
		}
		string const va(l[va_index]), vb(l[vb_index]), cdr3a(l[cdr3a_index]),
			cdr3b(l[cdr3b_index]);
		if ( !atcrdist.check_cdr3_ok(cdr3a) ) {
			cerr << "bad cdr3a: " << cdr3a << endl;
			exit(1);
		}
		if ( !btcrdist.check_cdr3_ok(cdr3b) ) {
			cerr << "bad cdr3b: " << cdr3b << endl;
			exit(1);
		}
		if ( !atcrdist.check_v_gene_ok(va)) {
			cerr << "bad va_gene: " << va << endl;
			exit(1);
		}
		if ( !btcrdist.check_v_gene_ok(vb)) {
			cerr << "bad vb_gene: " << vb << endl;
			exit(1);
		}

		tcrs.push_back(make_pair( atcrdist.create_distance_tcr_g(va, cdr3a),
				btcrdist.create_distance_tcr_g(vb, cdr3b)));
	}

	cout << "Read " << tcrs.size() << " paired tcrs from file " << filename << endl;
}

//////////////////////////////////// READ THE TCRS
void
read_single_chain_tcrs_from_tsv_file(
	string const filename,
	char const chain, // either 'A' or 'B'
	TCRdistCalculator const & tcrdist,
	vector< DistanceTCR_g > & tcrs
)
{
	runtime_assert( chain == 'A' || chain == 'B');

	ifstream data(filename.c_str());
	if ( !data.good() ) {
		cerr<< "unable to open " << filename << endl;
		exit(1);
	}

	string line;
	getline(data, line);
	strings const header(split_to_vector(line, "\t") ); // for debugging only
	string const vtags(chain == 'A' ? "v v_gene va va_gene" : "v v_gene vb vb_gene" );
	string const cdr3_tag(chain == 'A' ? "cdr3 cdr3a" : "cdr3 cdr3b" );
	Size const v_index(get_tsv_index(line, split_to_vector(vtags)));
	Size const cdr3_index(get_tsv_index(line, split_to_vector(cdr3_tag)));

	while ( getline( data, line ) ) {
		strings const l(split_to_vector(line, "\t"));
		if ( l.size() != header.size() ) {
			cerr << "bad line length: " << line << endl;
			exit(1);
		}
		string const v(l[v_index]), cdr3(l[cdr3_index]);
		if ( !tcrdist.check_cdr3_ok(cdr3) ) {
			cerr << "bad cdr3: " << cdr3 << endl;
			exit(1);
		}
		if ( !tcrdist.check_v_gene_ok(v)) {
			cerr << "bad v_gene: " << v << endl;
			exit(1);
		}

		tcrs.push_back(tcrdist.create_distance_tcr_g(v, cdr3));
	}

	cout << "Read " << tcrs.size() << " single chain tcrs from file " << filename << endl;
}


Sizes
read_groups_from_file( string const & filename )
{
	ifstream data(filename.c_str());
	if ( !data.good() ) {
		cerr<< "unable to open " << filename << endl;
		exit(1);
	}
	Sizes groups;
	Size g;
	while ( data.good() ) {
		data >> g;
		if ( !data.fail() ) {
			groups.push_back(g);
		} else {
			break;
		}
	}
	data.close();
	return groups;
}



#endif
