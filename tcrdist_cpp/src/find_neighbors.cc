#include "types.hh"
#include "tcrdist.hh"
#include <random>

// a paired tcr with gene-level (actually allele level) resolution
// DistanceTCR_g is defined in tcrdist.hh
typedef std::pair< DistanceTCR_g, DistanceTCR_g > PairedTCR;

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





int main(int argc, char** argv)
{
	try { // to catch tclap exceptions

		TCLAP::CmdLine cmd( "find_neighbors",' ', "0.1" );

 		TCLAP::ValueArg<Size> num_nbrs_arg("n","num_nbrs",
			"Number of nearest neighbors to find (not including self)", true,
			0, "integer", cmd);

		// path to database files
 		TCLAP::ValueArg<std::string> db_filename_arg("d","db_filename",
			"Database file with info for tcrdist calculation", true,
			"", "string",cmd);

 		TCLAP::ValueArg<std::string> outfile_prefix_arg("o","outfile_prefix",
			"Prefix for the knn_indices and knn_distances output files",true,
			"", "string",cmd);


 		TCLAP::ValueArg<string> tcrs_file_arg("f","tcrs_file","TSV (tab separated values) "
			"file containing TCRs for neighbor calculation. Should contain the 4 columns "
			"'va_gene' 'cdr3a' 'vb_gene' 'cdr3b' (or alt fieldnames: 'va' and 'vb')", true,
			"unk", "string", cmd);

		cmd.parse( argc, argv );

		string const db_filename( db_filename_arg.getValue() );
		Size const num_nbrs( num_nbrs_arg.getValue() );
		string const tcrs_file( tcrs_file_arg.getValue() );
		string const outfile_prefix( outfile_prefix_arg.getValue());

		TCRdistCalculator const atcrdist('A', db_filename), btcrdist('B', db_filename);

		vector< PairedTCR > tcrs;
		read_paired_tcrs_from_tsv_file(tcrs_file, atcrdist, btcrdist, tcrs);

		Size const num_tcrs(tcrs.size());

		Sizes dists(num_tcrs), sortdists(num_tcrs); // must be a better way to do this...
		Sizes knn_indices, knn_distances;
		knn_indices.reserve(num_nbrs);
		knn_distances.reserve(num_nbrs);

		minstd_rand0 rng(1); // seed
		Sizes shuffled_indices;
		for ( Size i=0; i<num_tcrs; ++i ) shuffled_indices.push_back(i);

		ofstream out_indices(outfile_prefix+"_knn_indices.txt");
		ofstream out_distances(outfile_prefix+"_knn_distances.txt");

		for ( Size ii=0; ii< num_tcrs; ++ii ) {
			// for ties, shuffle so we don't get biases based on file order
			shuffle(shuffled_indices.begin(), shuffled_indices.end(), rng);
			DistanceTCR_g const &atcr( tcrs[ii].first ), &btcr( tcrs[ii].second);
			{
				Size i(0);
				for ( PairedTCR const & other_tcr : tcrs ) {
					// NOTE we round down to an integer here!
					dists[i] = Size( 0.5 + atcrdist(atcr, other_tcr.first) + btcrdist(btcr, other_tcr.second) );
					++i;
				}
			}
			dists[ii] = 10000; // something big
			copy(dists.begin(), dists.end(), sortdists.begin());
			nth_element(sortdists.begin(), sortdists.begin()+num_nbrs-1, sortdists.end());
			Size const threshold(sortdists[num_nbrs-1]);
			Size num_at_threshold(0);
			for ( Size i=0; i<num_nbrs; ++i ) {
				// runtime_assert( sortdists[i] <= threshold ); // for debugging
				if ( sortdists[i] == threshold ) ++num_at_threshold;
			}
			// for ( Size i=num_nbrs; i< num_tcrs; ++i ) { // just for debugging
			// 	runtime_assert( sortdists[i] >= threshold );
			// }
			if ( ii%500==0 ) {
				cout << "threshold: " << threshold << " num_at_threshold: " <<
					num_at_threshold << " ii: " << ii << endl;
			}
			knn_distances.clear();
			knn_indices.clear();
			for ( Size i : shuffled_indices ) {
				if ( dists[i] < threshold ) {
					knn_indices.push_back(i);
					knn_distances.push_back(dists[i]);
				} else if ( dists[i] == threshold && num_at_threshold>0 ) {
					knn_indices.push_back(i);
					knn_distances.push_back(dists[i]);
					--num_at_threshold;
				}
			}
			runtime_assert(knn_indices.size() == num_nbrs);
			runtime_assert(knn_distances.size() == num_nbrs);
			// save to files:
			for ( Size j=0; j<num_nbrs; ++j ) {
				if (j) {
					out_indices << ' ';
					out_distances << ' ';
				}
				out_indices << knn_indices[j];
				out_distances << knn_distances[j];
			}
			out_indices << '\n';
			out_distances << '\n';
		}

		out_indices.close();
		out_distances.close();

	} catch (TCLAP::ArgException &e)  // catch any exceptions
		{ std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }

}
