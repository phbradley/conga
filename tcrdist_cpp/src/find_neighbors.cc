#include "types.hh"
#include "tcrdist.hh"
#include "io.hh"
#include <random>


int main(int argc, char** argv)
{
	try { // to catch tclap exceptions

		TCLAP::CmdLine cmd( "find_neighbors. Use either --num_nbrs or --threshold",' ', "0.1" );

 		TCLAP::ValueArg<Size> num_nbrs_arg("n","num_nbrs",
			"Number of nearest neighbors to find (not including self). Alternative to "
			"using --threshold.", false,
			0, "integer", cmd);

 		TCLAP::ValueArg<int> threshold_arg("t","threshold",
			"TCRdist threshold for neighborness (alternative to using --num_nbrs) -- should be an INTEGER", false,
			-1, "integer", cmd);

		// path to database files
 		TCLAP::ValueArg<std::string> db_filename_arg("d","db_filename",
			"Database file with info for tcrdist calculation", true,
			"", "string",cmd);

 		TCLAP::ValueArg<std::string> outfile_prefix_arg("o","outfile_prefix",
			"Prefix for the knn_indices and knn_distances output files",true,
			"", "string",cmd);

		TCLAP::SwitchArg only_tcrdists_arg("m","only_tcrdists", "Just write the matrix of tcrdists. "
			"Don't find neighbors. Matrix will be called <outfile_prefix>_tcrdists.txt", cmd, false);

 		TCLAP::ValueArg<string> tcrs_file_arg("f","tcrs_file","TSV (tab separated values) "
			"file containing TCRs for neighbor calculation. Should contain the 4 columns "
			"'va_gene' 'cdr3a' 'vb_gene' 'cdr3b' (or alt fieldnames: 'va' and 'vb')", true,
			"unk", "string", cmd);

 		TCLAP::ValueArg<string> agroups_file_arg("a","agroups_file","np.savetxt output "
			"(ie, one integer per line) ith the agroups information so we can exclude same-group neighbors", false,
			"", "string", cmd);

 		TCLAP::ValueArg<string> bgroups_file_arg("b","bgroups_file","np.savetxt output "
			"(ie, one integer per line) ith the bgroups information so we can exclude same-group neighbors", false,
			"", "string", cmd);

		cmd.parse( argc, argv );

		string const db_filename( db_filename_arg.getValue() );
		Size const num_nbrs( num_nbrs_arg.getValue() );
		int const threshold_int( threshold_arg.getValue() );
		bool const only_tcrdists( only_tcrdists_arg.getValue() );
		string const tcrs_file( tcrs_file_arg.getValue() );
		string const agroups_file( agroups_file_arg.getValue() );
		string const bgroups_file( bgroups_file_arg.getValue() );
		string const outfile_prefix( outfile_prefix_arg.getValue());

		runtime_assert( only_tcrdists || ( num_nbrs>0 && threshold_int==-1) || (num_nbrs==0 && threshold_int >=0 ) );

		TCRdistCalculator const atcrdist('A', db_filename), btcrdist('B', db_filename);

		vector< PairedTCR > tcrs;
		read_paired_tcrs_from_tsv_file(tcrs_file, atcrdist, btcrdist, tcrs);

		Size const num_tcrs(tcrs.size());

		Sizes agroups( agroups_file.size() ? read_groups_from_file(agroups_file) : Sizes() );
		Sizes bgroups( bgroups_file.size() ? read_groups_from_file(bgroups_file) : Sizes() );
		if ( agroups.empty() ) {
			for ( Size i=0; i<num_tcrs; ++i ) agroups.push_back(i);
		}
		if ( bgroups.empty() ) {
			for ( Size i=0; i<num_tcrs; ++i ) bgroups.push_back(i);
		}

		runtime_assert( agroups.size() == num_tcrs );
		runtime_assert( bgroups.size() == num_tcrs );
		// two different modes of operations

		Size const BIG_DIST(10000);

		if ( only_tcrdists ) {
			ofstream out(outfile_prefix+"_tcrdists.txt");
			cout << "making " << outfile_prefix+"_tcrdists.txt" << endl;
			for ( Size ii=0; ii< num_tcrs; ++ii ) {
				if ( ii && ii%100==0 ) cerr << '.';
				if ( ii && ii%5000==0 ) cerr << ' ' << ii << endl;

				DistanceTCR_g const &atcr( tcrs[ii].first ), &btcr( tcrs[ii].second);
				bool first(true);
				for ( PairedTCR const & other_tcr : tcrs ) {
					// NOTE we round down to an integer here!
					if ( first ) first=false;
					else out << ' ';
					out << Size( 0.5 + atcrdist(atcr, other_tcr.first) + btcrdist(btcr, other_tcr.second) );
				}
				out << '\n';
			}
			cerr << endl;
			out.close();

		} else if ( num_nbrs > 0 ) {
			// open the outfiles
			ofstream out_indices(outfile_prefix+"_knn_indices.txt");
			ofstream out_distances(outfile_prefix+"_knn_distances.txt");

			cout << "making " << outfile_prefix+"_knn_indices.txt" << " and " <<
				outfile_prefix+"_knn_distances.txt" << endl;

			Sizes dists(num_tcrs), sortdists(num_tcrs); // must be a better way to do this...
			Sizes knn_indices, knn_distances;
			knn_indices.reserve(num_nbrs);
			knn_distances.reserve(num_nbrs);

			minstd_rand0 rng(1); // seed
			Sizes shuffled_indices;
			for ( Size i=0; i<num_tcrs; ++i ) shuffled_indices.push_back(i);

			for ( Size ii=0; ii< num_tcrs; ++ii ) {
				if ( ii && ii%100==0 ) cerr << '.';
				if ( ii && ii%5000==0 ) cerr << ' ' << ii << endl;

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
				Size const a(agroups[ii]), b(bgroups[ii]);
				for ( Size jj=0; jj< num_tcrs; ++jj ) {
					if ( agroups[jj] == a || bgroups[jj] == b ) dists[jj] = BIG_DIST;
				}
				runtime_assert( dists[ii] == BIG_DIST );
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
				// if ( ii%500==0 ) {
				// 	cout << "threshold: " << threshold << " num_at_threshold: " <<
				// 		num_at_threshold << " ii: " << ii << endl;
				// }
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
			cerr << endl;
			// close the output files
			out_indices.close();
			out_distances.close();

		} else { // using threshold definition of nbr-ness
			// open the outfiles
			ofstream out_indices(outfile_prefix+"_nbr"+to_string(threshold_int)+"_indices.txt");
			ofstream out_distances(outfile_prefix+"_nbr"+to_string(threshold_int)+"_distances.txt");

			cout << "making " << outfile_prefix+"_nbr"+to_string(threshold_int)+"_indices.txt" << " and " <<
				outfile_prefix+"_nbr"+to_string(threshold_int)+"_distances.txt" << endl;

			runtime_assert( threshold_int >= 0 );
			Size const threshold(threshold_int);

			Sizes knn_indices, knn_distances;
			knn_indices.reserve(num_tcrs);
			knn_distances.reserve(num_tcrs);

			for ( Size ii=0; ii< num_tcrs; ++ii ) {
				if ( ii && ii%100==0 ) cerr << '.';
				if ( ii && ii%5000==0 ) cerr << ' ' << ii << endl;
				knn_indices.clear();
				knn_distances.clear();

				// for ties, shuffle so we don't get biases based on file order
				DistanceTCR_g const &atcr( tcrs[ii].first ), &btcr( tcrs[ii].second);
				Size const a(agroups[ii]), b(bgroups[ii]);
				for ( Size jj=0; jj< num_tcrs; ++jj ) {
					Size const dist( 0.5 + atcrdist(atcr, tcrs[jj].first) + btcrdist(btcr, tcrs[jj].second) );
					if ( dist <= threshold && agroups[jj] != a && bgroups[jj] != b ) {
						knn_indices.push_back(jj);
						knn_distances.push_back(dist);
					}
				}

				// save to file: note that these lines may be empty!!!
				for ( Size j=0; j<knn_indices.size(); ++j ) {
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
			cerr << endl;
			// close the output files
			out_indices.close();
			out_distances.close();
		}

	} catch (TCLAP::ArgException &e)  // catch any exceptions
		{ std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }

}
