#include "types.hh"
#include "tcrdist.hh"
#include "io.hh"
#include <random>



int main(int argc, char** argv)
{
	try { // to catch tclap exceptions

		TCLAP::CmdLine cmd( "tcrdist_distributions",' ', "0.1" );

 		TCLAP::ValueArg<Size> tcrdist_threshold_arg("t","tcrdist_threshold",
			"tcrdist threshold for matches", true,
			0, "integer", cmd);

		// path to database files
 		TCLAP::ValueArg<std::string> db_filename_arg("d","db_filename",
			"Database file with info for tcrdist calculation", true,
			"", "string",cmd);

 		TCLAP::ValueArg<string> tcrs_file_arg("f","tcrs_file","TSV (tab separated values) "
			"file containing TCRs for neighbor calculation. Should contain the 4 columns "
			"'va_gene' 'cdr3a' 'vb_gene' 'cdr3b' (or alt fieldnames: 'va' and 'vb')", true,
			"unk", "string", cmd);

 		TCLAP::ValueArg<std::string> outfile_prefix_arg("o","outfile_prefix",
			"Prefix for the knn_indices and knn_distances output files",true,
			"", "string",cmd);

 		TCLAP::ValueArg<string> chain_arg("c", "chain",
			"TCR chain (A or B)", true, "", "string", cmd);

		cmd.parse( argc, argv );

		string const db_filename( db_filename_arg.getValue() );
		Size const threshold_int(tcrdist_threshold_arg.getValue());
		Real const threshold(threshold_int+0.5);
		string const tcrs_file( tcrs_file_arg.getValue() );
		char const chain(chain_arg.getValue()[0]);
		runtime_assert(chain == 'A' || chain == 'B');
		string const outfile_prefix( outfile_prefix_arg.getValue());

		TCRdistCalculator const tcrdist(chain, db_filename);

		vector< DistanceTCR_g > tcrs;

		read_single_chain_tcrs_from_tsv_file(tcrs_file, chain, tcrdist, tcrs);

		string const
			indfile (outfile_prefix+"_nbr"+to_string(threshold_int)+"_indices.txt"),
			distfile(outfile_prefix+"_nbr"+to_string(threshold_int)+"_distances.txt");

		ofstream out_indices(indfile);
		ofstream out_distances(distfile);

		cout << "making " << indfile << " and " << distfile <<endl;

		Size const num_tcrs(tcrs.size());
		Sizes knn_indices, knn_distances;
		knn_indices.reserve(num_tcrs);
		knn_distances.reserve(num_tcrs);

		Real dist(0);

		for ( Size ii=0; ii< num_tcrs; ++ii ) {
			if (ii && ii%1000==0) cerr << '.';
			if (ii && ii%50000==0) cerr << endl;
			knn_indices.clear();
			knn_distances.clear();
			DistanceTCR_g const & tcr1(tcrs[ii]);
			for ( Size jj=0; jj< num_tcrs; ++jj ) {
				dist = tcrdist(tcr1, tcrs[jj]);
				if ( dist <= threshold ) {
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

		cout << "DONE" << endl;

	} catch (TCLAP::ArgException &e)  // catch any exceptions
		{ std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }

}


