#include "types.hh"
#include "tcrdist.hh"
#include "io.hh"
#include <random>



int main(int argc, char** argv)
{
	try { // to catch tclap exceptions

		TCLAP::CmdLine cmd( "tcrdist_distributions",' ', "0.1" );

 		TCLAP::ValueArg<Size> max_dist_arg("m","max_dist",
			"maximum tcrdist distance over which to compute background paired distance distributions", true,
			0, "integer", cmd);

		// path to database files
 		TCLAP::ValueArg<std::string> db_filename_arg("d","db_filename",
			"Database file with info for tcrdist calculation", true,
			"", "string",cmd);

 		TCLAP::ValueArg<std::string> outfile_arg("o","outfile",
			"Filename to use for writing the counts distribution",true,
			"", "string",cmd);

 		TCLAP::ValueArg<string> tcrs_file_arg("f","tcrs_file","TSV (tab separated values) "
			"file containing TCRs for neighbor calculation. Should contain the 4 columns "
			"'va_gene' 'cdr3a' 'vb_gene' 'cdr3b' (or alt fieldnames: 'va' and 'vb')", true,
			"unk", "string", cmd);

 		TCLAP::ValueArg<string> achains_file_arg("a", "achains_file",
			"TSV file with the background TCR alpha chains", true, "", "string", cmd);

 		TCLAP::ValueArg<string> bchains_file_arg("b", "bchains_file",
			"TSV file with the background TCR beta chains", true, "", "string", cmd);


		cmd.parse( argc, argv );

		string const db_filename( db_filename_arg.getValue() );
		Size const max_dist( max_dist_arg.getValue() );
		string const tcrs_file( tcrs_file_arg.getValue() );
		string const achains_file( achains_file_arg.getValue() );
		string const bchains_file( bchains_file_arg.getValue() );
		string const outfile( outfile_arg.getValue());

		TCRdistCalculator const atcrdist('A', db_filename), btcrdist('B', db_filename);

		vector< PairedTCR > tcrs;
		read_paired_tcrs_from_tsv_file(tcrs_file, atcrdist, btcrdist, tcrs);

		vector< DistanceTCR_g > achains, bchains;

		read_single_chain_tcrs_from_tsv_file(achains_file, 'A', atcrdist, achains);
		read_single_chain_tcrs_from_tsv_file(bchains_file, 'B', btcrdist, bchains);

		Size const num_tcrs(tcrs.size());

		Sizes acounts(max_dist+1), bcounts(max_dist+1);
		ofstream out(outfile);

		cout << "making " << outfile << endl;

		for ( Size ii=0; ii< num_tcrs; ++ii ) {
			if ( ii && ii%100==0 ) cerr << '.';
			if ( ii && ii%5000==0 ) cerr << ' ' << ii << endl;

			//DistanceTCR_g const &atcr( tcrs[ii].first ), &btcr( tcrs[ii].second);

			for ( Size r=0; r<2; ++r){
				Sizes & counts( r==0 ? acounts : bcounts);
				fill( counts.begin(), counts.end(), 0);
				DistanceTCR_g const &fg_tcr( r==0 ? tcrs[ii].first : tcrs[ii].second);
				vector<DistanceTCR_g> const & bg_tcrs( r==0 ? achains : bchains );
				TCRdistCalculator const & tcrdist( r==0 ? atcrdist : btcrdist );
				Size dist(0);
				for ( DistanceTCR_g const & bg_tcr : bg_tcrs ) {
					dist = Size( 0.5 + tcrdist(fg_tcr, bg_tcr));
					if ( dist <= max_dist ) ++counts[dist];
				}
			}

			// compute probability distribution for paired distances using convolution
			for ( Size d=0; d<= max_dist; ++d ) {
				Size count(0);
				for ( Size adist=0; adist<= d; ++adist ) {
					count += acounts[adist] * bcounts[d-adist];
				}
				if (d) out << ' ';
				out << count;
			}
			out << '\n';
		}

		cerr << endl;
		out.close();

	} catch (TCLAP::ArgException &e)  // catch any exceptions
		{ std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }

}
