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

 		TCLAP::ValueArg<Real> tcrdist_sigma_arg("s","tcrdist_sigma",
			"tcrdist std-dev for gaussian match score", true,
			0, "float", cmd);

		// path to database files
 		TCLAP::ValueArg<std::string> db_filename_arg("d","db_filename",
			"Database file with info for tcrdist calculation", true,
			"", "string",cmd);

 		TCLAP::ValueArg<string> tcrs_file1_arg("i","tcrs_file1","TSV (tab separated values) "
			"file containing TCRs for neighbor calculation. Should contain the 4 columns "
			"'va_gene' 'cdr3a' 'vb_gene' 'cdr3b' (or alt fieldnames: 'va' and 'vb')", true,
			"unk", "string", cmd);

 		TCLAP::ValueArg<string> tcrs_file2_arg("j","tcrs_file2","TSV (tab separated values) "
			"file containing TCRs for neighbor calculation. Should contain the 4 columns "
			"'va_gene' 'cdr3a' 'vb_gene' 'cdr3b' (or alt fieldnames: 'va' and 'vb')", true,
			"unk", "string", cmd);

		cmd.parse( argc, argv );

		string const db_filename( db_filename_arg.getValue() );
		Real const threshold(tcrdist_threshold_arg.getValue()+0.5 );
		Real const sigma(tcrdist_sigma_arg.getValue());
		string const tcrs_file1( tcrs_file1_arg.getValue() );
		string const tcrs_file2( tcrs_file2_arg.getValue() );

		TCRdistCalculator const atcrdist('A', db_filename), btcrdist('B', db_filename);

		vector< PairedTCR > tcrs1, tcrs2;
		read_paired_tcrs_from_tsv_file(tcrs_file1, atcrdist, btcrdist, tcrs1);
		read_paired_tcrs_from_tsv_file(tcrs_file2, atcrdist, btcrdist, tcrs2);

		bools matched1(tcrs1.size(), false), matched2(tcrs2.size(), false);
		Size total_matches(0);
		Real match_score(0), dist(0);
		Real const neginvar( -1.0/(sigma*sigma) );

		for ( Size ii=0; ii< tcrs1.size(); ++ii ) {
			if (ii && ii%1000==0) cerr << '.';
			if (ii && ii%50000==0) cerr << endl;
			DistanceTCR_g const &atcr1( tcrs1[ii].first ), &btcr1( tcrs1[ii].second);
			for ( Size jj=0; jj< tcrs2.size(); ++jj ) {
				dist = atcrdist(atcr1, tcrs2[jj].first) + btcrdist(btcr1, tcrs2[jj].second);
				match_score += exp( neginvar*dist*dist );
				if ( dist <= threshold ) {
					matched1[ii] = true;
					matched2[jj] = true;
					++total_matches;
				}
			}
		}
		cerr << endl;

		Size total_matched1(0), total_matched2(0);
		for ( bool b : matched1 ) total_matched1 += b;
		for ( bool b : matched2 ) total_matched2 += b;

		cout << "matches: " << total_matches <<
			" match_score: " << match_score <<
			" matched_file1: " << total_matched1 <<
			" matched_file2: " << total_matched2 <<
			" ntcrs1: " << tcrs1.size() <<
			" ntcrs2: " << tcrs2.size() << endl;
		cout << "DONE" << endl;

	} catch (TCLAP::ArgException &e)  // catch any exceptions
			{ std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }

}


