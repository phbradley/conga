#include "types.hh"
#include "tcrdist.hh"
#include "io.hh"


int main(int argc, char** argv)
{
	try { // to catch tclap exceptions

		TCLAP::CmdLine cmd( "compute a matrix of single-chain TCRdist "
			"distances",' ', "0.1" );

		// path to database files
 		TCLAP::ValueArg<std::string> db_filename_arg("d", "db_filename",
			"Database file with info for tcrdist calculation",
			true, "", "string", cmd);

 		TCLAP::ValueArg<string> chain_arg("c", "chain",
			"TCR chain (A or B)",
			true, "", "string", cmd);

 		TCLAP::ValueArg<string> outfile_arg("o", "outfile",
			"The name of the file to which the pairwise distance matrix "
			"will be written. In this file, the rows will correspond to the "
			"TCRs in tcrs_file1 and the columns will correspond to the TCRs "
			"in tcrs_file2",
			true, "", "string", cmd);

 		TCLAP::ValueArg<string> tcrs_file2_arg("j", "tcrs_file2",
			"TSV (tab separated values) file #2 containing TCRs "
			"for pairwise distance neighbor calculation. Should contain "
			"columns giving the V gene (valid headers v, va, vb, va_gene, vb_gene) "
			"and the CDR3 sequence (valid headers cdr3, cdr3a, cdr3b)",
			true, "unk", "string", cmd);

 		TCLAP::ValueArg<string> tcrs_file1_arg("i", "tcrs_file1",
			"TSV (tab separated values) file #1 containing TCRs "
			"for pairwise distance neighbor calculation. Should contain "
			"columns giving the V gene (valid headers v, va, vb, va_gene, vb_gene) "
			"and the CDR3 sequence (valid headers cdr3, cdr3a, cdr3b)",
			true, "unk", "string", cmd);


		cmd.parse( argc, argv );

		string const db_filename( db_filename_arg.getValue() );
		string const tcrs_file1( tcrs_file1_arg.getValue() );
		string const tcrs_file2( tcrs_file2_arg.getValue() );
		string const outfile( outfile_arg.getValue() );
		char const chain(chain_arg.getValue()[0]);
		runtime_assert(chain == 'A' || chain == 'B');

		TCRdistCalculator const tcrdist(chain, db_filename);

		vector< DistanceTCR_g > tcrs1, tcrs2;

		read_single_chain_tcrs_from_tsv_file(tcrs_file1, chain, tcrdist, tcrs1);
		read_single_chain_tcrs_from_tsv_file(tcrs_file2, chain, tcrdist, tcrs2);


		ofstream out(outfile.c_str());

		for ( Size ii=0; ii< tcrs1.size(); ++ii ) {
			if (ii && ii%1000==0) cerr << '.';
			if (ii && ii%50000==0) cerr << endl;
			DistanceTCR_g const & tcr1(tcrs1[ii]);
			for ( Size jj=0; jj< tcrs2.size(); ++jj ) {
				if (jj) out << ' ';
				out << (int)(tcrdist(tcr1, tcrs2[jj])+0.1);
			}
			out << '\n';
		}
		out.close();
		cerr << endl;

	} catch (TCLAP::ArgException &e)  // catch any exceptions
		{
			std::cerr << "error: " << e.error() << " for arg " << e.argId() <<
				std::endl;
		}

}


