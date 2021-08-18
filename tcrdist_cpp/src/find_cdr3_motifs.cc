#include "types.hh"
#include "tcrdist.hh"
#include "io.hh"
#include <random>
#include <boost/math/distributions/hypergeometric.hpp>

Size const MIN_FG_COUNT(5);
Real const MIN_CORE_ENRICH(1.05);
Real const MIN_PAIR_ENRICH(1.5);
Real const MIN_TRIPLE_ENRICH(2.0);
Real const MAX_CORE_PVAL_ADJ(100.);
Real const MAX_PAIR_PVAL_ADJ(1.);
Real const MAX_TRIPLE_PVAL_ADJ(1.);
Real const MAX_FINAL_PVAL_ADJ(1.);
Real const MAX_JACCARD(0.7);
Size const NUM_CORE_AAS(26); // 'A'-'Z'
Size const NUM_CORE_GAPS(3); // 0,1,2

Size const NUM_CORES(
	NUM_CORE_AAS*NUM_CORE_AAS*NUM_CORE_AAS*NUM_CORE_GAPS*NUM_CORE_GAPS);

Real const PSEUDOCOUNT(0.25);

typedef vector< char > chars;

inline
Size
get_core_index(
	char const aa1, // 'A' - 'Z'
	Size const gap1, // 0-2
	char const aa2, // 'A' - 'Z'
	Size const gap2, // 0-2
	char const aa3  // 'A' - 'Z'
)
{
	// 26 * 3 * 26 * 3 * 26
	return  (
		NUM_CORE_AAS*NUM_CORE_GAPS*NUM_CORE_AAS*NUM_CORE_GAPS*(aa1-'A')+
		NUM_CORE_AAS*NUM_CORE_GAPS*NUM_CORE_AAS*gap1+
		NUM_CORE_AAS*NUM_CORE_GAPS*(aa2-'A')+
		NUM_CORE_AAS*gap2+
		(aa3-'A'));
}


struct Motif {
	chars aas;
	Sizes gaps;
	//Size span;
	//Size core_index;

	Motif():
		aas(),
		gaps()
		//core_index(0)
	{}


	Motif(
		char const aa1,
		Size const gap1,
		char const aa2,
		Size const gap2,
		char const aa3
	){
		aas.resize(3);
		aas[0] = aa1;
		aas[1] = aa2;
		aas[2] = aa3;
		gaps.resize(2);
		gaps[0] = gap1;
		gaps[1] = gap2;
		//core_index = get_core_index(aa1,gap1,aa2,gap2,aa3);
		//span = gaps[0] + gaps[1] + 1; // gap between first and last position
	}

	// total gap between aa_i and aa_j
	Size
	gap(Size const i, Size const j) const
	{
		Size total(0);
		for ( Size k=i; k<j; ++k)
			total += gaps[k];
		return total+(j-i)-1;
	}

	string
	to_string() const
	{
		runtime_assert(aas.size() == gaps.size()+1);
		ostringstream out;
		for (Size k=0; k<aas.size(); ++k ){
			out << aas[k];
			if (k<gaps.size()) {
				for (Size g=0; g<gaps[k]; ++g) out << '.';
			}
		}
		return out.str();
	}
};


inline
Real
compute_overlap_pvalue(
	Size const overlap,
	Size const count1,
	Size const count2,
	Size const total
)
{
	using namespace boost::math;
	// whats the smallest possible overlap? max( 0, count1+count2-total )
	if ( overlap==0 || count1==total || count2==total ||
		count1+count2 >= total+overlap ) return 1.0;
	hypergeometric_distribution<> hgd( count1, count2, total );
	// need overlap-1 to be a valid overlap value:
	return cdf( complement( hgd, overlap-1 ) );
}


inline
void
unpack_core_index(
	Size ind,
	char & aa1,
	Size & gap1,
	char & aa2,
	Size & gap2,
	char & aa3
)
{
	Size aa3_i = ind%NUM_CORE_AAS;
	ind = (ind-aa3_i)/NUM_CORE_AAS;
	gap2 = ind%NUM_CORE_GAPS;
	ind = (ind-gap2)/NUM_CORE_GAPS;
	Size aa2_i = ind%NUM_CORE_AAS;
	ind = (ind-aa2_i)/NUM_CORE_AAS;
	gap1 = ind%NUM_CORE_GAPS;
	ind = (ind-gap1)/NUM_CORE_GAPS;
	Size aa1_i = ind%NUM_CORE_AAS;

	aa1 = char('A'+aa1_i);
	aa2 = char('A'+aa2_i);
	aa3 = char('A'+aa3_i);
}

inline
string
core_index_to_string(
	Size ind
)
{
	char aa1, aa2, aa3;
	Size gap1, gap2;

	unpack_core_index(ind, aa1, gap1, aa2, gap2, aa3);
	ostringstream out;
	out << aa1;
	for (Size g=0; g<gap1; ++g) out << '.';
	out << aa2;
	for (Size g=0; g<gap2; ++g) out << '.';
	out << aa3;
	return out.str();
}

inline
Motif
core_index_to_motif(
	Size ind
)
{
	char aa1, aa2, aa3;
	Size gap1, gap2;

	unpack_core_index(ind, aa1, gap1, aa2, gap2, aa3);
	return Motif(aa1, gap1, aa2, gap2, aa3);
}

inline
bool
motifs_overlap(
	Motif const & a,
	Motif const & b,
	Size & a1,
	Size & a2,
	Size & b1,
	Size & b2
)
{
	// NOTE -- there might be more than one match? if dup aas?
	//
	// a.aas[a1] == b.aas[b1]
	// a.aas[a2] == b.aas[b2]
	//

	for ( a1=0; a1<1; ++a1 ){
		for ( b1=0; b1<1; ++b1 ){
			if (a.aas[a1] == b.aas[b1] ) {
				for ( a2=a1+1; a2<2; ++a2 ){
					for ( b2=b1+1; b2<2; ++b2 ){
						if (a.aas[a2] == b.aas[b2] &&
							a.gap(a1,a2) == b.gap(b1,b2)){
							// actually they could be in conflict
							if (a1==b1 && a2==b2 && a.gaps[0]==b.gaps[0] &&
								a.gaps[1]==b.gaps[1]) return false; // unless they are identical
							return true; // and a1,a2,b1,b2 are set
						}
					}
				}
			}
		}
	}
	return false;
}

inline
bool
motifs_overlap(
	Motif const & a,
	Motif const & b
)
{
	static Size a1, a2, b1, b2;
	return motifs_overlap(a,b,a1,a2,b1,b2);
}


void
get_core_counts_for_cdr3s(
	strings const & cdr3s,
	Sizes & counts
)
{
	counts.clear();
	counts.resize(NUM_CORES);

	cout << "get_core_counts_for_cdr3s: num cdr3s= " << cdr3s.size() << endl;

	for (string const & cdr3 : cdr3s ) {
		for ( Size i=0; i<cdr3.size()-2; ++i ) {
			for ( Size j=i+1; j<cdr3.size()-1 && j<=i+NUM_CORE_GAPS; ++j ) {
				for ( Size k=j+1; k<cdr3.size() && k<=j+NUM_CORE_GAPS; ++k ) {
					++counts[get_core_index(cdr3[i], j-i-1, cdr3[j], k-j-1, cdr3[k])];
				}
			}
		}
	}


}

//typedef map<Size, set<Size>> AllOccs;
typedef map<Size, Sizes> AllOccs;

void
get_core_occurrences_for_cdr3s(
	strings const & cdr3s,
	AllOccs & all_occs
)
{
	cout << "get_core_occurrences_for_cdr3s: num cdr3s= " << cdr3s.size() << endl;

	AllOccs::iterator it;

	Size ind(0);
	for (Size icdr3=0; icdr3< cdr3s.size(); ++icdr3) {
		string const & cdr3(cdr3s[icdr3]);
		for ( Size i=0; i<cdr3.size()-2; ++i ) {
			for ( Size j=i+1; j<cdr3.size()-1 && j<=i+NUM_CORE_GAPS; ++j ) {
				for ( Size k=j+1; k<cdr3.size() && k<=j+NUM_CORE_GAPS; ++k ) {
					ind = get_core_index(cdr3[i], j-i-1, cdr3[j], k-j-1, cdr3[k]);
					it = all_occs.find(ind);
					if (it != all_occs.end()) it->second.push_back(icdr3);
				}
			}
		}
	}


}

// assumes that a and b are sorted!!!
Real
compute_jaccard(
	Sizes const & a,
	Sizes const & b
)
{
	Sizes overlap(max(a.size(),b.size()));

	Sizes::iterator const ite(
		set_intersection(a.begin(), a.end(), b.begin(), b.end(), overlap.begin()));
	Size overlap_size(ite - overlap.begin());
	return Real(overlap_size)/(a.size() + b.size() - overlap_size);
}


int main(int argc, char** argv)
{
	try { // to catch tclap exceptions

		TCLAP::CmdLine cmd( "tcrdist_distributions",' ', "0.1" );

		// path to database files
 		TCLAP::ValueArg<std::string> db_filename_arg("d","db_filename",
			"Database file with info for tcrdist calculation", true,
			"", "string",cmd);

 		TCLAP::ValueArg<string> fg_tcrs_file_arg("f","fg_tcrs_file","TSV (tab separated values) "
			"file containing TCRs for neighbor calculation. Should contain the 4 columns "
			"'va_gene' 'cdr3a' 'vb_gene' 'cdr3b' (or alt fieldnames: 'va' and 'vb')", true,
			"unk", "string", cmd);

 		TCLAP::ValueArg<string> bg_tcrs_file_arg("b","bg_tcrs_file","TSV (tab separated values) "
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
		string const fg_tcrs_file( fg_tcrs_file_arg.getValue() );
		string const bg_tcrs_file( bg_tcrs_file_arg.getValue() );
		char const chain(chain_arg.getValue()[0]);
		runtime_assert(chain == 'A' || chain == 'B');
		string const outfile_prefix( outfile_prefix_arg.getValue());

		TCRdistCalculator const tcrdist(chain, db_filename);

		// tmp hacking-- we really only need the CDR3, right?
		vector< DistanceTCR_g > fg_tcrs, bg_tcrs;

		read_single_chain_tcrs_from_tsv_file(fg_tcrs_file, chain, tcrdist, fg_tcrs);
		read_single_chain_tcrs_from_tsv_file(bg_tcrs_file, chain, tcrdist, bg_tcrs);

		strings fg_cdr3s, bg_cdr3s;
		for (auto t : fg_tcrs) {
			fg_cdr3s.push_back("B"+t.cdr3.substr(3,t.cdr3.size()-5)+"Z");
		}
		//cout <<fg_cdr3s.front() << ' ' << fg_cdr3s.back() << endl;

		for (auto t : bg_tcrs) {
			bg_cdr3s.push_back("B"+t.cdr3.substr(3,t.cdr3.size()-5)+"Z");
		}

		Real const fg_bg_ratio(Real(fg_cdr3s.size())/bg_cdr3s.size());

		//
		Sizes fg_counts, bg_counts;

		get_core_counts_for_cdr3s(fg_cdr3s, fg_counts);
		get_core_counts_for_cdr3s(bg_cdr3s, bg_counts);


		Size total_core_pvals(0);
		for ( Size ind=0; ind<NUM_CORES; ++ind ) {
			if ( fg_counts[ind] >= MIN_FG_COUNT ) ++total_core_pvals;
		}

		vector<pair<Real, Size> > sortl;
		//vector<Motif> motifs;
		map<Size,Motif> motifs;
		AllOccs all_fg_occs, all_bg_occs;
		Reals all_pvals_adj(NUM_CORES, Real(NUM_CORES));

		for ( Size ind=0; ind<NUM_CORES; ++ind ) {
			if ( fg_counts[ind] < MIN_FG_COUNT ) continue;
			Size const fg_count(fg_counts[ind]), bg_count(bg_counts[ind]);
			Real const expected(max(PSEUDOCOUNT,Real(bg_count))*fg_bg_ratio),
				enrich(fg_count/expected);
			if (enrich>=MIN_CORE_ENRICH) {
				Real const pval(compute_overlap_pvalue(fg_count, fg_count+bg_count,
						fg_cdr3s.size(), fg_cdr3s.size() + bg_cdr3s.size())),
					pval_adj(pval*total_core_pvals);
				//Real chisq((fg_count-expected)*(fg_count-expected)/expected);
				//if (chisq>5){ // tmp hack
				if (pval_adj <= MAX_CORE_PVAL_ADJ) {
					cout << "core_pval: " << core_index_to_string(ind) <<
						//" pval: " << pval <<
						" pval_adj: " << pval_adj <<
						" enrich: " << enrich <<
						" counts: " << fg_count << ' ' << bg_count <<
						" total_pvals: " << total_core_pvals << endl;
					sortl.push_back(make_pair(pval*total_core_pvals, ind));
					motifs[ind] = core_index_to_motif(ind);
					all_fg_occs[ind]; // create empty list
					all_bg_occs[ind]; // create empty list
					all_pvals_adj[ind] = pval_adj;
				}
			}
		}

		// get occurrence list for all good cores
		get_core_occurrences_for_cdr3s(fg_cdr3s, all_fg_occs);
		get_core_occurrences_for_cdr3s(bg_cdr3s, all_bg_occs);


		sort(sortl.begin(), sortl.end());

		//Size a,b,c,d;
		Sizes fg_overlap(fg_cdr3s.size()), bg_overlap(bg_cdr3s.size());


		Size total_pair_pvals(0);
		vector<pair<Real, Size> > pair_sortl;
		map<Size, Sizes> all_combos; // indexed by "index"

		for (Size jj=1; jj<sortl.size(); ++jj) { // jj is the lower-scoring motif
			Size const ind2(sortl[jj].second);
			Motif const & m2(motifs[ind2]);
			Sizes const & m2_fg_occs(all_fg_occs.find(ind2)->second);
			Sizes const & m2_bg_occs(all_bg_occs.find(ind2)->second);
			for (Size ii=0; ii<jj; ++ii) { // ii is the higher-scoring motif
				Size const ind1(sortl[ii].second);
				Motif const & m1(motifs[ind1]);
				if (motifs_overlap(m1,m2)) {
					// look at intersection
					Sizes const & m1_fg_occs(all_fg_occs.find(ind1)->second);
					Sizes::iterator const fg_ite(
						set_intersection(m1_fg_occs.begin(), m1_fg_occs.end(),
							m2_fg_occs.begin(), m2_fg_occs.end(), fg_overlap.begin()));
					Size const fg_overlap_size(fg_ite - fg_overlap.begin());
					if ( fg_overlap_size >= MIN_FG_COUNT ){
						Sizes const & m1_bg_occs(all_bg_occs.find(ind1)->second);
						Sizes::iterator const bg_ite(
							set_intersection(m1_bg_occs.begin(), m1_bg_occs.end(),
								m2_bg_occs.begin(), m2_bg_occs.end(), bg_overlap.begin()));
						Size const bg_overlap_size(bg_ite - bg_overlap.begin());
						++total_pair_pvals;

						Real const
							expected(max(PSEUDOCOUNT,Real(bg_overlap_size))*fg_bg_ratio),
							enrich(fg_overlap_size/expected);
						if (enrich>=MIN_PAIR_ENRICH) {
							Real const
								pval(compute_overlap_pvalue(fg_overlap_size, fg_overlap_size+
										bg_overlap_size, fg_cdr3s.size(), fg_cdr3s.size() +
										bg_cdr3s.size())),
								pval_adj(pval*(total_core_pvals+total_pair_pvals));
							if (pval_adj <= MAX_PAIR_PVAL_ADJ) {
								cout << "pair_pval: " << core_index_to_string(ind1) << '+' <<
									core_index_to_string(ind2) << ' ' <<
									//" pval: " << pval <<
									" pval_adj: " << pval_adj <<
									" enrich: " << enrich <<
									" counts: " << fg_overlap_size << ' ' << bg_overlap_size <<
									" total_pvals: " << total_core_pvals + total_pair_pvals <<
									endl;
								// create a new index
								Size const index(fg_counts.size());
								// extend the arrays
								runtime_assert(index == bg_counts.size());
								fg_counts.push_back(fg_overlap_size);
								bg_counts.push_back(bg_overlap_size);
								all_pvals_adj.push_back(pval_adj);
								// add index to the dictionaries
								all_fg_occs[index].resize(fg_overlap_size);
								copy(fg_overlap.begin(), fg_ite,
									all_fg_occs.find(index)->second.begin());
								all_bg_occs[index].resize(bg_overlap_size);
								copy(bg_overlap.begin(), bg_ite,
									all_bg_occs.find(index)->second.begin());
								all_combos[index] = Sizes({ind1, ind2});
								pair_sortl.push_back(make_pair(pval_adj, index));
							}
						}
					}
				}
			}
		}

		sort(pair_sortl.begin(), pair_sortl.end());
		// now go through and look for pair+single combos
		//
		// what data structures do we have?
		// all_fg_occs, all_bg_occs  keys are "index"
		//
		// all_combos  keys are "index"
		//
		// motifs  keys are "index" only valid for core indices
		//
		Size total_triple_pvals(0);
		for ( Size ij_sum=0; ij_sum <pair_sortl.size()+sortl.size(); ++ij_sum) {
			for ( Size ii=0; ii<ij_sum && ii<pair_sortl.size(); ++ii) {
				Size const jj(ij_sum-ii);
				if (jj>=sortl.size()) continue;
				Size const ii_index(pair_sortl[ii].second);
				Size const jj_index(sortl[jj].second);
				Sizes const & ii_combos(all_combos[ii_index]);
				Motif const & m1(motifs[ii_combos[0]]);
				Motif const & m2(motifs[ii_combos[1]]);
				Motif const & m3(motifs[jj_index]);
				runtime_assert(motifs_overlap(m1,m2)); // HACKING REMOVE LATER
				if (motifs_overlap(m1,m3) || motifs_overlap(m2,m3)) {
					Sizes const & m1m2_fg_occs(all_fg_occs[ii_index]),
						& m1m2_bg_occs(all_bg_occs[ii_index]),
						& m3_fg_occs(all_fg_occs[jj_index]),
						& m3_bg_occs(all_bg_occs[jj_index]);

					Sizes::iterator const fg_ite(
						set_intersection(m1m2_fg_occs.begin(), m1m2_fg_occs.end(),
							m3_fg_occs.begin(), m3_fg_occs.end(), fg_overlap.begin()));
					Size const fg_overlap_size(fg_ite - fg_overlap.begin());
					if ( fg_overlap_size >= MIN_FG_COUNT ){
						Sizes::iterator const bg_ite(
							set_intersection(m1m2_bg_occs.begin(), m1m2_bg_occs.end(),
								m3_bg_occs.begin(), m3_bg_occs.end(), bg_overlap.begin()));
						Size const bg_overlap_size(bg_ite - bg_overlap.begin());
						++total_triple_pvals;

						Real const
							expected(max(PSEUDOCOUNT,Real(bg_overlap_size))*fg_bg_ratio),
							enrich(fg_overlap_size/expected);
						if (enrich>=MIN_TRIPLE_ENRICH) {
							Real const
								pval(compute_overlap_pvalue(fg_overlap_size, fg_overlap_size+
										bg_overlap_size, fg_cdr3s.size(), fg_cdr3s.size() +
										bg_cdr3s.size())),
								pval_adj(pval*(total_core_pvals+total_pair_pvals+
										total_triple_pvals));
							if (pval_adj <= MAX_TRIPLE_PVAL_ADJ) {
								cout << "triple_pval: " <<
									core_index_to_string(ii_combos[0]) << '+' <<
									core_index_to_string(ii_combos[1]) << '+' <<
									core_index_to_string(jj_index) << ' ' <<
									//" pval: " << pval <<
									" pval_adj: " << pval_adj <<
									" enrich: " << enrich <<
									" counts: " << fg_overlap_size << ' ' << bg_overlap_size <<
									" total_pvals: " << total_core_pvals + total_pair_pvals +
									total_triple_pvals << endl;

								// new motif, create a new index
								Size const index(fg_counts.size());
								// extend the arrays
								runtime_assert(index == bg_counts.size());
								fg_counts.push_back(fg_overlap_size);
								bg_counts.push_back(bg_overlap_size);
								all_pvals_adj.push_back(pval_adj);
								// add index to the dictionaries
								all_fg_occs[index].resize(fg_overlap_size);
								copy(fg_overlap.begin(), fg_ite,
									all_fg_occs.find(index)->second.begin());
								all_bg_occs[index].resize(bg_overlap_size);
								copy(bg_overlap.begin(), bg_ite,
									all_bg_occs.find(index)->second.begin());
								all_combos[index] = Sizes(
									{ii_combos[0], ii_combos[1], jj_index}); // pair then single
							} // pval <
						} // enrich >
					} // fg_overlap_size >
				}
			} // ii
		} // ij_sum, actually ii+jj


		// now sort all the motifs by pval, look for overlap
		vector<pair<Real,Size> > final_sortl;
		for (Size i=0; i<all_pvals_adj.size(); ++i ) {
			if (all_pvals_adj[i] <= MAX_FINAL_PVAL_ADJ) {
				final_sortl.push_back(make_pair(all_pvals_adj[i], i));
			}
		}
		sort(final_sortl.begin(), final_sortl.end());

		bools is_good(final_sortl.size(), true); // indexed by rank in final_sortl

		ofstream out_info(outfile_prefix+"_motif_info.tsv");
		ofstream out_members(outfile_prefix+"_motif_members.txt");
		ofstream out_distances(outfile_prefix+"_motif_distances.txt");

		out_info << "motif_string\tpval_adj\tenrich\tmax_jaccard\tnum_positions\tfg_count\tbg_count\ttotal_pvals\n";

		vector<Reals> all_jdists;
		Size const total_pvals(
			total_core_pvals + total_pair_pvals + total_triple_pvals);

		for ( Size ii=0; ii<final_sortl.size(); ++ii ) {
			Size const ii_index(final_sortl[ii].second);
			// check if too close to previous motif
			Real max_jaccard(0);
			Reals jdists;
			for ( Size jj=0; jj<ii; ++jj ) {
				if (is_good[jj]) {
					Size const jj_index(final_sortl[jj].second);
					Real const jaccard(compute_jaccard(
							all_fg_occs[ii_index], all_fg_occs[jj_index]));
					jdists.push_back(1.0-jaccard);
					max_jaccard = max(max_jaccard, jaccard);
					if ( jaccard > MAX_JACCARD ) {
						//cout<< "too close: " << jaccard << ' ' << ii << ' ' << jj << endl;
						is_good[ii] = false;
						break;
					}
				}
			}
			if ( is_good[ii] ) {
				// update the distance matrix
				jdists.push_back(0.0); // self dist
				for ( Size jj=0; jj<all_jdists.size(); ++jj ) { // make it square
					all_jdists[jj].push_back(jdists[jj]);
					// runtime_assert(all_jdists[jj].size() == jdists.size());
				}
				all_jdists.push_back(jdists);

				// not redundant by jaccard
				cout << "final_pval: ";
				// show the motif
				if (ii_index < NUM_CORES) {
					cout << motifs[ii_index].to_string();
					out_info << motifs[ii_index].to_string() << '\t';
				} else {
					bool first(true);
					for ( Size kk : all_combos[ii_index] ) {
						if (!first) {
							cout << '+';
							out_info << '+';
						} else first=false;
						cout << motifs[kk].to_string();
						out_info << motifs[kk].to_string();
					}
					out_info << '\t';
				}
				// this is not quite right: triple could only specify 4 positions...
				Size const num_positions
					(ii_index<NUM_CORES ? 3 : 2+all_combos[ii_index].size());
				Real const
					expected(max(PSEUDOCOUNT,Real(bg_counts[ii_index]))*fg_bg_ratio),
					enrich(fg_counts[ii_index]/expected);
				cout << " pval_adj: " << all_pvals_adj[ii_index] <<
					" enrich: " << enrich <<
					" max_jaccard: " << max_jaccard <<
					" num_positions: " << num_positions <<
					" counts: " << fg_counts[ii_index] << ' ' << bg_counts[ii_index] <<
					" total_pvals: " << total_pvals << endl;

				out_info << all_pvals_adj[ii_index] << '\t' <<
					enrich << '\t' <<
					max_jaccard << '\t' <<
					num_positions << '\t' <<
					fg_counts[ii_index] << '\t' <<
					bg_counts[ii_index] << '\t' <<
					total_pvals << '\n';
				{ // write members
					runtime_assert(all_fg_occs[ii_index].size() == fg_counts[ii_index]);
					bool first(true);
					for ( Size m : all_fg_occs[ii_index] ) {
						if (first) first=false;
						else out_members << ' ';
						out_members << m;
					}
					out_members << '\n';
				}
			} // good!
		}
		out_info.close();
		out_members.close();

		// write jaccard distance matrix:
		for ( Reals const & jdists : all_jdists ) {
			runtime_assert(jdists.size() == all_jdists.size());
			bool first(true);
			for ( Real d : jdists ) {
				if ( first ) first=false;
				else out_distances << ' ';
				out_distances << d;
			}
			out_distances << '\n';
		}
		out_distances.close();


		cout << "DONE" << endl;

	} catch (TCLAP::ArgException &e)  // catch any exceptions
		{
			std::cerr << "error: " << e.error() << " for arg " << e.argId() <<
				std::endl;
		}

}


