import argparse

parser = argparse.ArgumentParser(description="""Generate tcrdist database files for the tcrdist_cpp code.
Only necessary to re-run if the supported organisms change or new genes are added to the tcrdist database
in <path_to_conga>/conga/tcrdist/db/combo_xcr.tsv""")

#type is str by default
parser.add_argument('--outdir', required=True, help="""The directory in which to put the database files. Should
probably be something like <path_to_conga>/tcrdist_cpp/db/ since that is where the conga code will look.""")
args = parser.parse_args()

if args.outdir[-1] != '/':
    args.outdir += '/'

import sys
import os
conga_dir = os.path.dirname( os.path.dirname( os.path.abspath(__file__) ) )
sys.path.append(conga_dir) # in order to import conga package
import conga
from conga.tcrdist.tcr_distances import compute_all_v_region_distances
from conga.tcrdist.tcr_distances_blosum import bsd4
from conga.tcrdist.all_genes import all_genes
from conga.tcrdist.amino_acids import amino_acids

for organism in all_genes:
    outfile = '{}tcrdist_info_{}.txt'.format(args.outdir, organism)
    print('making:', outfile)

    out = open(outfile, 'w')

    for aa in amino_acids:
        out.write('AAdist {} {}\n'.format(aa, ' '.join(['{:.3f}'.format( bsd4[(aa,x)]) for x in amino_acids])))

    rep_dists = compute_all_v_region_distances(organism)
    for ab in 'AB':

        #v_genes = [ x for x in rep_dists.keys() if x[2] == ab ]
        v_genes = sorted([ x for x,y in all_genes[organism].items() if y.chain == ab and y.region == 'V' ])
        print('num V{} genes {}'.format(ab, len(v_genes)))

        out.write('num_V{}_genes {}\n'.format( ab, len(v_genes) ))

        for v1 in v_genes:
            out.write('V{}dist {} {}\n'.format(ab, v1, ' '.join(['{:.3f}'.format(rep_dists[v1][x]) for x in v_genes])))
    out.close()
