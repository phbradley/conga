''' This script reads TCR information from a TSV file and makes alpha and beta chain
logos as PNG files
'''
import argparse
from os.path import exists
import sys
import os

required_cols = 'va ja cdr3a cdr3a_nucseq vb jb cdr3b cdr3b_nucseq'.split()

parser = argparse.ArgumentParser(
    description='Make TCR logos from a TSV file containing TCR gene and CDR3 information',
    epilog = f'''
    Required columns in the TSV file: {' '.join(required_cols)}

    Should generate SVG and PNG formatted files. Sometimes the default SVG-->PNG
    conversion process does not work well, in which case you could try converting the
    SVG files with something like Inkscape.

    Example first two lines from a TSV file:
va	ja	cdr3a	cdr3a_nucseq	vb	jb	cdr3b	cdr3b_nucseq
TRAV1-2*01	TRAJ33*01	CAVSDSNYQLIW	tgtgctgtgagtgatagcaactatcagttaatctgg	TRBV6-4*01	TRBJ2-1*01	CASSDGQPNNEQFF	tgtgccagcagtgatggacagcctaacaatgagcagttcttc

    Example command:

python3 ~/gitrepos/conga/scripts/make_tcr_logos.py --tcrs_tsvfile tcrs.tsv --outfile_prefix tcrs_test --organism human

    ''',
    formatter_class=argparse.RawDescriptionHelpFormatter,
)

parser.add_argument('--tcrs_tsvfile', required=True,
                    help='Tab-separated-values format file with TCR information, '
                    'required columns= '+' '.join(required_cols))
parser.add_argument('--outfile_prefix', required=True,
                    help='Prefix for the output SVG and PNG filenames')
parser.add_argument('--organism', choices=['mouse', 'human'], required=True)
#maybe these could work? choices=['mouse', 'human', 'mouse_gd', 'human_gd', 'human_ig')

args = parser.parse_args()

# debug the input file
import pandas as pd
tcrs_df = pd.read_table(args.tcrs_tsvfile)
for col in required_cols:
    assert col in tcrs_df, f'Need column {col} in {args.tcrs_tsvfile}'

# now the slower imports
conga_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if conga_dir not in sys.path:
    sys.path.append(conga_dir)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import conga
from conga.tcrdist.tcr_distances import TcrDistCalculator
from conga.tcrdist.make_tcr_logo import make_tcr_logo_for_tcrs

# the logo-making code expects a list of tuples, not a pandas dataframe:
tcrs = [((t.va, t.ja, t.cdr3a, t.cdr3a_nucseq), (t.vb, t.jb, t.cdr3b, t.cdr3b_nucseq))
        for t in tcrs_df.itertuples()]
print(f'Read {len(tcrs)} paired TCRs from {args.tcrs_tsvfile}')

tcrdist_calculator = TcrDistCalculator(args.organism)

for ab in 'AB':
    pngfile = f'{args.outfile_prefix}_tcr_logo_{ab}.png'
    make_tcr_logo_for_tcrs(
        tcrs, ab, args.organism, pngfile, tcrdist_calculator=tcrdist_calculator)
    print('made:', pngfile)

