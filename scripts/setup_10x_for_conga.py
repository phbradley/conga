import argparse
from os.path import exists
import sys
import os

parser = argparse.ArgumentParser()

parser.add_argument('--clones_file')
parser.add_argument('--organism', choices=['mouse', 'human'], required=True)
#parser.add_argument('--n_components', type=int, default=50)
parser.add_argument('--filtered_contig_annotations_csvfile', required=True)
parser.add_argument('--consensus_annotations_csvfile')

args = parser.parse_args()
assert exists(args.filtered_contig_annotations_csvfile)
if args.consensus_annotations_csvfile is not None:
    assert exists(args.consensus_annotations_csvfile)

# put this after arg parsing because it's so dang slow
sys.path.append( os.path.dirname( os.path.dirname( os.path.abspath(__file__) ) ) ) # so we can import conga
import conga
from compute_tcrdist_kpcs_from_clones_file import compute_tcrdist_kpcs_from_clones_file

##################################### END IMPORTS AND ARG PARSING #########################################


default_clones_file = args.filtered_contig_annotations_csvfile+'_tcrdist_clones.tsv'

#### first make the clones file ###########################################################################
python2 = conga.util.PYTHON2_EXE
script = conga.util.TCRDIST_REPO+'make_10x_clones_file.py'
if not exists(script):
    print('cant find the tcr-dist script:', script)
    print('please install tcr-dist repo from github and set the location in conga/conga/util.py at the top')
    sys.exit(1)

extra_args = ''
if args.consensus_annotations_csvfile:
    assert exists(args.consensus_annotations_csvfile)
    extra_args += f' --consensus_annotations_csvfile {args.consensus_annotations_csvfile} '

clones_file = default_clones_file if args.clones_file is None else args.clones_file
cmd = '{} {} --organism {} --clones_file {} --filtered_contig_annotations_csvfile {} --clobber --stringent {}'\
    .format( python2, script, args.organism, clones_file, args.filtered_contig_annotations_csvfile, extra_args )
print('run:', cmd)
sys.stdout.flush()
conga.util.run_command(cmd)

if not exists(clones_file):
    print('failed to make the clones_file!\ncmd= {}'.format(cmd))
    sys.exit(1)

#### Now compute the kernal PCs #############################################################3
compute_tcrdist_kpcs_from_clones_file( clones_file, args.organism )

print(f'If this all worked you should be able to pass {clones_file} as the --clones_file argument to run_conga.py')
