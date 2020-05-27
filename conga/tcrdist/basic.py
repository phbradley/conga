## LAZY-- some basic imports and shared params
##
##
from glob import glob
from os import popen, system, chdir, remove, getcwd, mkdir
from os.path import exists, isdir, isfile
import math
from math import floor,sqrt
from sys import stderr,argv,exit
import random
#from . import paths
#from blargs import Parser
#from parse_tsv import *

# yes, hardcoding for alpha beta right now
db_file = 'alphabeta_db.tsv'


## naming scheme for the gene segment types, occasionally useful for iterating
segtypes_uppercase = ['VA','JA','VB','JB']
segtypes_lowercase = ['va','ja','vb','jb']


############################################################################################
############################################################################################



## you could modify this function if you have a different cmdline tool for converting svg to png
## like inkscape or cairosvg
##
def convert_svg_to_png( svgfile, pngfile, verbose=True, allow_missing=False, allow_failure=True ):
    if not isfile(svgfile):
        errmsg = 'Error: convert_svg_to_png: svgfile does not exist: {}'.format(svgfile)
        print(errmsg)
        Log( errmsg )
        if allow_missing:
            return
        else:
            exit()
    cmd = 'convert {} {}'.format( svgfile, pngfile )
    if verbose:
        print(cmd)
    system(cmd)

    if isfile( pngfile ):
        ## success
        return

    ## cmdline inkscape
    cmd = 'inkscape --export-png {} {}'.format( pngfile, svgfile )
    if verbose:
        print(cmd)
    system(cmd)

    if isfile( pngfile ):
        ## success
        return

    ## this is probably a long-shot, but in case inkscape is installed on mac
    inkscape_exe = '/Applications/Inkscape.app/Contents/Resources/bin/inkscape'
    if isfile( inkscape_exe ):
        from os.path import abspath
        svgfile_full = abspath( svgfile )
        pngfile_full = abspath( pngfile )

        cmd = '{} --export-png {} {}'.format( inkscape_exe, pngfile_full, svgfile_full )
        if verbose:
            print(cmd)
        system(cmd)

        if isfile( pngfile ):
            ## success
            return


    ## another possibility
    cmd = 'rsvg-convert {} -o {}'.format( svgfile, pngfile )
    if verbose:
        print(cmd)
    system(cmd)

    if isfile( pngfile ):
        ## success
        return


    ## this might also occur if the svgfile were empty...
    errmsg = 'Error: convert command failed: cmd="{}" -- is the "convert" cmdline tool (Imagemagick) installed?'\
             .format( cmd )
    print(errmsg)
    Log( errmsg )
    if not allow_failure:
        exit()

# def get_mean_and_sdev( l ):
#     N = len(l)
#     assert N>0
#     mean = float( sum( l ) ) / N
#     sdev = 0.0
#     for x in l: sdev += ( x- mean )**2

#     sdev = sqrt( float(sdev)/N )
#     return mean, sdev

# def get_median(l_in):
#     l = l_in[:] ##make copy
#     l.sort()
#     n = len(l)
#     if n%2: return l[n/2]  ## n is odd
#     else: return 0.5 * ( l[n/2] + l[n/2 - 1 ] ) ## n is even


# def Log(s): ## silly legacy helper function
#     stderr.write(s)
#     if s and not s.endswith('\n'):
#         stderr.write('\n')


