from os import system
from os.path import isfile
from sys import stderr, exit

# modify this if you have command line inkscape installed in a different place
# note that we also try just 'inkscape' first so if it's in your path you should be good.
#
PATH_TO_INKSCAPE='/Applications/Inkscape.app/Contents/Resources/bin/inkscape'


## you could modify this function if you have a different cmdline tool for converting svg to png
## like cairosvg
##

def convert_svg_to_png(
        svgfile,
        pngfile,
        verbose=False,
        allow_missing=False,
        allow_failure=True
):
    if not isfile(svgfile):
        errmsg = 'Error: convert_svg_to_png: svgfile does not exist: {}\n'.format(svgfile)
        print(errmsg)
        stderr.write( errmsg )
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
    if isfile( PATH_TO_INKSCAPE ):
        from os.path import abspath
        svgfile_full = abspath( svgfile )
        pngfile_full = abspath( pngfile )

        cmd = '{} --export-png {} {}'.format( PATH_TO_INKSCAPE, pngfile_full, svgfile_full )
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
    errmsg = 'Error: convert command failed: cmd="{}" -- is the "convert" cmdline tool (Imagemagick) installed?\n'\
        .format( cmd )
    print(errmsg)
    stderr.write( errmsg )
    if not allow_failure:
        exit()
    return

