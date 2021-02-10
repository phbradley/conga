from os import system
from os.path import isfile
from sys import stderr, exit

# modify this (ie PATH_TO_INKSCAPE) if you have command line inkscape installed
#   in a different place
# note that we also try just 'inkscape' first so if it's in your path you should be good.
#
PATH_TO_INKSCAPE='/Applications/Inkscape.app/Contents/MacOS/inkscape'
ALT_PATH_TO_INKSCAPE='/Applications/Inkscape.app/Contents/Resources/bin/inkscape'

# if the fonts in the tcr logos look bad, try uncommenting one of the other three lines below this one
#MONOSPACE_FONT_FAMILY = 'monospace'
#MONOSPACE_FONT_FAMILY = 'courier'
MONOSPACE_FONT_FAMILY = 'DejaVu Sans Mono'

if isfile(ALT_PATH_TO_INKSCAPE) and not isfile(PATH_TO_INKSCAPE):
    PATH_TO_INKSCAPE = ALT_PATH_TO_INKSCAPE

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
    from shutil import which

    if not isfile(svgfile):
        errmsg = 'Error: convert_svg_to_png: svgfile does not exist: {}\n'.format(svgfile)
        print(errmsg)
        stderr.write( errmsg )
        if allow_missing:
            return
        else:
            exit()

    if which('convert') is not None:
        cmd = 'convert {} {}'.format( svgfile, pngfile )
        if verbose:
            print(cmd)
        system(cmd)

        if isfile( pngfile ):
            ## success
            return
        else:
            print(f'conga.convert_svg_to_png: convert failed, command="{cmd}"')

    if which('inkscape') is not None:
        ## cmdline inkscape: new options
        cmd = 'inkscape -o {} {}'.format( pngfile, svgfile )
        if verbose:
            print(cmd)
        system(cmd)

        if isfile( pngfile ):
            ## success
            return
        else:
            print(f'conga.convert_svg_to_png: inkscape failed, command="{cmd}"')

        cmd = 'inkscape --export-png {} {}'.format( pngfile, svgfile )
        if verbose:
            print(cmd)
        system(cmd)

        if isfile( pngfile ):
            ## success
            return
        else:
            print(f'conga.convert_svg_to_png: inkscape failed, command="{cmd}"')

    ## this is probably a long-shot, but in case inkscape is installed on mac
    if isfile( PATH_TO_INKSCAPE ):
        from os.path import abspath
        svgfile_full = abspath( svgfile )
        pngfile_full = abspath( pngfile )

        # this is the new syntax
        cmd = '{} -o {} {}'.format( PATH_TO_INKSCAPE, pngfile_full, svgfile_full )
        if verbose:
            print(cmd)
        system(cmd)

        if isfile( pngfile ):
            ## success
            return

        # this is the old syntax
        cmd = '{} --export-png {} {}'.format( PATH_TO_INKSCAPE, pngfile_full, svgfile_full )
        if verbose:
            print(cmd)
        system(cmd)

        if isfile( pngfile ):
            ## success
            return
        else:
            print(f'conga.convert_svg_to_png: inkscape failed, command="{cmd}"')


    if which('rsvg-convert') is not None:
        ## another possibility
        cmd = 'rsvg-convert {} -o {}'.format( svgfile, pngfile )
        if verbose:
            print(cmd)
        system(cmd)

        if isfile( pngfile ):
            ## success
            return
        else:
            print(f'conga.convert_svg_to_png: 3rd try failed, command="{cmd}"')

    if which('cairosvg') is not None:
        ## another possibility
        cmd = f'cairosvg -f png -o {pngfile} {svgfile}'
        if verbose:
            print(cmd)
        system(cmd)

        if isfile( pngfile ):
            ## success
            return
        else:
            print(f'conga.convert_svg_to_png: 3rd try failed, command="{cmd}"')

    if which('magick') is not None:
        ## another possibility
        cmd = 'magick convert {} {}'.format( svgfile, pngfile )
        if verbose:
            print(cmd)
        system(cmd)

        if isfile( pngfile ):
            ## success
            return
        else:
            print(f'conga.convert_svg_to_png: 4th try failed, command="{cmd}"')


    ## this might also occur if the svgfile were empty...
    errmsg = 'Error: conga.convert_svg_to_png failed to convert svg file to png file\nIs the "convert" cmdline tool (ImageMagick) installed, or Inkscape?\n'
    print(errmsg)
    stderr.write( errmsg )
    if not allow_failure:
        exit()
    return
