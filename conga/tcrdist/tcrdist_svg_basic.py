from os import system
from os.path import exists
import base64 ## png encoding
import math
from . import basic ## convert_svg_to_png
from .basic import MONOSPACE_FONT_FAMILY
from .html_colors import CB_RED, CB_GREEN, CB_BLUE, CB_ORANGE, CB_PURPLE


## these are not actually defined by the language...
all_text_height = {20:15.5, 100:75.}
all_text_width  = {20:12.5, 100:60.}

def rgb_from_fraction(fraction, cmap=None):

    assert fraction >=0 and fraction<=1

    fraction = float(fraction)

    if cmap is None:
        if fraction<0.5:
            i = 2*fraction*256

            #green = min(200,2*i)
            green = int( min(200,(200*i)/128) )
            blue = max(0,int( min(255,510-2*i) ))
            red = 0
        else:
            i = 2*(fraction-0.5)*256
            red = int( min(255,2*i) )
            green = int( max(0,min(200,400 - (200*i)/128)) )
            #green = min(200,510-2*i)
            blue = 0
    else:
        color_tuple = cmap(fraction)[:3]
        red, green, blue = [min(255,int(x*255)) for x in color_tuple]

    assert 0<=red  <=255
    assert 0<=green<=255
    assert 0<=blue <=255

    return '#{:02x}{:02x}{:02x}'.format(red,green,blue)


    # return '<rect x="{}" y="{}" height="{}" width="{}" fill="{}" stroke="{}"/>\n'\
    #     .format( upper_left[0], upper_left[1], lower_right[1]-upper_left[1], lower_right[0]-upper_left[0],
    #              fill, stroke )

def rectangle( upper_left, lower_right, fill, stroke, stroke_width=1, dashed=False ):
    stroke_dasharray_style= ';stroke-dasharray:5,5' if dashed else ''
    return '<rect x="{}" y="{}" height="{}" width="{}" style="fill:{};stroke:{};stroke-width:{}{}" />\n'\
        .format( upper_left[0], upper_left[1], lower_right[1]-upper_left[1], lower_right[0]-upper_left[0],
                 fill, stroke, stroke_width, stroke_dasharray_style )

def create_file( cmds, width, height, filename, create_png = False, background_color=None, use_xlink=False ):
    out = open( filename,'w')
    extra = '' if not use_xlink else 'xmlns:xlink="http://www.w3.org/1999/xlink"'
    out.write('<svg width="{}" height="{}" xmlns="http://www.w3.org/2000/svg" version="1.1" {} >\n'\
              .format(int(width),int(height),extra))
    if background_color:
        out.write(rectangle( (0,0), (width,height), background_color, 'white', 0 ) )
    out.write('\n'.join(cmds)+'\n')
    out.write('</svg>\n')
    out.close()

    if create_png:
        assert filename.endswith('.svg')
        pngfile = filename[:-4]+'.png'
        basic.convert_svg_to_png( filename, pngfile )


## this will work if oldfile was created using create_file
## returns: cmds, old_width, old_height
## cmds do not include the newline from the original file
## ie, they could go right into create_file
##
## matplotlib seems to generate files that are sized in pt units, which is a bit of a nuisance
## converting between them is something like #px = 1.25 * #pt
## pt = 1/72 inch, so inch = 72pt = 72*1.25=72*5/4 = 18*5 = 90px
## but somewhere on the web somebody says the new standard is 1inch=96px...
## at least that's what I get when I run convert on a matplotlib generated svg to make a png
## ( 576pt x 864pt goes to 720px x 1080px )
##
def embed_file( oldfile, x_shift, y_shift ):
    cmds = [] ##
    cmds = ['<g transform="translate({},{})" >'.format(x_shift,y_shift) ]
    for line in open(oldfile,'r'):
        if line.startswith("<svg"):
            l = line.split()
            assert l[1].startswith('width="')
            assert l[2].startswith('height="')
            old_width  = float( l[1].split('"')[1] ) ## will fail if pt or px or something in size
            old_height = float( l[2].split('"')[1] )
        elif line.startswith('</svg'):
            break
        else:
            cmds.append(line[:-1])
    cmds.append( '</g>')
    return cmds, old_width, old_height

## returns a single line
## note that we have to set use_xlink=True when we call create_file in order for the link to work here:
##
def embed_pngfile( pngfile, width, height, x0, y0, aspect=None ):
    assert exists( pngfile )
    encoded = base64.b64encode(open(pngfile, "rb").read()).decode('ascii')

    asp = '' if aspect==None else 'preserveAspectRatio="{}"'.format(aspect)

    cmd = '<image width="{}" height="{}" x="{}" y="{}" {} xlink:href="data:image/png;base64,{}" />'\
                     .format( width, height, x0, y0, asp, encoded )
    return cmd


def make_text( text, lower_left, fontsize,
               extra_tag = None, font_family = MONOSPACE_FONT_FAMILY, color = "black", font_weight = "normal" ):
    assert font_weight in ['normal','bold']
    cmd = '<text x="{:.3f}" y="{:.3f}" font-size="{}" font-weight="{}" font-family="{}" fill="{}" xml:space="preserve">{}</text>\n'\
        .format( lower_left[0], lower_left[1], fontsize, font_weight, font_family, color, text )
    return cmd


class SVG_tree_plotter:
    def __init__( self, cmap=None ):
        self.cmds = []
        self.cmap = cmap

    def make_line( self, p0, p1, line_width, normalized_score, extra_tag=None, color=None ):

        if normalized_score==None:
            #color='black'
            color = '#aaaaaa' # gray
        elif color==None:
            color = rgb_from_fraction(normalized_score, cmap=self.cmap)

        if p0[0] == p1[0]: ## vertical line, get width exactly right
            x = p0[0] ; y0 = min(p0[1],p1[1]) ; y1 = max(p0[1],p1[1] )

            upper_left = ( x-line_width/2.0, y0 )
            lower_right = ( x+line_width/2.0, y1 )
            cmd = rectangle( upper_left, lower_right, color, color )

        elif p0[1] == p1[1]: ## horizontal line, get width exactly right
            y = p0[1] ; x0 = min(p0[0],p1[0]) ; x1 = max(p0[0],p1[0] )

            upper_left  = ( x0, y-line_width/2.0 )
            lower_right = ( x1, y+line_width/2.0 )
            cmd = rectangle( upper_left, lower_right, color, color )
        else:
            cmd = '<line x1="{:.3f}" y1="{:.3f}" x2="{:.3f}" y2="{:.3f}" stroke-width="{}" stroke="{}"/>\n'\
                .format( p0[0], p0[1], p1[0], p1[1], 0.5*line_width,color )
        self.cmds.append( cmd )

    def make_text( self, text, p, fontsize, extra_tag=None ):
        color = 'black'
        ## shift so that p is the lower right
        textwidth = fontsize*0.6 * len(text)
        cmd = '<text x="{:.3f}" y="{:.3f}" font-size="{}" font-family="{}" fill="{}" xml:space="preserve">{}</text>\n'\
            .format( p[0]-textwidth, p[1], fontsize, MONOSPACE_FONT_FAMILY, color, text )
        self.cmds.append( cmd )

    def write( self, out ):
        out.write('\n'.join( self.cmds )+'\n' )



def color_stack( upper_left, lower_right, letters, colors, values ):
    fontsize= 20
    fontwidth  = all_text_width [fontsize]
    fontheight = all_text_height[fontsize]
    width  = lower_right[0] - upper_left[0]
    height = lower_right[1] - upper_left[1]

    ## draw them proportional to their values
    total = sum(values)

    ## draw from the top down
    height_sum = 0.
    lines = []

    for (letter,(color,value)) in zip( letters, list(zip( colors, values )) ):
        my_height = height * float(value)/total
        height_sum += my_height
        my_x = upper_left[0]
        my_y = upper_left[1] + height_sum

        ## we are going to have to scale
        x_scale = float( width )     / fontwidth
        y_scale = float( my_height ) / fontheight

        #print letter, 'my_height:',my_height, 'y_scale:',y_scale,'fontheight:',fontheight,'my_xy:',my_x,my_y

        ## we need to translate so that
        lines.append( '<text x="{:.6f}" y="{:.6f}" font-size="{}" font-family="{}" fill="{}" transform="scale({:.6f},{:.6f})">{}</text>\n'\
                      .format( my_x / x_scale, my_y / y_scale, fontsize, MONOSPACE_FONT_FAMILY, color,
                               x_scale, y_scale, letter ) )

        #lines.append( '<text x="{:.6f}" y="{:.6f}" font-size="{}" font-family="monospace" fill="{}" >{}</text>\n'\
        #              .format( my_x , my_y , fontsize, color, letter ) )
    return ''.join( lines )


def text_in_box( upper_left, lower_right, text, color ):
    fontsize= 100

    fontwidth  = all_text_width [fontsize]
    fontheight = all_text_height[fontsize]

    width  = lower_right[0] - upper_left[0]
    height = lower_right[1] - upper_left[1]

    ## we are going to have to scale
    x_scale = float( width )     / (fontwidth * len(text))
    y_scale = float( height ) / fontheight

    #print text, 'height:',height, 'y_scale:',y_scale,'fontheight:',fontheight,'my_xy:'

    ## we need to translate so that
    return '<text x="{:.6f}" y="{:.6f}" font-size="{}" font-family="{}" fill="{}" transform="scale({:.6f},{:.6f})">{}</text>\n'\
        .format( float(upper_left [0]) / x_scale,
                 float(lower_right[1]) / y_scale,
                 fontsize, MONOSPACE_FONT_FAMILY, color, x_scale, y_scale, text )



def protein_logo( upper_left, lower_right, pwm, scale={} ):## scale[pos] should be in range [0,1]

    aacolor = {}
    for aa in "G"       :aacolor[aa] = CB_ORANGE #"orange" ## special for glycine
    for aa in "STYC"    :aacolor[aa] = CB_GREEN #"green"
    for aa in "NQ"      :aacolor[aa] = CB_PURPLE #"purple"
    for aa in "KRH"     :aacolor[aa] = CB_BLUE #"blue"
    for aa in "DE"      :aacolor[aa] = CB_RED #"red"
    for aa in "P"       :aacolor[aa] = "black" ## for right now
    for aa in "AWFLIMV" :aacolor[aa] = "black" ## do we want to make W,F something different??
    for aa in "J":aacolor[aa] = "green"
    for aa in ".-":aacolor[aa] = "black"

    width  = lower_right[0] - upper_left[0]
    height = lower_right[1] - upper_left[1]

    N = len(list(pwm.keys()))

    letter_width = float(width) / N

    cmds = []

    for pos in range(N):
        l = [(y,x) for x,y in pwm[pos].items() ]
        l.sort()
        l.reverse()
        ## start from the top
        col_scale = scale.get(pos,1.0)
        assert col_scale>=0 and col_scale<=1.00001
        scaled_height = col_scale * height
        upper_y = lower_right[1] - scaled_height

        totfreq=0.0
        for freq,aa in l:
            if freq * scaled_height>1: ## at least one pixel high
                y0 = upper_y + totfreq * scaled_height
                y1 = upper_y + (totfreq+freq) * scaled_height
                x0 = upper_left[0] + pos*letter_width
                x1 = upper_left[0] + (pos+1)*letter_width
                cmds.append( text_in_box( (x0,y0), (x1,y1), aa, aacolor[aa] ) )
            totfreq += freq
            #print pos,aa,freq,totfreq,sum(pwm[pos].values()),pwm[pos]
        assert abs(totfreq-1.0)<1e-2

    return '\n'.join(cmds)+'\n'

def generic_logo( upper_left, lower_right, pwm ):

    color = 'black'

    width  = lower_right[0] - upper_left[0]
    height = lower_right[1] - upper_left[1]

    N = len(list(pwm.keys()))

    column_width = float(width) / N

    cmds = []

    for pos in range(N):
        l = [(y,x) for x,y in pwm[pos].items() ]
        l.sort()
        l.reverse()
        ## start from the top
        totfreq=0
        for freq,aa in l:
            if freq<1e-3:continue
            y0 = upper_left[1] + totfreq * height
            y1 = upper_left[1] + (totfreq+freq) * height
            x0 = upper_left[0] +  pos   *column_width
            x1 = upper_left[0] + (pos+1)*column_width
            cmds.append( text_in_box( (x0,y0), (x1,y1), aa, color ) )
            totfreq += freq
    return '\n'.join(cmds)+'\n'


def make_stack( upper_left, lower_right, l ):
    x0 = upper_left [0]
    x1 = lower_right[0]

    width  = lower_right[0] - upper_left[0]
    height = lower_right[1] - upper_left[1]

    l.sort()
    l.reverse()

    total = sum([x[0] for x in l])
    totval = 0
    cmds = []
    for vwc in l:
        if len(vwc) == 2:
            val,word = vwc
            color = 'black'
        else:
            assert len(vwc) == 3
            val,word,color = vwc
        y0 = upper_left[1] + height * float(totval) / total
        y1 = upper_left[1] + height * float(totval+val) / total
        cmds.append( text_in_box( (x0,y0), (x1,y1), word, color ) )
        totval+=val

    return '\n'.join(cmds)+'\n'



def enrichment_glyph_marker_old(marker_id):
    cmd = """<marker id="{}"
      viewBox="0 0 10 10" refX="10" refY="5"
      markerUnits="strokeWidth"
      markerWidth="4" markerHeight="3"
      orient="auto">
      <path d="M 0 0 L 10 5 L 0 10 z" />
    </marker>
    """.format(marker_id)
    return cmd

def enrichment_glyph_old( center, arrow_length, arrow_width, fontsize, marker_id, enrichment ):
    line_text_sep = 2.
    if enrichment>1: ## up arrow, text at the top
        line_p0 = [ center[0], center[1] + arrow_length/2. ]
        line_p1 = [ center[0], center[1] - arrow_length/2. ]
        text_lower_y = line_p1[1] - line_text_sep
        enrichment_text = '{:.1f}'.format(enrichment)
    else:
        line_p0 = [ center[0], center[1] - arrow_length/2. ]
        line_p1 = [ center[0], center[1] + arrow_length/2. ]
        text_lower_y = line_p1[1] + line_text_sep + 0.75 * fontsize
        enrichment_text = '{:.1f}'.format(1.0/enrichment)
    text_lower_x = center[0] - 0.6 * fontsize * 0.5 * len(enrichment_text)



    ## define the marker
    cmds = []
    cmds.append(
        '<line x1="{:.3f}" y1="{:.3f}" x2="{:.3f}" y2="{:.3f}" marker-end="url(#{})" stroke="black" stroke-width="{}"/>'\
        .format( line_p0[0], line_p0[1], line_p1[0], line_p1[1], marker_id, arrow_width ) )
    cmds.append( make_text( enrichment_text, [ text_lower_x, text_lower_y ], fontsize ) )

    return cmds

# eg_marker_id_prefix = 'eg_'
# eg_max_arrow_heads = 6

# def enrichment_glyph_markers():
#     global eg_marker_id_prefix
#     global eg_max_arrow_heads

#     for heads in range(1,eg_max_arrow_heads+1):
#         markerid = '{}{}'.format( eg_marker_id_prefix, heads )

#     cmd = """<marker id="{}"
#       viewBox="0 0 10 10" refX="10" refY="5"
#       markerUnits="strokeWidth"
#       markerWidth="4" markerHeight="3"
#       orient="auto">
#       <path d="M 0 0 L 10 5 L 0 10 z" />
#     </marker>
#     """.format(marker_id)
#     return cmd

def enrichment_glyph_cmds( center, arrow_length, arrow_width, enrichment, add_rectangle = False, box_color = 'gold' ):
    cmds = []

    num_heads = int(math.floor( abs( math.log(enrichment,2.) ) ))
    if num_heads<1:
        return cmds

    head_sep = 3. * arrow_width
    head_width = 3. * arrow_width
    head_slant = 1.

    if enrichment>1: ## up arrow, text at the top
        line_p0 = [ center[0], center[1] + arrow_length/2. ]
        line_p1 = [ center[0], center[1] - arrow_length/2. ]
        head_step = 1 * head_sep
    else:
        line_p0 = [ center[0], center[1] - arrow_length/2. ]
        line_p1 = [ center[0], center[1] + arrow_length/2. ]
        head_step = -1 * head_sep

    if add_rectangle:
        ## put a nice rounded rectangle around the outside
        rect_x0 = center[0] - head_width - arrow_width
        rect_y0 = center[1] - arrow_length/2 - arrow_width
        rect_width = 2*(head_width+arrow_width)
        rect_height = arrow_length + 2*arrow_width
        rect_round = 3.0 * arrow_width

        cmds.append( '<rect x="{:.3f}" y="{:.3f}" width="{:.3f}" height="{:.3f}" rx="{:.3f}" ry="{:.3f}" stroke="{}" fill="{}"/>'\
                     .format( rect_x0, rect_y0, rect_width, rect_height, rect_round, rect_round,
                              box_color, box_color ))


    ## make the line
    cmds.append( '<line x1="{:.3f}" y1="{:.3f}" x2="{:.3f}" y2="{:.3f}" stroke="black" stroke-width="{}"/>'\
                 .format( line_p0[0], line_p0[1], line_p1[0], line_p1[1], arrow_width ) )

    ## now make the heads

    for head in range(num_heads):
        for xsign in [1,-1]:
            x1 = line_p0[0]
            x2 = x1 + xsign * head_width
            y1 = line_p1[1] + head * head_step
            y2 = y1 + head_slant * head_step

            cmds.append( '<line x1="{:.3f}" y1="{:.3f}" x2="{:.3f}" y2="{:.3f}" stroke="black" stroke-width="{}"/>'\
                         .format( x1,y1,x2,y2, arrow_width ) )

            if False: ## add shifted white line
                y_shift = arrow_width * head_step / head_sep
                cmds.append( '<line x1="{:.3f}" y1="{:.3f}" x2="{:.3f}" y2="{:.3f}" stroke="white" stroke-width="{}"/>'\
                             .format( x1,y1+y_shift,x2,y2+y_shift, arrow_width ) )


    return cmds


if __name__ == '__main__':
    cmds = []
    # marker_id = 'tri'
    # cmds.append( enrichment_glyph_marker( marker_id ) )
    center = [50,50]
    arrow_length = 40.
    arrow_width = 3.
    enrichment = 2.01
    cmds += enrichment_glyph_cmds( center, arrow_length, arrow_width, enrichment )

    center = [100,50]
    enrichment = 45.01
    cmds += enrichment_glyph_cmds( center, arrow_length, arrow_width, enrichment )

    center = [250,50]
    enrichment = 0.49
    cmds += enrichment_glyph_cmds( center, arrow_length, arrow_width, enrichment )

    svgfile = 'tmp.svg'
    width = 300
    height = 100

    create_file( cmds, width, height, svgfile, create_png=True, background_color='red' )





    # out = open('tmp3.svg','w')

    # out.write( '<svg width="1000" height="1000">\n' )

    # texts = [ ( (0,0), (500,300), 'DUDE' ),\
    #               ( (100,350), (700,400), 'TEST' ) ]


    # for (upper_left,lower_right,text) in texts:


    #     out.write( '<rect x="{}" y="{}" height="{}" width="{}" fill="none" stroke="black"/>\n'\
    #                    .format( upper_left[0], upper_left[1], lower_right[1]-upper_left[1], lower_right[0]-upper_left[0] ) )

    #     out.write( text_in_box( upper_left, lower_right, text ) )


    # out.write( '</svg>\n' )
