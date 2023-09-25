#!/usr/bin/env python
"""
This code belongs to the repositories:
https://gitlab.hpc.cineca.it/ceottogroup/draw_molecular_geometries

Author: Michele Gandolfi
Date:   March 2022

Program to draw simple molecules with tex + asymtote
"""

# standard packages
import os
import argparse
from argparse import RawTextHelpFormatter

# extra required packages for efficiency
import numpy  as np
import pandas as pd

# Michele's custom packages
import vmd


# INITIALIZE PARENT ARGUMENT PARSER
parser = argparse.ArgumentParser(
          description="""
                      Draw a molecular structure with asymptote.
                      """
         #,formatter_class= argparse.ArgumentDefaultsHelpFormatter
         ,formatter_class=RawTextHelpFormatter
         ,fromfile_prefix_chars='+'
         ,epilog='Further arguments can be given in a file using +FILENAME'
         )

# GROUP GLOBAL ARGUMENTS
global_args = parser.add_argument_group( title='global arguments', description=None)
global_args.add_argument( 'xyz', help='input xyz file geometry')
global_args.add_argument( '--debug', action='store_true', help='enter program in debug mode')
global_args.add_argument( '-v', '--verbose', action='count', default=0, help='write output verbosely (repeat to increase verbosity, e.g. -vvv)')
global_args.add_argument( '-O', '--output', default='picture.tex', help='name of the output file')

# GROUP STYLE ARGUMENTS
style_args = parser.add_argument_group( title='style arguments', description=None)
#style_args.add_argument( '--color_scheme', default=None, help='filename for custom color scheme')
style_args.add_argument( '--sph_s', default='0.30', help='atom sphere scale')
style_args.add_argument( '--cyl_s', default='0.08', help='bond cylinder radius scale')
style_args.add_argument( '--proj_type', default='orthographic', help='type of projection')
style_args.add_argument( '--proj', nargs=3, default=[0,0,0], help='projection')

topo_args = parser.add_argument_group( title='topology arguments', description=None)
topo_args.add_argument( '-T', '--topo', default=None, help="""topology space separated value file with predefined headers. As:
id1 id2 bond
1   2   [val12]
1   3   [val13]
2   3   [val23]""")
topo_args.add_argument( '--topo_sep', default='\s+', help="character (or regex) separator of topology file")
topo_args.add_argument( '--topo_title', action='store_true', help="topology file has (1) header line")
topo_args.add_argument( '--numlab', action='store_true', help='display label on bond')
topo_args.add_argument( '--colorlab', action='store_true', help='color bond with grayscale depend on label')
topo_args.add_argument( '--sizelab', action='store_true', help='display larger bonds for large labels')
topo_args.add_argument( '--maxval', default=8, type=int, help='maximum bin id for colorlab and sizelab')
topo_args.add_argument( '--num_bins', default=15, type=int, help='number of bins for colorlab and sizelab')
topo_args.add_argument( '--scalefact', default=1.e6, type=float, help='scale factor for bonds labels')
topo_args.add_argument( '--labelsize', default=2, type=int, help='label font size in px (usually 2 or 3 is fine)')
topo_args.add_argument( '--labeldigits', default=1, type=int, help='number of digits displayed for the labels')
topo_args.add_argument( '--whichatomlabel', nargs='+', default=["all",], type=int, help='indexes of atoms whose connections are labeled (if labels are prescribed)')

render_args = parser.add_argument_group( title='rendering arguments', description=None)
render_args.add_argument( '--render', default=4, type=int, help='rendering value (higher for more quality). If ftype=pdf, then render should be 0, else 5 should be way more than enough')
render_args.add_argument( '--ftype' , default='png', choices=['png', 'pdf'], help='format of embedded figure')

###############################################################################

# PARSE ARGUMENTS FROM STDIN
args = parser.parse_args()

if args.topo is not None:
    assert os.path.isfile( args.topo), f'No such file {args.topo}. Aborting'
    header = 1 if args.topo_title else 0
    df = pd.read_csv( args.topo, comment='#', sep=args.topo_sep, header=header)
    assert df.shape[1] in {2,3}, 'Something wrong with topology file. It should have 2 or 3 named columns'
    if df.shape[1] == 2:
        df['bond'] = args.bond_s
    # ensure correct column names
    df = df.set_axis( ['id1' ,'id2', 'bond'], axis='columns')
    df.astype( {'id1': 'str', 'id2':'str', 'bond':'float64'})
    #if df['bond'].max() < 1.0:
    #    df['bond'] = ( df['bond'] - df['bond'].min() ) / ( df['bond'].max() - df['bond'].min() ) + 0.05

    unique_bonds = set( df['bond'])

    if args.whichatomlabel[0] == 'all':
        args.whichatomlabel = list( range( 0, max(
                                                    max( df['id1']),
                                                    max( df['id2'])
                                                 )))
    else:
        # reduce index of whichatomlabel by 1
        args.whichatomlabel = [wal - 1 for wal in args.whichatomlabel]



with open( args.output, 'w') as f:
    # default arguments
    f.write( """% pdflatex file.tex
% asy file-1.asy
% pdflatex file.tex
\n""")
    f.write( r"""
\documentclass{standalone}
\usepackage{asymptote}
\begin{document}
\begin{asy}""")

    f.write( f"""
import three;
import palette;
import graph;
settings.render={args.render};
settings.prc=false;
size(10cm);
""")
    f.write( f'currentprojection = {args.proj_type}(({args.proj[0]},{args.proj[1]},{args.proj[2]}));\n')
    f.write( f'material Hcol = material(diffusepen=gray(0.8), emissivepen=gray(0.1),specularpen=white);\n')
    f.write( f'material Ocol = material(diffusepen=red, emissivepen=gray(0.1),specularpen=white);\n')
    f.write( f'material Ncol = material(diffusepen=blue, emissivepen=gray(0.1),specularpen=white);\n')
    f.write( f'material Ccol = material(diffusepen=gray(0.1), emissivepen=gray(0.1),specularpen=white);\n')
    f.write( f'material Scol = material(diffusepen=yellow, emissivepen=gray(0.1),specularpen=white);\n')
    f.write( f'material cylcolor = material(diffusepen=gray(0.8), emissivepen=gray(0.3),specularpen=mediumgray);\n')
    #f.write( f'material cylcolor = material(emissivepen=gray(0.6));\n')
    
    f.write( f'\nreal cylRadius={args.cyl_s};\n')
    f.write( f'real sphereRadius={args.sph_s};\n\n')

    f.write( """
triple cameradirection(triple pt, projection
    P=currentprojection) {
        if (P.infinity) {
        return unit(P.camera);
    } else {
        return unit(P.camera - pt);
    }
    }
    """)
    f.write( """
triple towardcamera( triple pt, real distance=1, projection
    P=currentprojection) {
    return pt + distance * cameradirection(pt, P);
}
    """)

    f.write( f'material cylcolor0 = material(diffusepen=gray(0.5),specularpen=gray(0.30),emissivepen=gray(0.40));\n')

    if args.numlab:
        f.write( f"""
void drawRod(triple a, triple b, string lab) {{
  surface rod = extrude(scale(cylRadius)*unitcircle, axis=length(b-a)*Z);
  triple orthovector = cross(Z, b-a);
  pair parrvec = ( 1.0/orthovector.x, 1.0/orthovector.y );
  if (length(orthovector) > .01) {{
    real angle1 = aCos(dot(Z, b-a) / length(b-a));
    rod = rotate(angle1, orthovector) * rod;
  }}
  real angle1 = aCos( dot( a,b) / length(a) / length(b));
  // must project onto X-Y plane and compute angle //
  //draw(a -- b, L=rotate(angle1)*Label(lab, align=NoAlign, fontsize(2pt), position=MidPoint));
  //draw(a -- b, L=Label(lab, align=N, fontsize(2pt), position=Relative(0.5)));
  draw(shift(a)*rod, surfacepen=cylcolor0);
  triple mid = (a + b) / 2.0 + (0,0,1);
  //label( rotate( angle1, z=(0,0))*Label( lab, align=N, fontsize(2pt)), position=towardcamera( mid));
  label( Label( lab, align=N, fontsize({args.labelsize}pt)), position=towardcamera( mid));
}}
        """)
    elif args.colorlab:
        # first scale
        if args.num_bins % 2 != 0:
            args.num_bins += 1
        for i, v in enumerate( np.linspace( 0.2, 1.0, args.num_bins), start=0):
            #f.write( f'material cylcolor{args.num_bins-i}= material(specularpen=gray({v}),emissivepen=gray({v}));\n')
            # blue scale
            #f.write( f'material cylcolor{i}= material(specularpen=rgb(1-{v},1-{v},1),emissivepen=rgb(1-{v},1-{v},1));\n')
            # gray scale
            f.write( f'material cylcolor{i}= material(emissivepen=rgb(1-{v}*0.7,1-{v},1-{v}));\n')
            #f.write( f'filldraw(box((0,{i}), (2,{i+1})), red);;\n')
        # second scale
        #for i, v in enumerate( np.linspace( 0.0, 1.0, args.num_bins//2), start=args.num_bins//2):
        #    #f.write( f'material cylcolor{args.num_bins-i}= material(specularpen=gray({v}),emissivepen=gray({v}));\n')
        #    # red scale
        #    f.write( f'material cylcolor{i}= material(emissivepen=rgb(0.3,{v},{v}));\n')
        #    # greenish gray scale
        #    #f.write( f'material cylcolor{i}= material(emissivepen=rgb(0.85-{v},1-{v},1-{v}));\n')
        #    #f.write( f'filldraw(box((0,{i}), (2,{i+1})),blue);;\n')
        for n in range( args.num_bins):
            f.write( f"""
void drawRod{n}(triple a, triple b, string lab) {{
  real cylRadius1 = cylRadius;
  surface rod = extrude(scale(cylRadius1)*unitcircle, axis=length(b-a)*Z);
  triple orthovector = cross(Z, b-a);
  pair parrvec = ( 1.0/orthovector.x, 1.0/orthovector.y );
  if (length(orthovector) > .01) {{
    real angle1 = aCos(dot(Z, b-a) / length(b-a));
    rod = rotate(angle1, orthovector) * rod;
  }}
  draw(shift(a)*rod, surfacepen=cylcolor{n});
  triple mid = (a + b) / 2.0 + (0,0,1);
  label( Label( lab, align=N, fontsize({args.labelsize}pt)), position=towardcamera( mid));
}}
            """)
    elif args.sizelab:
        for n in range( args.num_bins):
            f.write( f"""
void drawRod{n}(triple a, triple b, string lab) {{
  real cylRadius1 = cylRadius * {n+1} * 0.3;
  surface rod = extrude(scale(cylRadius1)*unitcircle, axis=length(b-a)*Z);
  triple orthovector = cross(Z, b-a);
  pair parrvec = ( 1.0/orthovector.x, 1.0/orthovector.y );
  if (length(orthovector) > .01) {{
    real angle1 = aCos(dot(Z, b-a) / length(b-a));
    rod = rotate(angle1, orthovector) * rod;
  }}
  draw(shift(a)*rod, surfacepen=cylcolor0);
  triple mid = (a + b) / 2.0 + (0,0,1);
  label( Label( lab, align=N, fontsize({args.labelsize}pt)), position=towardcamera( mid));
}}
            """)
    else:
        f.write( """
void drawRod(triple a, triple b, string lab) {
  surface rod = extrude(scale(cylRadius)*unitcircle, axis=length(b-a)*Z);
  triple orthovector = cross(Z, b-a);
  pair parrvec = ( 1.0/orthovector.x, 1.0/orthovector.y );
  if (length(orthovector) > .01) {
    real angle1 = aCos(dot(Z, b-a) / length(b-a));
    rod = rotate(angle1, orthovector) * rod;
  }
  draw(shift(a)*rod, surfacepen=cylcolor);
}
        """)
        
    f.write( """
void dH(triple center) {
     draw(shift(center)*scale3(sphereRadius)*unitsphere, surfacepen=Hcol);
}
    """)

    f.write( """
void dO(triple center) {
     draw(shift(center)*scale3(sphereRadius*1.4)*unitsphere, surfacepen=Ocol);
}
    """)

    f.write( """
void dC(triple center) {
     draw(shift(center)*scale3(sphereRadius*1.3)*unitsphere, surfacepen=Ccol);
}
    """)

    f.write( """
void dN(triple center) {
     draw(shift(center)*scale3(sphereRadius*1.35)*unitsphere, surfacepen=Ncol);
}
    """)

    f.write( """
void dS(triple center) {
     draw(shift(center)*scale3(sphereRadius*1.5)*unitsphere, surfacepen=Scol);
}
    """)

    f.write( """
void label( Label L, triple pos);
    """)

    atoms, xyz = vmd.read_vmd( args.xyz)
    # write down all the atoms
    for i, x in enumerate( xyz, start=1):
        f.write( f"""
triple x{i} = ({x[0]}, {x[1]}, {x[2]});
        """)

    #f.write( 'defaultpen(fontsize(1cm));')

    # DRAW BONDS
    if args.topo is not None:
        df.reset_index()
        f.write( '\n// add specific bonds\n')
        m, M = df['bond'].min(), df['bond'].max()
        for index, row in df.iterrows():
            a1 = xyz[int(row["id1"])]
            a2 = xyz[int(row["id2"])]
            #b = round( row["bond"]*1.e6, 1)
            if int( row["id1"]) in args.whichatomlabel:
                b = round( row["bond"]*args.scalefact, args.labeldigits)
            else:
                b = ""
            # interval id
            if M != m:
                bin_id = int( (row['bond'] - m) / (M - m) * args.num_bins - 1.e-10)
            else:
                bin_id = 4
            if bin_id > args.maxval:
                bin_id = args.maxval
            #f.write( f'label( "hello", ({a1[0]},{a1[1]}));\n')
            if args.colorlab or args.sizelab:
                f.write( f'drawRod{bin_id}(x{int(row["id1"])+1}, x{int(row["id2"])+1}, "{b}");\n')
            else:
                f.write( f'drawRod(x{int(row["id1"])+1}, x{int(row["id2"])+1}, "{b}");\n')
            #f.write( f'drawRod(x{index}, x{index});\n')

    # DRAW ATOMS
    f.write( '\n// drawing atom spheres\n')
    for i,x in enumerate( xyz):
        f.write( f'd{atoms[i]}( ({x[0]},{x[1]},{x[2]}));\n')

    #if 'num' in args.displ:
    #    # draw bond labels
    #    for index, row in df.iterrows():
    #        f.write( '\n// add bond labels\n')
    #        f.write( f'label add Bonds 0/{int(row["id1"])} 0/{int(row["id2"])}\n')

    f.write( '\\end{asy}\n\\end{document}')

