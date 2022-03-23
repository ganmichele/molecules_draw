#!/usr/bin/env python
"""
This code belongs to the repository:
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
import pandas as pd

# Michele's custom packages
from PES import vmd


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
style_args.add_argument( '--color_scheme', default=None, help='filename for custom color scheme')
style_args.add_argument( '--sph_s', default='0.30', help='atom sphere scale')
style_args.add_argument( '--cyl_s', default='0.04', help='bond cylinder radius scale')
style_args.add_argument( '--proj_type', default='orthographic', help='type of projection')
style_args.add_argument( '--proj', nargs=3, default=[0,0,0], help='projection')

topo_args = parser.add_argument_group( title='topology arguments', description=None)
topo_args.add_argument( '-T', '--topo', default=None, help="""topology space separated value file with predefined headers. As:
id1 id2 bond
1   2   [val12]
1   3   [val13]
2   3   [val23]""")
topo_args.add_argument( '--numlab', action='store_true', help='display label on bond')
topo_args.add_argument( '--colorlab', action='store_true', help='display label on bond')
topo_args.add_argument( '--sizelab', action='store_true', help='display label on bond')

render_args = parser.add_argument_group( title='rendering arguments', description=None)
render_args.add_argument( '--render', default='auto', help='rendering value (higher for more quality). If ftype=pdf, then render should be 0, else 5 should be enough')
render_args.add_argument( '--ftype' , default='png', choices=['png',], help='format of embedded figure')

###############################################################################

# PARSE ARGUMENTS FROM STDIN
args = parser.parse_args()


if args.topo is not None:
    assert os.path.isfile( args.topo), f'No such file {args.topo}. Aborting'
    df = pd.read_csv( args.topo, comment='#', sep='\s+')
    assert df.shape[1] in {2,3}, 'Something wrong with topology file. It should have 2 or 3 named columns'
    if df.shape[1] == 2:
        df['bond'] = args.bond_s
    # ensure correct column names
    df.set_axis( ['id1' ,'id2', 'bond'], axis='columns', inplace=True)
    df.astype( {'id1': 'str', 'id2':'str', 'bond':'float64'})
    #if df['bond'].max() < 1.0:
    #    df['bond'] = ( df['bond'] - df['bond'].min() ) / ( df['bond'].max() - df['bond'].min() ) + 0.05

    unique_bonds = set( df['bond'])


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

    f.write( """
import three;
settings.render=10;
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

    if args.numlab:
        f.write( """
void drawRod(triple a, triple b, string lab) {
  surface rod = extrude(scale(cylRadius)*unitcircle, axis=length(b-a)*Z);
  triple orthovector = cross(Z, b-a);
  pair parrvec = ( 1.0/orthovector.x, 1.0/orthovector.y );
  if (length(orthovector) > .01) {
    real angle1 = aCos(dot(Z, b-a) / length(b-a));
    rod = rotate(angle1, orthovector) * rod;
  }
  real angle1 = aCos( dot( a,b) / length(a) / length(b));
  // must project onto X-Y plane and compute angle //
  //draw(a -- b, L=rotate(angle1)*Label(lab, align=NoAlign, fontsize(2pt), position=MidPoint));
  //draw(a -- b, L=Label(lab, align=N, fontsize(2pt), position=Relative(0.5)));
  draw(shift(a)*rod, surfacepen=cylcolor);
  triple mid = (a + b) / 2.0 + (0,0,1);
  //label( rotate( angle1, z=(0,0))*Label( lab, align=N, fontsize(2pt)), position=towardcamera( mid));
  label( Label( lab, align=N, fontsize(2pt)), position=towardcamera( mid));
}
        """)
    elif args.colorlab:
        f.write( f'material cylcolor14= material(specularpen=gray(0.10),emissivepen=gray(0.10));\n')
        f.write( f'material cylcolor13= material(specularpen=gray(0.15),emissivepen=gray(0.15));\n')
        f.write( f'material cylcolor12= material(specularpen=gray(0.20),emissivepen=gray(0.20));\n')
        f.write( f'material cylcolor11= material(specularpen=gray(0.25),emissivepen=gray(0.25));\n')
        f.write( f'material cylcolor10= material(specularpen=gray(0.30),emissivepen=gray(0.30));\n')
        f.write( f'material cylcolor9 = material(specularpen=gray(0.35),emissivepen=gray(0.35));\n')
        f.write( f'material cylcolor8 = material(specularpen=gray(0.40),emissivepen=gray(0.40));\n')
        f.write( f'material cylcolor7 = material(specularpen=gray(0.45),emissivepen=gray(0.45));\n')
        f.write( f'material cylcolor6 = material(specularpen=gray(0.50),emissivepen=gray(0.50));\n')
        f.write( f'material cylcolor5 = material(specularpen=gray(0.55),emissivepen=gray(0.55));\n')
        f.write( f'material cylcolor4 = material(specularpen=gray(0.60),emissivepen=gray(0.60));\n')
        f.write( f'material cylcolor3 = material(specularpen=gray(0.65),emissivepen=gray(0.65));\n')
        f.write( f'material cylcolor2 = material(specularpen=gray(0.70),emissivepen=gray(0.70));\n')
        f.write( f'material cylcolor1 = material(specularpen=gray(0.75),emissivepen=gray(0.75));\n')
        f.write( f'material cylcolor0 = material(specularpen=gray(0.80),emissivepen=gray(0.80));\n')
        for n in range( 15):
            f.write( f"""
void drawRod{n}(triple a, triple b, string lab) {{
  //real cylRadius1 = cylRadius * {n+1} * 0.4;
  real cylRadius1 = cylRadius;
  surface rod = extrude(scale(cylRadius1)*unitcircle, axis=length(b-a)*Z);
  triple orthovector = cross(Z, b-a);
  pair parrvec = ( 1.0/orthovector.x, 1.0/orthovector.y );
  if (length(orthovector) > .01) {{
    real angle1 = aCos(dot(Z, b-a) / length(b-a));
    rod = rotate(angle1, orthovector) * rod;
  }}
  draw(shift(a)*rod, surfacepen=cylcolor{n});
}}
            """)
    elif args.sizelab:
        f.write( f'material cylcolor0 = material(diffusepen=gray(0.5),specularpen=gray(0.30),emissivepen=gray(0.40));\n')
        for n in range( 15):
            f.write( f"""
void drawRod{n}(triple a, triple b, string lab) {{
  real cylRadius1 = cylRadius * {n+1} * 0.4;
  surface rod = extrude(scale(cylRadius1)*unitcircle, axis=length(b-a)*Z);
  triple orthovector = cross(Z, b-a);
  pair parrvec = ( 1.0/orthovector.x, 1.0/orthovector.y );
  if (length(orthovector) > .01) {{
    real angle1 = aCos(dot(Z, b-a) / length(b-a));
    rod = rotate(angle1, orthovector) * rod;
  }}
  draw(shift(a)*rod, surfacepen=cylcolor0);
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
    for i, x in enumerate( xyz):
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
            b = round( row["bond"]*1.e6, 1)
            # interval id
            int_id = int( (row['bond'] - m) / (M - m) * 20 - 1.e-10)
            if int_id > 8:
                int_id = 8
            #f.write( f'label( "hello", ({a1[0]},{a1[1]}));\n')
            if args.colorlab or args.sizelab:
                f.write( f'drawRod{int_id}(({a1[0]}, {a1[1]}, {a1[2]}), ({a2[0]}, {a2[1]}, {a2[2]}), "{b}");\n')
            else:
                f.write( f'drawRod(({a1[0]}, {a1[1]}, {a1[2]}), ({a2[0]}, {a2[1]}, {a2[2]}), "{b}");\n')
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

