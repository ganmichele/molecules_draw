#!/usr/bin/env python
"""
This code belongs to the repository:
https://gitlab.hpc.cineca.it/ceottogroup/draw_molecular_geometries

Author: Michele Gandolfi
Date:   March 2022

Program to draw simple molecules with vmd
"""

# standard packages
import os
import argparse
from argparse import RawTextHelpFormatter

# extra required packages for efficiency
import pandas as pd


# INITIALIZE PARENT ARGUMENT PARSER
parser = argparse.ArgumentParser(
          description="""
                      Draw a single molecule with vmd.
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
global_args.add_argument( '-O', '--output', default='vmd_commands.vmd', help='name of the output file')

# GROUP STYLE ARGUMENTS
style_args = parser.add_argument_group( title='style arguments', description=None)
style_args.add_argument( '--draw', default='CPK', choices=['CPK', 'lines',], help='name of drawing style')
style_args.add_argument( '--color_scheme', default='Name', choices=['Name', 'Atom',], help='name of atom coloring scheme')
style_args.add_argument( '--sph_s', default='1.0', help='atom sphere scale')
style_args.add_argument( '--bond_s', default='0.3', help='bond cylinder radius scale')

topo_args = parser.add_argument_group( title='topology arguments', description=None)
topo_args.add_argument( '-T', '--topo', default=None, help="""topology space separated value file with predefined headers. As:
id1 id2 bond
1   2   [val12]
1   3   [val13]
2   3   [val23]""")
topo_args.add_argument( '--displ', nargs='*', choices=['num', 'size', 'no'], help='how to displayt the topology')

render_args = parser.add_argument_group( title='rendering arguments', description=None)
render_args.add_argument( '--render', default='no', choices=['no', 'povray', 'tachyon',], help='rendering program')

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
    if df['bond'].max() < 1.0:
        df['bond'] = ( df['bond'] - df['bond'].min() ) / ( df['bond'].max() - df['bond'].min() ) + 0.05

    unique_bonds = set( df['bond'])


with open( args.output, 'w') as f:
    # default arguments
    f.write( '# remove menu and axes\nmenu main off\naxes location off\n')
    f.write( '# orthographics proj\ndisplay proection orthographic\n')
    f.write( '# white background, no transparency\ncolor Display Background white\ndisplay depthcue off\n')
    f.write( '# Black labels\ncolor Labels Bonds black\n')
    f.write( f'\n# load molecule\nmol new {args.xyz}\n')
    f.write( '\n# set default color and repr style\n')
    f.write( f'mol color {args.color_scheme}\nmol addrep 0\nmol modselect 1 0 all\n')
    f.write( f'mol modstyle 1 0 {args.draw} {args.sph_s} {args.bond_s}  30.0  30.0\n')
    
    if 'size' in args.displ:
        f.write( '\n# create custom selections\n')
        # assuming representations with ids 0 and 1 are already defined, start from id 2
        for n, ub in enumerate( unique_bonds, start=2):
            atom_indexes = set( df[ df['bond'] == ub][['id1','id2'] ].values.flatten() )
            atom_indexes = ' '.join( str(a) for a in atom_indexes)
            f.write( f'# selection for bonds size {ub}\n')
            f.write( 'mol addrep 0\n')
            f.write( f'mol modselect {n} 0 index 0 {atom_indexes}\n')
            f.write( f'mol modstyle {n} 0 {args.draw} {args.sph_s} {ub} 30.0 30.0\n')
    
    df.reset_index()
    f.write( '\n# add specific bonds\n')
    for index, row in df.iterrows():
        f.write( f'topo addbond {int(row["id1"])} {int(row["id2"])}\n')

    if 'num' in args.displ:
        # draw bond labels
        for index, row in df.iterrows():
            f.write( '\n# add bond labels\n')
            f.write( f'label add Bonds 0/{int(row["id1"])} 0/{int(row["id2"])}\n')

    if args.render != 'no':
        #f.write( f'\nrender {render_name} {args.render}')
        f.write( f'\npovray +W1016 +H434 -Ivmdscene.pov -Ovmdscene.pov.tga +D +X +A +FT\n')
