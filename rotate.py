#!/usr/bin/env python
"""
Tools to rotate and plot small molecules
Michele Gandolfi 
June 2019
"""
import argparse
import numpy as np
from math import sin, cos, acos
import re
import os
import matplotlib.pyplot as plt
import itertools as it
import pandas as pd


parser = argparse.ArgumentParser(
        description="""
                    Helps on draw single molecules and prepare input for a tikz picture.
                    The initial configuration is read from an input xyz file.
                    The final configuration is optionally saved in an equivalent xyz file.
                    A .tex file is optionally generated with the rotated geometry
                    """,
        formatter_class= argparse.ArgumentDefaultsHelpFormatter,
        fromfile_prefix_chars='+',
        epilog='Further arguments caan be given in a file using +FILENAME')

parser.add_argument( 'input'                                  , help='an xyz geometry file for an input geometry')
parser.add_argument( '-R', '--rotation'  , nargs=3            , help='degrees of rotation around the normal axes (e.g. 60 0 180)', default=None)
parser.add_argument( '-o', '--output'                         , help='output filename')
parser.add_argument( '-c', '--center'    , nargs='?'          , help='a label of the atom that will be at the origin of the coordinate system')
parser.add_argument(       '--order'     , default='zyx'      , help='defines the order according to which the rotations are performed')
parser.add_argument(       '--plot'      , action='store_true', help='plot the resulting picture on the xy surface')
parser.add_argument(       '--tikz'                           , help='defines the output of a tikz ppicture (.tex file)')
parser.add_argument(       '--depth'     , action='store_true', help='add depths cues to bonds (not always accurate)')
parser.add_argument(       '--angles', nargs='+',default=None , help="list of angles to represent on the figure (e.g. [('atom1', 'vertex_atom', 'atom2')])")
parser.add_argument(       '--asize'     , default=0.4        , help="size of the angle arcs")
parser.add_argument(       '--connect'                        , help="connectivity csv file if the label (third) column is not provided, the labels will be bond distances")
parser.add_argument( '-A', '--atomlabel' , action='store_true', help='display atom symbol in label')
parser.add_argument(       '--numat'     , action='store_true', help='display atom number in label')
parser.add_argument( '-B', '--bondlabel' , action='store_true', help='display bond label')
parser.add_argument(       '--digits'    , default=3          , help='number of digits in bondlength and anglesize')
parser.add_argument(       '--opaque'    , default=None       , help='draw additional darker layers on atoms based on depth')

args = parser.parse_args()


def read_geo( fpath):
    line_regex = '\s*\w+\d*(\s+\-?\d+\.\d+(e|E\-?\d+)?){3}'
    at_regex = '\w+\d*'
    cor_regex = r'(-?\d+\.\d+){1}(e|E-\d+)?'
    with open( fpath, 'r') as f:
        atoms, coor = [], []
        for l in f:
            if re.match( line_regex, l):
                atom = re.search( at_regex, l).group()
                x, y, z = re.findall( cor_regex, l)
                x, y, z = clean(x), clean(y), clean(z)
                atoms.append( atom)
                coor.append( [x, y, z] )
    print( 'input read from \t{0}'.format( fpath))
    coor = np.array( coor).astype( float)
    atom_index = [ a + str(i) for i, a in enumerate( atoms)]
    geom = pd.DataFrame( coor, columns=['x', 'y', 'z'], index=atom_index)
    return atoms, coor, geom


def clean( x):
    a, b = x[0], x[1]
    if b == '':
        return a
    elif 'e' in b.lower():  # should be obsolete
        b = b.lower()
        b = b.replace('e', '')
        a = round( float(a) * 10**float(b), 8)
        return a


def generate_mat( a, axis='z'):
    m = np.array( [[  cos(a), sin(a),    0],
                   [ -sin(a), cos(a),    0],
                   [  0     , 0     ,    1]]
                )
    a_to_print = a / np.pi * 180
    if axis=='x':
        print( 'rotate around x-axis {0} degrees'.format(a_to_print))
        perm = [2, 0, 1]
    elif axis=='y':
        print( 'rotate around y-axis {0} degrees'.format(a_to_print))
        perm = [1, 2, 0]
    else: # axis=='z'
        print( 'rotate around z-axis {0} degrees'.format(a_to_print))
        perm = [0, 1, 2]
    m = m[:, perm]
    m = m[perm, :]
    return m


def rot_geo( mat, xdeg=0.0, ydeg=0.0, zdeg=0.0, order='zyx', t_mat=None, center=None):
    d2r = lambda x: x / 360 * 2 * np.pi
    xrad, yrad, zrad = d2r( xdeg), d2r( ydeg), d2r( zdeg)
    if center is not None:
        print( 'centering molecule to vector center')
        if type( center) is int:
            center = mat[ center, :]
        num_atoms = len( mat)
        center = np.array( center.tolist() * num_atoms)
        center.shape = (num_atoms, 3)
        mat -= center
    if t_mat:
        return mat @ t_mat
    else:
        order_dic = {'x':0, 'y':1, 'z':2}
        order2 = [order_dic[l] for l in order]
        x_mat = generate_mat( xrad, 'x')
        y_mat = generate_mat( yrad, 'y')
        z_mat = generate_mat( zrad, 'z')
        rots = np.array( [x_mat, y_mat, z_mat])
        rots = rots[ order2]
        t_mat = rots[0] @ rots[1] @ rots[2]
    return mat @ t_mat


def write_out( fpath, coor):
    coor = np.round( coor, 8)
    with open( fpath, 'w') as f:
        f.write( str( len(coor)) + '\n')
        f.write( '0.0\n')
        for i, r in enumerate( coor):
            line = atoms[i] + ' '
            line += ' '.join( str(num) for num in r)
            f.write( line + '\n')
    print( 'output written to \t{0}'.format( fpath))


def write2template( fpath, tpath, connect=None, ang=None, rot_note=None, colors_dic=None, bond_depth=False, asize=0.4):
    # get geometries and atoms
    atoms, coor, geom = read_geo( fpath)

    # define atoms characteristics. Add as necessary
    #colors_dic = {'H': 'white', 'C': 'gray', 'D':'pink', 'N':'blue', 'O':'red', 'S': 'yellow'}
    colors_dic = {'H': 'white', 'C': 'green!80!black', 'D':'pink', 'N':'blue', 'O':'red', 'S': 'yellow'} #FIXME
    scale_dic = {'H': 1, 'C': 1.8, 'D': 1.1, 'N': 2.0, 'O':2.2, 'S':2.6}
    #scale_dic = {'H': 0.8, 'C': 1.6, 'D': 1.1, 'N': 2.0, 'O':1.8, 'S':2.2} #FIXME

    # scale sizes according to depth (on z axiz)
    scale3 = coor[:,-1]
    mi, ma = 0.7, 1.3
    scale3 = (scale3 - scale3.min()) / (scale3.max() - scale3.min()) * (ma - mi) + mi

    # write to file
    with open( tpath, 'w') as t:
        #for i, l in enumerate( t):
        t.write("""\\documentclass[crop,tikz]{standalone}\n
\\usepackage{xcolor,amsmath,mathtools,amsfonts,amssymb,graphicx,pdfpages}
\\usetikzlibrary{automata,positioning,calc}\n
\\begin{document}\n
\\newcommand{\\tikzAngleOfLine}{\\tikz@AngleOfLine}
  \\def\\tikz@AngleOfLine(#1)(#2)#3{%
  \\pgfmathanglebetweenpoints{%
    \\pgfpointanchor{#1}{center}}{%
    \\pgfpointanchor{#2}{center}}
  \\pgfmathsetmacro{#3}{\\pgfmathresult}%
  }\n
\\begin{tikzpicture}\n""")
        t.write('[' + '\n')

        # define atom styles
        for a in set( atoms):
            t.write( '{0}/.style={{circle, fill={1}, minimum size=3*{2}mm, inner sep=1, draw=black}},\n'.format(a, colors_dic[a], scale_dic[a]))
        t.write( ']' + '\n'*2)

        # takes note of rotation used to get this
        t.write( '% rotations: {0}\n\n'.format( rot_note))

        # draw atom coordinates and scales
        for i, a in enumerate( atoms):
            c = coor[i]
            t.write( '\\node[{0}, scale={4}] ({1}) at ( {2}, {3} ) {{  }};\n'.format(a, a + str(i), c[0], c[1], scale3[i]))
        t.write('\n')

        # draw bonds and annotate lengths
        if connect:
            for k in connect:
                try:
                    connect[k] = round( connect[k], int( args.digits))
                except Exception as e:
                    pass
                if args.bondlabel:
                    t.write( '\\draw[very thick, draw=black!70!white] ({0}) -- node[sloped, anchor=center, above, scale=0.4] {{ {1} }} ({2});\n'.format(k[0], connect[k], k[1])) 
                else:
                    t.write( '\\draw[very thick, draw=black!70!white] ({0}) -- ({1});\n'.format(k[0], k[1])) 
        t.write('\n')

        if bond_depth:
            zz = coor[:,-1]
            cue_connect = { k: zz[ int(k[0][-1])] - zz[ int(k[1][-1])] for k in connect}
            sizescale = 0.12
            cue_points = list( cue_connect.keys())
            cues = []
            for p in cue_points:
                a, b = p
                avec = geom.loc[ a, ['x','y']].values
                bvec = geom.loc[ b, ['x','y']].values
                size = abs( geom.loc[ a, 'z'] - geom.loc[ b, 'z'] ) * sizescale
                u, d = cue_size( bvec, avec, size=size)
                if geom.loc[a,'z'] > geom.loc[b,'z']:
                    r = (a, b, u, d)
                else:
                    #u, d = cue_size( bvec, avec, size)
                    r = (b, a, d, u)
                cues.append( r)
                t.write( '\\fill[fill=black!70!white] ($ ({0}.center) + {2} $) -- ($ ({0}.center) + {3} $) -- ({1}.center);\n'.format(r[0], r[1], r[2], r[3])) 
        
        # draw and annotate angles
        if ang:
            for k in ang:
                t.write("""
\\tikzAngleOfLine({1})({0}){{\\AngleStart}}
\\tikzAngleOfLine({1})({2}){{\\AngleEnd}}
\\draw[black,<->,very thin] ({1})+(\\AngleStart:{4}cm) arc (\\AngleStart:\\AngleEnd:{4}cm);
\\node[circle, scale=0.4] at ($({1})+({{(\\AngleStart+\\AngleEnd)/2}}:{5}cm)$) {{${3}$}};
""".format(k[0], k[1], k[2], round( ang[k], int(args.digits)), asize, float(asize)+0.1))
            t.write('\n')
        
        # draw atoms again, to be on top of all
        order = np.argsort( geom.loc[:, 'z'])
        atoms = np.array( atoms)
        atoms = atoms[ order]
        coor = coor[order, :]
        geom = geom.iloc[order, :]
        scale3 = scale3[ order]
        regex_chem = '^[A-Za-z]+'
        regex_num = '[0-9]+$'
        for i, a in enumerate( geom.index):
            c = tuple( geom.loc[a, ['x', 'y']].values)
            chem = re.match( regex_chem, a).group()
            num_chem = re.search( regex_num, a).group()
            opaque = float( args.opaque) - ( i/len(geom.index) * float( args.opaque) )
            if args.numat and args.atomlabel:
                t.write( '\\node[{0}, scale={3}] ({1}) at {2} {{ {4} }};\n'.format(chem, a, c, scale3[i], '$\mathrm{{ {0} }}_{{ {1} }}$'.format( chem, num_chem) )) 
            elif args.atomlabel:
                t.write( '\\node[{0}, scale={3}] ({1}) at {2} {{ {4} }};\n'.format(chem, a, c, scale3[i], '$\mathrm{{ {0} }}$'.format( chem) )) 
            else:
                t.write( '\\node[{0}, scale={3}] ({1}) at {2} {{ }};\n'.format(chem, a, c, scale3[i] )) 

            if args.opaque is not None: # draw additional layer of opaque atoms
                t.write( '\\node[{0}, scale={3}, fill=black!50!white, opacity={4}] ({1}) at {2} {{ }};\n'.format(chem, a, c, scale3[i], opaque )) 

        t.write( '\n\\end{tikzpicture}\n')
        t.write( '\n\\end{document}\n')


def cue_size( beg, end, size=0.1):
    v = end - beg
    if np.allclose( v, 0, atol=1e-05):
        vh = np.array( (0, 0))
    elif np.round( v[0], 5)==0:
        vh = np.array( (1, 0))
    elif np.round( v[1], 5)==0:
        vh = np.array( (0, 1))
    else:
        mh = - v[0] / v[1]
        vh = np.array( (1, mh))
        vh = vh / np.sqrt(vh @ vh) * 0.5 * size
    #endup   = end + vh
    #enddown = end - vh
    endup   = + vh
    enddown = - vh
    return tuple( endup), tuple( enddown)


def distance_matrix( coor, metric='euclidean', conversion=None):
    if conversion is not None:
        print( 'rescaling distances with {0} factor'.format( conversion))
    if metric=='euclidean':
        d = lambda x, y: np.sqrt( (x - y) @ (x - y)) * conversion
    natoms = len( coor)
    dist = np.zeros( (natoms, natoms))
    for i in range( natoms):
        for j in range( i+1, natoms):
            dist[i,j] = d( coor[i], coor[j])
            dist[j,i] = dist[i,j]
    return dist


def calc_angles( atoms, coor, selection=None):
    # should be deprecated
    car = coor[0]
    a_dic = {}
    for c1i, c1 in enumerate( coor):
        for c2i, c2 in enumerate( coor):
            if c1i==0 or c2i==0 or c1i==c2i or c2i < c1i:
                continue
            norm_c1 = np.sqrt( c1 @ c1)
            norm_c2 = np.sqrt( c2 @ c2)
            a = np.arccos( c1 @ c2 / norm_c1 / norm_c2)
            a = a / np.pi * 180
            atom1 = atoms[c1i] + str(c1i)
            atom2 = atoms[c2i] + str(c2i)
            a_dic[ atom1, atom2] = a
    if selection not in {None, 'all'}:
        a_dic = {k: a_dic[k] for k in selection}
    return a_dic


def calc_angles2( geom, selection=None):
    atoms = list( geom.index)
    coor = geom.values
    if selection == 'all':
        selection = []
        for vi, v in enumerate( atoms):
            for xi in range( 0, len(atoms)):
                if vi==xi: continue
                for yi in range( xi+1, len(atoms)):
                    if vi==yi: continue
                    selection.append( (atoms[xi], v, atoms[yi]))
    a_dic = {}
    for s in selection:
        v1 = geom.loc[s[0],:] - geom.loc[s[1],:]
        v2 = geom.loc[s[2],:] - geom.loc[s[1],:]
        norm = lambda x: np.sqrt(x @ x)
        a_dic[s] = acos( (v1 @ v2) / norm(v1) / norm(v2) ) * 180 / np.pi
    return a_dic


if __name__ == '__main__':
    atoms, coor, geom = read_geo( args.input)
    atoms_lab = [a + str(i) for i, a in enumerate( atoms)]
    degs = np.array( ['0', '0', '0']).astype( float) if not args.rotation else np.array( args.rotation).astype( float)
    c = 0 if not args.center else int( args.center)
    new_coor = rot_geo( coor, xdeg=degs[0], ydeg=degs[1], zdeg=degs[2], order=args.order, center=c)
    if not args.output and args.tikz:
        args.output = 'example.xyz'
    d = distance_matrix( coor, conversion=0.529177249)
    d = pd.DataFrame( d, index=atoms_lab, columns=atoms_lab)
    if args.output:
        write_out( args.output, new_coor)
    if args.angles:
        if args.angles[0] == 'all':
            print( 'drawing all angles')
            angles = args.angles[0]
        else:
            regex_angs = '[A-Za-z]+[0-9]+'
            angles = list( map( lambda x: tuple( re.findall( regex_angs, x)), args.angles))
        angles = calc_angles2( geom, selection=angles)
    else:
        angles=args.angles
    if args.connect:
        adj = np.genfromtxt( args.connect, delimiter=',', dtype=str)
        if adj.shape[1] == 3:
            print( 'using custom bond labels')
            adj_mat = { tuple((row[0], row[1])): row[2] for row in adj}
        else:
            print( 'using distances as bond labels')
            adj_mat = { tuple((row[0], row[1])): d.loc[row[0], row[1]] for row in adj}
    else:
        adj_mat = None
    print( '\nUSING CONNECTIVITY', adj_mat, '\n')
    if args.tikz is not None:
        print( 'writing a tikz picture in {0}'.format( args.tikz))
        write2template( args.output, args.tikz, connect=adj_mat, ang=angles, rot_note=[degs, args.order], bond_depth=args.depth, asize=args.asize)
    if args.plot:
        size = new_coor[:,2]
        size[0] *= 5
        ma, mi = max(size), min(size)
        tma, tmi = 50, 10
        scaled_size = (size - mi) / (ma - mi) * (tma - tmi) + tmi
        fig, ax = plt.subplots()
        ax.scatter( new_coor[:,0], new_coor[:,1], scaled_size)
        #n = ['C', 'H1', 'H2', 'H3', 'H4']
        for i, txt in enumerate(atoms_lab):
            ax.annotate(txt, (new_coor[i,0], new_coor[i,1]))
        plt.xlabel( 'x axis')
        plt.ylabel( 'y axis')
        plt.show()

        #################################################################
        # CODE TO PLOT 3D
        #from matplotlib import pyplot
        #from mpl_toolkits.mplot3d import Axes3D

        #fig = pyplot.figure()
        #ax = Axes3D(fig)

        #sequence_containing_x_vals = new_coor[:,0]
        #sequence_containing_y_vals = new_coor[:,1]
        #sequence_containing_z_vals = new_coor[:,2]

        #ax.scatter(sequence_containing_x_vals, sequence_containing_y_vals, sequence_containing_z_vals)
        #pyplot.xlabel( 'x-axis')
        #pyplot.ylabel( 'y-axis')
        #pyplot.zlabel( 'z-axis')
        #pyplot.show()
