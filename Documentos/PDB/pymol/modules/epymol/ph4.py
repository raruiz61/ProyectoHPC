'''
Simple Pharmacophore CGO library for PyMOL

(c) 2014 Schrodinger, Inc.
Author: Thomas Holder <thomas.holder@schrodinger.com>

'''

import os.path
import pymol
from pymol import  cgo
from chempy import cpv

color_A = [1., .5, .5]
color_D = [.4, .7, 1.]

def _name_from_filename(filename, ext=''):
    name = os.path.basename(filename)
    if ext and name.endswith(ext):
        name = name[:-len(ext)]
    return name

def _cgo_arrow(position=(0., 0., 0), direction=(1., 0., 0.), length=1., radius=.2,
        offset=0., color='red'):
    '''
    Generate and return an arrow CGO from (position + offset * direction) to
    (position + ((length + offset) * direction)).
    '''
    if isinstance(color, str):
        color = list(pymol.cmd.get_color_tuple(color))

    hlength = radius * 2.5
    hradius = radius * 2.0

    normal = cpv.normalize(direction)
    xyz1 = cpv.add(position, cpv.scale(normal, offset))
    xyz2 = cpv.add(position, cpv.scale(normal, length + offset))
    xyz3 = cpv.add(xyz2, cpv.scale(normal, hlength))

    return [cgo.CYLINDER] + xyz1 + xyz2 + [radius] + color + color + \
           [cgo.CONE] + xyz2 + xyz3 + [hradius, 0.0] + color + color + \
           [1.0, 0.0]

def _cgo_sphere(center=(0., 0., 0.), radius=1.0, color=''):
    '''
    Generate and return a sphere CGO with given center, radius and color.
    '''
    if color and isinstance(color, str):
        color = list(pymol.cmd.get_color_tuple(color))
    obj = []
    if color:
        obj.append(cgo.COLOR)
        obj.extend(color)
    obj.append(cgo.SPHERE)
    obj.extend(center)
    obj.append(radius)
    return obj

def load_hypothesis_xyz(filename, name='', ref=''):
    '''
    Load pharmacophores from an hypothesis.xyz file.

    Directions may be derived from from Q features present in "filename"
    or "ref" (reference xyz file with all features). Without Q features,
    Rings will have no normals and acceptors/donors will be simple spheres.
    '''
    import numpy

    if not name:
        name = _name_from_filename(filename, '.xyz')

    colors = {
        'A': color_A,   # acceptor
        'D': color_D,   # donor
        'N': 'red',     # negative
        'P': 'blue',    # positive
        'H': 'forest',  # hydrophobic
        'X': 'gray',    # custom
        'Y': 'gray',    # custom
        'Z': 'gray',    # custom
        'Q': 'pink',    # target
        'R': 'orange',  # ring
    }

    sites = []
    Qs_coords = []

    # read hypotheses from file
    for line in open(ref or filename):
        a = line.split()
        coords = numpy.asfarray(map(float, a[2:5]))
        sites.append([a[1], coords, [], int(a[0])])
        if a[1] == 'Q':
            Qs_coords.append(coords)

    # match Qs to closest non-Q
    for coord in Qs_coords:
        site = min(sites,
                key=lambda site: (numpy.inf if site[0] == 'Q'
                    else ((coord - site[1])**2).sum()))
        # add direction to site
        site[2].append(coord - site[1])

    # if both filename and ref are provided, filter ref based on filename
    if ref:
        incl = set(int(line.split()[0]) for line in open(filename))
        sites = [site for site in sites if site[3] in incl]

    obj = []
    for site in sites:
        typ, center, dirs = site[:3]
        color = colors.get(typ, 'yellow')

        if typ == 'R':
            if not dirs:
                print 'error: no normal for ring'
                dirs = [(1., 0., 0.)]
            obj += cgo.torus(center, dirs[0], radius=.5, cradius=.2, color=color)
        else:
            obj += _cgo_sphere(center, color=color)
            for d in dirs:
                obj += _cgo_arrow(center, d, length=.7, offset=.9, color=color)

    pymol.cmd.load_cgo(obj, name)

def load_hypothesis_xvol(filename, name=''):
    '''
    Load excluded volumes from a hypothesis.xvol file.
    '''
    from epymol import m2io
    from chempy import Atom, models

    if not name:
        name = _name_from_filename(filename, '.xvol')

    reader = m2io.M2IOReader(filename)

    try:
        table = reader.find_block('f_m_table').m_row
        col_radius = table.keys.index('r_m_phase_radius')
        cols_xyz = map(table.keys.index, ['r_m_phase_x', 'r_m_phase_y', 'r_m_phase_z'])
    except (AttributeError, IndexError):
        print ' Error: could not find table'
        return

    model = models.Indexed()

    for row in table:
        atom = Atom()
        atom.coord = [row[i] for i in cols_xyz]
        atom.vdw = row[col_radius]
        atom.trgb = 0x00FFFF00
        atom.visible = 0x2 # spheres
        model.add_atom(atom)

    pymol.cmd.load_model(model, name)
