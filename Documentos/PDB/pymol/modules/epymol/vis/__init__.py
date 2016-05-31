"""
PyMOL Support for vis files.

Requires:
  - h5py

"""

import epymol
from pymol import cmd, cgo, commanding, CmdException

@epymol.extendapi
def load_vis_surface(filename, name=None, mimic=1, quiet=1, _self=cmd):
    '''
DESCRIPTION

    Load a surface from a vis file.

ARGUMENTS

    filename = string: filename

    name = string: name of surface or group prefix, if ending with
    a period. {default: "Dataset Name" property from file}

    mimic = 0/1: apply styles from vis file {default: 1}
    '''
    mimic, quiet = int(mimic), int(quiet)

    data = get_data(filename)

    try:
        visibility = data.attrs.get("Visibility", 1)
        style = data.attrs.get("Style", 0)
        transparency = data.attrs.get("Transparency", 0) / 100.0
        color_enum = data.attrs.get("Color Enum", 2)
        ds_name = data.attrs.get("Dataset Name") or "unnamed"

        vertices = data[u'Coordinates of Vertices'].value.reshape((-1,3))
        normals = data[u'Normals of Vertices'].value.reshape((-1,3))
        patches = data[u'Patches'].value.reshape((-1,3))
    except:
        raise CmdException('could not find surface in vis file')

    obj = [
        cgo.BEGIN,
    ]

    if style == 1:
        obj.append(cgo.LINES)
        for line in data["Lines"].value.reshape((-1,2)):
            for vertex in vertices[line].tolist():
                obj.append(cgo.VERTEX)
                obj.extend(vertex)
    else:
        obj.append(cgo.TRIANGLES)
        for patch in patches:
            for vertex, normal in zip(vertices[patch].tolist(),
                                      normals[patch].tolist()):
                obj.append(cgo.NORMAL)
                obj.extend(normal)
                obj.append(cgo.VERTEX)
                obj.extend(vertex)

    obj.append(cgo.END)

    if not name:
        name = _self.get_unused_name(ds_name, 0)
    elif name.endswith("."):
        name = _self.get_unused_name(name + ds_name, 0)

    _self.load_cgo(obj, name)

    if not mimic:
        return name

    if transparency > 0:
        _self.set("cgo_transparency", transparency, name)

    if not visibility:
        _self.disable(name)

    apply_color(color_enum, name, quiet, _self)

    return name

def load_vis_map(filename, name=None, mimic=1, quiet=1, _self=cmd):
    '''
DESCRIPTION

    Load a map (a.k.a. volume or brick) from a vis file.

ARGUMENTS

    filename = string: filename

    name = string: name of map or group prefix, if ending with
    a period. {default: "Dataset Name" property from file}

    mimic = 0/1: load isosurfaces {default: 1}
    '''
    from chempy.brick import Brick

    mimic, quiet = int(mimic), int(quiet)

    data = get_data(filename)

    try:
        if data.attrs['Rank'] != 3:
            raise ValueError('Rank %s' % data.attrs['Rank'])
        if tuple(data.attrs['Dimensions']) != data.value.shape:
            raise ValueError('Dimensions %s' % data.attrs['Dimensions'])
        if data.attrs['Grid Type'] != 1:
            raise ValueError('Unknown Grid Type: %s' % data.attrs['Grid Type'])

        ds_name = data.attrs.get("Dataset Name") or "unnamed"

        obj = Brick.from_numpy(data.value,
                data.attrs["Space Ratio"],
                data.attrs["Origin"])

    except Exception as e:
        raise CmdException('could not find map in vis file')

    if not name:
        name = _self.get_unused_name(ds_name, 0)
    elif name.endswith("."):
        name = _self.get_unused_name(name + ds_name, 0)

    _self.load_brick(obj, name)

    try:
        cell = data.attrs['Unit Cell ABC'].tolist() + \
               data.attrs['Unit Cell Angle'].tolist()
    except KeyError:
        cell = [False]

    if all(cell):
        cell.append('P 1')
        _self.set_symmetry(name, *cell)

    if not mimic:
        return name

    isofunc = {1: _self.isomesh, 2: _self.isodot}.get(
            0, _self.isosurface)

    for color, isolevel in zip(data.attrs.get('Contour Colors', []),
                               data.attrs.get('Contour Isovalues', [])):
        mesh_name = _self.get_unused_name(name + '_contour', 0)
        isofunc(mesh_name, name, isolevel)
        apply_color(color, mesh_name, quiet, _self)

    return name

def load_vis(filename, *args, **kwargs):
    '''
DESCRIPTION

    Load a surface or map from a vis file.

SEE ALSO

    load, load_vis_surface, load_vis_map
    '''
    data = get_data(filename)

    func = load_vis_map if hasattr(data, 'value') else \
           load_vis_surface

    return func(data, *args, **kwargs)

def get_data(handle):
    if not isinstance(handle, basestring):
        return handle

    try:
        import h5py
    except ImportError:
        raise CmdException("h5py module not available, cannot load vis files")

    try:
        handle = h5py.File(handle, "r")
    except IOError:
        raise CmdException("Failed to open file " + str(handle))

    return handle.values()[0].values()[0]

def apply_color(color_enum, name, quiet=1, _self=cmd):
    try:
        from epymol.mae import dotmae

        if not dotmae._color_dict:
            dotmae.prime_colors()

        rgb = dotmae._color_dict[color_enum]
        _self.color("0x%06x" % rgb, name)
    except Exception as e:
        if not quiet:
            print "Warning: applying color %s failed" % (color_enum)
