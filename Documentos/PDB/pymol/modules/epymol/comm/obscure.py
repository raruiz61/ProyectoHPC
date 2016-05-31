#
# Obscure an object to protect property rights.
#
import pymol
from pymol import cmd
from chempy import cpv

def perturb(v):
    """perturbs a vector"""
    return cpv.add(v, cpv.scale(cpv.random_vector(),0.4))

def obscure(objSel, hiding="medium", keep=None):
    """
    Given an object or selection, usually a small molecule, obscure it
    to protect its exact identity.

    PARAMS

    objSel
        Any object or selection to hide

    hiding
        (string; default="medium") level to which PyMOL obscures the object

    keep
       (boolean) by default, PyMOL removes the obscured object from your file, this flag
       will keep that object in the file.  Be careful!

   RETURNS
       None

    NOTES
       Use a small molecule, large ones are very, very slow.
    """
    # if the user puts erroneous hiding settings, make it medium
    if hiding not in ["low", "medium", "high"] : hiding = medium

    # these parameters are fine-tuned for a good visual experience
    hiding_levels = { "low" : { "gauss" : 2.00, "grid" : 0.18, "level" : 2.0 },
                   "medium" : { "gauss" : 3.75, "grid" : 0.25, "level" : 2.5 },
                   "high"   : { "gauss" : 4.00, "grid" : 0.33, "level" : 2.0 } }
    ga=hiding_levels[hiding]["gauss"]
    gr=hiding_levels[hiding]["grid"]
    le=hiding_levels[hiding]["level"]

    # detect if we're hiding a subset of a molecule, then add one bond, so ensure that what sticks out looks good
    t = "_target"
    cnt_sel = cmd.count_atoms(objSel)
    cnt_obj = cmd.count_atoms( "bo. " + objSel )
    if cnt_sel != cnt_obj:
        cmd.select(t,  objSel + " extend 1")
    else: 
        cmd.select(t, objSel)

    # get a new name for the map/surf
    newMap = cmd.get_unused_name("obsc_map")
    newSurf = cmd.get_unused_name("obsc_surf")

    # randomize the coordinates; should this call get_state?
    if keep==None:
        cmd.alter_state(1, t, "(x,y,z)=%s.perturb([x,y,z])" % __name__ )
        # clear out charge and b-factor
        cmd.alter(t, "b=10.0")
        cmd.alter(t, "q=1.0")

    # make the gaussian map and draw its surface
    cmd.set("gaussian_resolution", ga)
    cmd.map_new(newMap, "gaussian", gr, t)
    cmd.isosurface(newSurf, newMap, le)
    
    if keep==None:
        cmd.remove(objSel)

cmd.extend("obscure", obscure)
