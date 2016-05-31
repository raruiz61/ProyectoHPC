

from pymol import cmd
import operator

def color_by(objs):
    cmd.spectrum("p['r_mmod_Potential_Energy-OPLS-2005']", "red blue", ' '.join(objs))
