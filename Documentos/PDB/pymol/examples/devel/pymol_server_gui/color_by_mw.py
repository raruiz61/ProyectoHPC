

from pymol import cmd
import operator

def color_by(objs):
    cmd.spectrum("p['MW']", "red blue", ' '.join(objs))

