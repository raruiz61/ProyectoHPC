
from pymol import cmd
import operator

def prime_spectrum_energy(objs):
    cmd.spectrum("p['r_psp_Prime_Energy']", "white blue", ' '.join(objs))
    
