#!/usr/bin/python
#
#   Usage: superimposition.py STATIC MOBILE
#   Utrecht @ 2009
#   Joao Rodrigues
#
 
import __main__
__main__.pymol_argv = [ 'pymol', '-qc'] # Quiet and no GUI
 
import sys, time, os
import pymol
 
pymol.finish_launching()
 
##
# Read User Input
staticStructurePath = os.path.abspath(sys.argv[1])
staticStructureName = staticStructurePath.split('/')[-1].split('.')[0]
mobileStructurePath = os.path.abspath(sys.argv[2])
mobileStructureName = mobileStructurePath.split('/')[-1].split('.')[0]
 
# Load Structures
 
pymol.cmd.load(mobileStructurePath, mobileStructureName)
pymol.cmd.load(staticStructurePath, staticStructureName)
 
# CEAlign STATIC, MOBILE
# CEAlign produces same alignment as the complex SUPER below
pymol.cmd.do("run /home/joao/Software/cealign-0.9/cealign.py") # Import Module
pymol.cmd.do("cealign %s, %s" %(staticStructureName, mobileStructureName))
time.sleep(1) # Dunno why, but if I don't wait, structures do not align properly..
# Save Superimposition
# save(file, selection, state (0 default), format)
pymol.cmd.save("%s_%s.pdb" %(mobileStructureName, staticStructureName), mobileStructureName, 0, 'pdb')
 
## SUPER - old
#pymol.cmd.super((staticStructureName and (resn ZN around 5 and (resn CYS or resn HIS))), (mobileStructureName and (resn ZN around 5 and (resn CYS or resn HIS))))
#pymol.cmd.save("%s_%s_SUPER.pdb" %(mobileStructureName, staticStructureName), mobileStructureName, 0, 'pdb')
 
# Get out!
pymol.cmd.quit()
