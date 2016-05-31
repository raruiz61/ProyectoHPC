#! /usr/bin/python

# This script will set appropriate path to Bepipred executable.

import os
import re
import subprocess

def bepipred():
    mainpath = os.path.abspath('.')
    srcpath = os.path.join(mainpath,'src')
    bepipredpath = os.path.join(srcpath,'bepipred-1.0b')
    
    p = subprocess.Popen(["which","gawk"], stdout=subprocess.PIPE)
    gawk = p.communicate()[0].strip()
        
    fname = os.path.join(bepipredpath, 'bepipred.template')
    fname_out = os.path.join(bepipredpath, 'bepipred')
 
    fin = open(fname, 'r')
    fout = open(fname_out, 'w')   
    
    content = fin.readlines()
    for line in content:
        fout.write(line.replace('/bepipred-path', srcpath).replace('/gawk-path', gawk))
        
    fin.close()
    fout.close()
    
    print '* Bepipred script has been configured.\n|-%s' % fname_out
    print '* gawk path: %s\n' % gawk
    
    os.chmod(fname_out, 0755)

if __name__ == '__main__':
    bepipred()
