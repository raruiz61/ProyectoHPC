#! /usr/bin/env python
# Copyright (c) 2008 Robert L. Campbell (rlc1@queensu.ca)

from __future__ import print_function
import math, sys, re, gzip
import numpy as N

# import my stats.py statistical functions
import stats

# blank for spacing output (obviously)
blank=''
Mass_table={
  ' H':   1.0079,
  'HE':   4.0026,
  'LI':   6.9410,
  'BE':   9.0122,
  ' B':  10.8110,
  ' C':  12.0110,
  ' N':  14.0067,
  ' O':  15.9994,
  ' F':  18.9984,
  'NE':  20.1797,
  'NA':  22.9898,
  'MG':  24.3050,
  'AL':  26.9815,
  'SI':  28.0855,
  ' P':  30.9738,
  ' S':  32.0660,
  'CL':  35.4527,
  'AR':  39.9480,
  ' K':  39.0983,
  'CA':  40.0780,
  'SC':  44.9559,
  'TI':  47.8800,
  ' V':  50.9415,
  'CR':  51.9961,
  'MN':  54.9380,
  'FE':  55.8470,
  'CO':  58.9332,
  'NI':  58.6934,
  'CU':  63.5460,
  'ZN':  65.3900,
  'GA':  69.7230,
  'GE':  72.6100,
  'AS':  74.9216,
  'SE':  78.9600,
  'BR':  79.9040,
  'KR':  83.8000,
  'RB':  85.4678,
  'SR':  87.6200,
  ' Y':  88.9059,
  'ZR':  91.2240,
  'NB':  92.9064,
  'MO':  95.9400,
  'TC':  97.9072,
  'RU': 101.0700,
  'RH': 102.9055,
  'PD': 106.4200,
  'AG': 107.8682,
  'CD': 112.4110,
  'IN': 114.8180,
  'SN': 118.7100,
  'SB': 121.7600,
  'TE': 127.6000,
  ' I': 126.9045,
  'XE': 131.2900,
  'CS': 132.9054,
  'BA': 137.3270,
  'LA': 138.9055,
  'CE': 140.1150,
  'PR': 140.9076,
  'ND': 144.2400,
  'PM': 144.9127,
  'SM': 150.3600,
  'EU': 151.9650,
  'GD': 157.2500,
  'TB': 158.9253,
  'DY': 162.5000,
  'HO': 164.9303,
  'ER': 167.2600,
  'TM': 168.9342,
  'YB': 173.0400,
  'LU': 174.9670,
  'HF': 178.4900,
  'TA': 180.9479,
  ' W': 183.8400,
  'RE': 186.2070,
  'OS': 190.2300,
  'IR': 192.2200,
  'PT': 195.0800,
  'AU': 196.9665,
  'HG': 200.5900,
  'TL': 204.3833,
  'PB': 207.2000,
  'BI': 208.9804,
  'PO': 208.9824,
  'AT': 209.9871,
  'RN': 222.0176,
  'FR': 223.0197,
  'RA': 226.0254,
  'AC': 227.0278,
  'TH': 232.0381,
  'PA': 231.0359,
  ' U': 238.0289,
  'NP': 237.0480,
  'PU': 244.0642,
  'AM': 243.0614,
  'CM': 247.0703,
  'BK': 247.0703,
  'CF': 251.0796,
  'ES': 252.0830,
  'FM': 257.0951,
  'MD': 258.1000,
  'NO': 259.1009,
  'LR': 262.1100,
}

amino_acid_names_list = ['ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU',
                         'MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR']
solvent_names_list = ['HOH','H2O','WAT','SOL','TIP']
aa_one_letter_list = ['A','C','D','E','F','G','H','I','K','L',
                      'M','N','P','Q','R','S','T','V','W','Y']
# {{{ cross_product
def cross_product(a,b):
  """
  Cross product of two Numpy 3-d vectors.
  """
#left over from numeric (is type needed?)
#  cross = N.zeros(3,a.typecode())
  cross = N.zeros(3)
  cross[0] = a[1]*b[2]-a[2]*b[1]
  cross[1] = a[2]*b[0]-a[0]*b[2]
  cross[2] = a[0]*b[1]-a[1]*b[0]
  return cross
# }}}

# {{{ torsion
def torsion(a,b,c,d):
  """
  Calculate a torsion angle from 4 vectors (a,b,c,d)
  """
  p = a - b
  q = c - b
  r = b - c
  s = d - c
#  print("a,b,p",a,b,p)
#  print("c,b,q",c,b,q)
#  print("b,c,r",b,c,r)
#  print("d,c,s",d,c,s)
  u = cross_product(p,q)
  v = cross_product(r,s)
#  print("p,q,u",p,q,u)
#  print("r,s,v",r,s,v)
#  print("u,v ",u,v)
#  print("dot(u,v)",N.dot(u,v))
#  print("mag(u)*mag(v)",mag(u),mag(v),mag(u)*mag(v), "N.pi",N.pi)
  torsion = N.arccos(N.dot(u,v)/(mag(u)*mag(v)))*180/N.pi
#  print("torsion",torsion,"\n")
  w = cross_product(u,v)
  if N.dot(r,w) > 0:
    torsion = -torsion
  return torsion
# }}}

# {{{ mag
def mag(x):
  """
  Calculate the magnitude of vector x
  """
  return math.sqrt(N.dot(x,x))
# }}}

# {{{ dist
def dist(x1,x2):
  """
  Calculate the distance between two coordinates.
  Returns a float
  """
  return N.sqrt((x1[0]-x2[0])**2 + (x1[1]-x2[1])**2 + (x1[2]-x2[2])**2)
# }}}

# {{{ distx2
def distx2(x1,x2):
  """
  Calculate the square of the distance between two coordinates.
  Returns a float
  """
  distx2 = (x1[0]-x2[0])**2 + (x1[1]-x2[1])**2 + (x1[2]-x2[2])**2
  return distx2
# }}}

# {{{ distx2_a
def distx2_a(x1,x2):
  """
  Calculate the square of the distance between two coordinates.
  Returns a float
  This version uses the numpy method
    distx2 = N.sum(N.subtract(x1,x2)**2)
  which is equivalent to
    distx2 = N.add.reduce(N.subtract(x1,x2)**2)
  and is considerably slower than calculating directly:
    distx2 = (x1[0]-x2[0])**2 + (x1[1]-x2[1])**2 + (x1[2]-x2[2])**2
  """
# not needed
#  x1_a = N.asarray(x1)
#  x2_a = N.asarray(x2)

  distx2 = N.sum(N.subtract(x1,x2)**2)
  return distx2
# }}}

# {{{ coordiff
def coordiff(a,b):
  """
  Calculate the distance between two coordinate sets (two Protein class objects) a and b.
  Returns a list of floats
  """
  if len(a.atnum) != len(b.atnum):
    print("Sorry, unequal number of atoms -- I don't know how to deal with that!")
    return None

  diff=[]
  for i in range(len(a.atnum)):
    diff.append(dist(a.coord[i],b.coord[i]))

  return diff
# }}}

# {{{ rmsdiff
def rmsdiff(a,b):
  """
  Calculate the RMS difference between two coordinate sets (two Protein
  class objects) containing identical sequences or at least identical sets
  of atoms.
  Returns two dictionaries and a float:
    rms[chain][resnum],
    res_dict[chain][resnum],
    totalrms
  """
  if len(a.atnum) != len(b.atnum):
    print("Sorry, unequal number of atoms -- I don't know how to deal with that!")
    return None

# initialize dictionaries
  sumdist2 = {}
  natom = {}
# res_dict for keeping track of residue names when passing back the results.
  res_dict = {}
  totdist2 = 0
  totnum = 0
# create a dictionary of rmsdifferences indexed by chain and resnumber: rms[chain][resnum]
  rms = {}

  for i in range(len(a.atnum)):
    dist2 = distx2(a.coord[i],b.coord[i])
    if a.chain[i] in sumdist2:
      if a.resnum[i] in sumdist2[a.chain[i]]:
        sumdist2[a.chain[i]][a.resnum[i]] += dist2
        natom[a.chain[i]][a.resnum[i]] += 1
#        res_dict[a.chain[i]][a.resnum[i]] = a.resname[i]
#        print("Has key: a.resnum",a.resnum[i])
      else:
        sumdist2[a.chain[i]][a.resnum[i]] = dist2
        natom[a.chain[i]][a.resnum[i]] = 1
        res_dict[a.chain[i]][a.resnum[i]] = a.resname[i]
#        print("Did not have key: a.resnum",a.resnum[i])
    else:
      sumdist2[a.chain[i]] = {a.resnum[i]:dist2}
      natom[a.chain[i]] = {a.resnum[i]:1}
      res_dict[a.chain[i]] = {a.resnum[i]:a.resname[i]}
#      print("did not have key: a.chain",a.chain[i])

    totdist2 += dist2
    totnum += 1

  for ch in sumdist2.keys():
    for res in sumdist2[ch].keys():
      if natom[ch][res] != 0:
        rmsdist = N.sqrt(sumdist2[ch][res]/float(natom[ch][res]))
      if ch in rms:
        rms[ch][res] = rmsdist
      else:
        rms[ch] = {res:rmsdist}

  if totnum != 0:
    totrms = N.sqrt(totdist2/float(totnum))
  return rms,res_dict,totrms
# }}}

# {{{ rmsd
def rmsd(x,y):
  """
  Calculate the RMS difference between two coordinate sets (two Protein
  class objects) containing identical numbers of atoms.
  Returns two dictionaries and a float:
    rms[chain][resnum],
    res_dict[chain][resnum],
    totalrms
  """
  if len(x) != len(y):
    print("Sorry, unequal number of atoms -- I don't know how to deal with that!")
    return None

# initialize sum and counter
  sumdist2 = 0
  natom = 0

  for i in range(len(x)):
    sumdist2 += distx2(x[i],y[i])
    natom += 1

  rmsdist = N.sqrt(sumdist2/float(natom))

  return rmsdist
# }}}

# {{{ radius_gyration
def radius_gyration(P,mass_weighted=1):
  """
  MyPDB.radius_gyration(P), where P is a MyPDB.Protein instance

  Returns the radius of gyration of P
  R(g) = SQRT ( SUM ( m<sub>i . (x<sub>i - centre)**2 ) / SUM (m<sub>i))
  """
  centre,stdev,median,max_coord,min_coord = stats.stats(P.coord)
  #return N.sqrt(N.sum(N.multiply(N.sum(N.subtract(P.coord,centre)**2,dim=1),P.mass))/N.sum(P.mass))
  if mass_weighted:
    return N.sqrt(N.sum(N.multiply(N.sum(N.subtract(P.coord,centre)**2,axis=1),P.mass))/N.sum(P.mass))
  else:
    return N.sqrt(N.sum(N.sum(N.subtract(P.coord,centre)**2,axis=1))/P.atcounter)
# }}}

# {{{ mw
def mw(P):
  """
  MyPDB.mw(P), where P is a MyPDB.Protein instance

  Returns mw of P (molecular weight)
  """
  mw = 0
  for m in P.mass:
    mw += m

  return mw
# }}}

# {{{ coordstats
def coordstats(P):
  """
  MyPDB.coordstats(P), where P is a MyPDB.Protein instance

  Returns centre,stdev,median,max_coord,min_coord of P
  """
  centre,stdev,median,max_coord,min_coord = stats.stats(P.coord)
  return centre,stdev,median,max_coord,min_coord
# }}}

# {{{ apply_rt_matrix
def apply_rt_matrix(P,rt_matrix):
  rt_P = Protein()
  for h in P.header:
    rt_P.header.append(h)

  rt_P.header.append('REMARK rotated and translated according to')
  for i in range(len(rt_matrix)):
    rt_P.header.append('REMARK rot-trans_matrix %s' % str(rt_matrix[i]))

  for i in range(len(P.atnum)):
    rt_P.atnum.append(i)
    rt_P.atname.append(P.atname[i])
    rt_P.resname.append(P.resname[i])
    rt_P.resnum.append(P.resnum[i])
    rt_P.coord.append(N.dot(rt_matrix[0:3,0:3],P.coord[i]) + rt_matrix[:3,3])
    rt_P.velocity.append(P.velocity[i])
    rt_P.label.append(P.label[i])
    rt_P.atalt.append(P.atalt[i])
    rt_P.resext.append(P.resext[i])
    rt_P.chain.append(P.chain[i])
    rt_P.elem.append(P.elem[i])
    rt_P.b.append(P.b[i])
    rt_P.occ.append(P.occ[i])

  return rt_P
# }}}

class Protein:
# blank for spacing output (obviously)
  blank=''

# {{{ __init__
  def __init__(self):
    """
    Initialize Protein class object (creates empty lists for coordinates etc. to
    be filled upon calling the readPDB method
    """
    self.title=''
    self.atcounter = 0
    self.allatom_line = []

    self.header = []
    self.footer = []

    self.label = []
    self.atnum = []
    self.elem = []
    self.mass = []
    self.atname = []
    self.atalt = []
    self.resname = []
    self.chain = []
    self.resnum = []
    self.resext = []
    self.coord = []
    self.occ = []
    self.b = []
    self.sequence = {}
    self.box = []
    self.velocity = []

# make some "PyMOL"-like synonyms
    self.name = self.atname
    self.resn = self.resname
    self.resi = self.resnum

    self.ca_line = []
    self.ca_x = []

    self.neighbours = {}
# }}}

# {{{ concat
  def concat(self,a,b):
    for i in range(len(a.atnum)):
      self.atnum.append(i)
      self.atname.append(a.atname[i])
      self.resname.append(a.resname[i])
      self.resnum.append(a.resnum[i])
      self.coord.append(a.coord[i])
    for j in range(len(b.atnum)):
      self.atnum.append(j)
      self.atname.append(b.atname[j])
      self.resname.append(b.resname[j])
      self.resnum.append(b.resnum[j])
      self.coord.append(b.coord[j])
# }}}

# {{{ __add__
  def __add__(a,b):
    new=Protein()
    new.title = a.title  + " + " + b.title
    new.atcounter = a.atcounter + b.atcounter
    for i in range(len(a.atnum)):
      new.atnum.append(i)
      new.atname.append(a.atname[i])
      new.resname.append(a.resname[i])
      new.resnum.append(a.resnum[i])
      new.coord.append(a.coord[i])
      new.velocity.append(a.velocity[i])
      new.label.append(a.label[i])
      new.atalt.append(a.atalt[i])
      new.resext.append(a.resext[i])
      new.chain.append(a.chain[i])
      new.elem.append(a.elem[i])
      new.b.append(a.b[i])
      new.occ.append(a.occ[i])

    for j in range(len(b.atnum)):
      new.atnum.append(j)
      new.atname.append(b.atname[j])
      new.resname.append(b.resname[j])
      new.resnum.append(b.resnum[j])
      new.coord.append(b.coord[j])
      new.velocity.append(b.velocity[j])
      new.label.append(b.label[j])
      new.atalt.append(b.atalt[j])
      new.resext.append(b.resext[j])
      new.chain.append(b.chain[j])
      new.elem.append(b.elem[j])
      new.b.append(b.b[j])
      new.occ.append(b.occ[j])

    centre,stdev,median,max_coord,min_coord = coordstats(new)
    new.box = max_coord-min_coord

    return new
# }}}

# {{{ readgmx
  def readgmx(self,file_in,debug=0):
    lines = open(file_in).readlines()
    self.readgmxlines(lines,debug)
# }}}

# {{{ readgmxlines
  def readgmxlines(self,lines,debug=0):
    """
    Lines contain the following information (top to bottom):

    * title string (free format string, optional time in ps after 't=')
    * number of atoms (free format integer)
    * one line for each atom (fixed format, see below)
    * box vectors (free format, space separated reals), values: v1(x) v2(y) v3(z) v1(y) v1(z) v2(x) v2(z) v3(x) v3(y),
      the last 6 values may be omitted (they will be set to zero). Gromacs only supports boxes with v1(y)=v1(z)=v2(z)=0.

    * residue number (5 positions, integer)
    * residue name (5 characters)
    * atom name (5 characters)
    * atom number (5 positions, integer)
    * position (in nm, x y z in 3 columns, each 8 positions with 3 decimal places)
    * velocity (in nm/ps (or km/s), x y z in 3 columns, each 8 positions with 4 decimal places)

  C format
    "%5d%5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f"

    """
    self.title = lines[0][:-1]
    self.atcounter = int(lines[1][:-1])
#    for i in range(self.atcounter):
    if debug:
      print(self.title)
      print(self.atcounter)
      print(lines[2][:-1], len(lines[2]))
    for line in lines[2:self.atcounter + 3][:-1]:
      self.resnum.append(int(line[0:5]))
      self.resname.append(line[5:10].strip())
      self.atname.append(line[10:15].strip())
      self.atnum.append(int(line[15:20]))
      first_decimal = line.index('.')
      second_decimal = line[first_decimal+1:].index('.')
      incr = second_decimal + 1
      self.coord.append((float(line[20:20+incr])*10., float(line[28:20+2*incr])*10., float(line[36:20+3*incr])*10.))
      if len(line) == 68:
        self.velocity.append((float(line[44:20+4*incr])*10., float(line[52:5*incr])*10., float(line[60:6*incr])*10.))
      else:
        self.velocity.append((0.0,0.0,0.0))

# add things that might be in PDB format, but aren't stored in .gro format (in case we want
# to write out a PDB file

      self.label.append('ATOM')
      self.atalt.append('')
      self.resext.append('')
      self.chain.append('')
      self.elem.append('')
      self.b.append(0.0)
      self.occ.append(0.0)

    self.box = [float(x)*10. for x in lines[self.atcounter+2].split()]
    if debug:
      return len(self.atnum),self.atcounter
# }}}

# {{{ writegmxlines
  def writegmxlines(self):
    lines = []
    lines.append(self.title)
    lines.append(self.atcounter)
    for i in range(self.atcounter):
      lines.append(self.gmxlineout(i))
    lines.append(self.box)
    return lines
# }}}

# {{{ writegmx
  def writegmx(self,file_out):
# check if passed filename string or a file descriptor
    if type(file_out) == type(sys.stdout):
      out = file_out
    else:
      out = open(file_out,'w')

    out.write(self.title + "\n")
    out.write("%5d\n" % self.atcounter)

# use sequential residue and atom numbers, not what was read in
# check if residue number is same as previous one, increment counter if not
# initial with really absurd residue number
    last_resnum = -999999
    self.rescounter = 0
    for i in range(self.atcounter):
      if self.resnum[i] != last_resnum:
        last_resnum = self.resnum[i]
        self.rescounter += 1
      out.write(self.gmxlineout(i) + "\n")
    for j in range(len(self.box)):
      out.write("%10.5f" % (self.box[j]/10.))
    out.write("\n")
# }}}

# {{{ gmxlineout
  def gmxlineout(self,i):
    return '%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f' % (self.rescounter,self.resname[i],self.atname[i],i+1,
               self.coord[i][0]/10.,self.coord[i][1]/10.,self.coord[i][2]/10.,
               self.velocity[i][0]/10.,self.velocity[i][1]/10.,self.velocity[i][2]/10.)
# }}}

# {{{ readPDB
  def readPDB(self,file_in,debug=0):
    """
    Read a PDB-format file into Protein class object that has been initialized by
    p = MyPDB.Protein()
    """
# check if passed filename string or a file descriptor
#    if type(file_in) == type(sys.stdin):
#      readfile = file_in
#    else:
#      readfile = open(file_in,'r')
#    lines = readfile.readlines()
    try:
      lines = file_in.readlines()
    except AttributeError:
      if file_in[-3:] == '.gz':
        lines = gzip.open(file_in).readlines()
      else:
        lines = open(file_in).readlines()

    self.readPDBlines(lines,debug)
# }}}

# {{{ readPDBlines
  def readPDBlines(self,lines,debug=0):
    """
    Reads a list of PDB-format file lines in to the Protein class object.
    Thus can be called by another routine that already has the lines in a list.
    Returns the number of atoms read in.
    """
    i = 0
#    j = 0
    atom_hetatm = re.compile('(ATOM  |HETATM)')
    head = re.compile('^(HEADER|COMPND|REMARK|SEQRES|CRYST1|SCALE|ORIG)')
    title = re.compile('^TITLE')
    foot = re.compile('(CONECT |TER  |MASTER|END)')
    element = re.compile('[A-Za-z ][A-Za-z]')
    for line in lines:
      #print(line)
      if atom_hetatm.match(line):
        # not sure why I wrote it this way rather than using ".strip()"
        #        line = line[:-1]
        line = line.strip()
        self.allatom_line.append(line)
        self.label.append(line[0:6])
        self.atnum.append(int(line[6:12]))
        self.atname.append(line[12:16])
        self.atalt.append(line[16:17])
        self.resname.append(line[17:21])
        self.chain.append(line[21])
        self.resnum.append(int(line[22:26]))
        self.resext.append(line[27])
        self.coord.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
        self.occ.append(float(line[54:60]))
        self.b.append(float(line[60:66]))
        self.velocity.append(0.0)
        if element.match(line[76:78]):
          self.elem.append(line[76:78])
        else:
          self.elem.append(line[12:14])
        try:
          self.mass.append(Mass_table[self.elem[i]])
        except KeyError:
          self.mass.append(0.0)
          if debug:
            sys.stderr.write('Problem determining mass of element: %s\n' % self.elem[i])
        self.atcounter += 1

# make list of only C-alpha's for speeding up certain calculations
        if self.elem[i] == ' C' and self.atname[i] == ' CA ':
          self.ca_line.append(line)
          self.ca_x.append([self.coord[i][0],self.coord[i][1],self.coord[i][2]])

#          j = j + 1

        i = i + 1
      elif head.match(line):
        self.header.append(line[:-1])
      elif foot.match(line):
        self.footer.append(line[:-1])
      elif title.match(line):
        self.title = self.title + line[10:-1].strip() + " "

    self.sequence = self.get_sequence()

    if debug:
      return len(self.atnum),self.atcounter
# }}}


# writePDB {{{
  def writePDB(self,file_out,header=1,hetatm=1):
    """
    Saves the Protein class object to a PDB-format file
      if header = 0, no header (REMARK etc.) or footer lines are written
      if hetatm = 0, no hetero atoms are written
    """
# check if passed filename string or a file descriptor
    if type(file_out) == type(sys.stdout):
      out = file_out
    else:
      out = open(file_out,'w')

    if header == 1:
      out.write('%-10s%s\n' % ('TITLE',self.title))
      for i in range(len(self.header)):
        out.write("%s\n" % self.header[i])

    # put this "if" outside the loop (for speed?)
    if hetatm == 1:
      for i in range(len(self.resnum)):
        self.writePDBline(out,i)
    else:
      for i in range(len(self.resnum)):
        if self.label[i] == 'ATOM  ':
          self.writePDBline(out,i)

    if header == 1:
      for i in range(len(self.footer)):
        out.write("%s\n" % self.footer[i])
    out.close()
# }}}

# {{{ writePDBline
  def writePDBline(self,FD,i):
    """
    Writes a single line of data in the PDB-format
    called by writePDB
    """
# removed the writing of self.atnum in favour of "i" == sequential atom numbers
#    FD.write('%-6s%5i %4s%1s%4s%1s%4i%1s   %8.3f%8.3f%8.3f%6.2f%6.2f%10s%2s\n' %
#    (self.label[i],self.atnum[i],self.atname[i],self.atalt[i],
#     self.resname[i],self.chain[i], self.resnum[i],self.resext[i],
#     self.coord[i][0],self.coord[i][1],self.coord[i][2],self.occ[i],self.b[i],blank,self.elem[i]))
    FD.write('%-6s%5i %-4s%1s%-4s%1s%4i%1s   %8.3f%8.3f%8.3f%6.2f%6.2f%10s%2s\n' %
    (self.label[i],i+1,self.atname[i],self.atalt[i],
     self.resname[i],self.chain[i], self.resnum[i],self.resext[i],
     self.coord[i][0],self.coord[i][1],self.coord[i][2],self.occ[i],self.b[i],blank,self.elem[i]))
# }}}

# {{{ make_dist_matrix
  def make_dist_matrix(self):
    """
    Creates an attribute of the class Protein that contains the
    complete (!) inter-atomic distance matrix

    Warning: this can take some time!!!

    """
    print("Creating distance_matrix")
    self.dist_matrix=N.zeros((self.atcounter,self.atcounter))
    print("... zeroed %dx%d array" % (self.atcounter,self.atcounter))
    for i in range(self.atcounter):
      if i%100 == 0:
        print("Creating %dth row." % i)
      for j in range(i,self.atcounter):
        self.dist_matrix[i,j] = dist(self.coord[i],self.coord[j])
# }}}

# {{{ make_CA_dist_matrix
  def make_CA_dist_matrix(self):
    """
    Creates an attribute of the class Protein
    that contains the inter-atomic C-alpha distance matrix

    Note that if you import the NumTut view method, you can view the intra-Calpha
    matrix with view(p.CA_dist_matrix), where "p" is the name of your Protein instance.
    Black is close and white is far away.  Cool!

    Calculate a difference distance matrix:
    create p1.CA_dist_matrix and p2.CA_dist_matrix for two proteins p1 and p2
    then,

      diff = p1.CA_dist_matrix - p2.CA_dist_matrix

    but first check that p1.CA_dist_matrix.shape == p2.CA_dist_matrix.shape

    you can do e.g.:

      diff = p1.CA_dist_matrix[0:320,0:320] - p2.CA_dist_matrix[3:323,3:323]

    """
    print("Creating C-alpha distance_matrix")
    ca_count = len(self.ca_x)
    self.CA_dist_matrix=N.zeros((ca_count,ca_count))
    print("... zeroed %dx%d array" % (ca_count,ca_count))
    for i in range(ca_count):
      if i%20 == 0:
        print("Creating %dth row." % i)
      for j in range(i,ca_count):
        self.CA_dist_matrix[i,j] = dist(self.ca_x[i],self.ca_x[j])
        self.CA_dist_matrix[j,i] = self.CA_dist_matrix[i,j]
# }}}

# {{{ make_contact_map
  def make_contact_map(self,close_limit=2,cutoff=8):
    """
    Creates an attribute of the class Protein
    that contains the inter-atomic C-alpha contact map

    close_limit represents a limit such that C-alphas whose indices are less than or equal to
    close_limit apart give rise to a Cmap value of 0

    cutoff is the distance cutoff to consider C-alphas to be in contact

    """
    print("Creating C-alpha distance_matrix")
    ca_count = len(self.ca_x)
    self.CA_cmap=N.zeros((ca_count,ca_count))
    print("... zeroed %dx%d array" % (ca_count,ca_count))
    for i in range(ca_count):
      if i%20 == 0:
        print("Creating %dth row." % i)
      for j in range(i,ca_count):
        if abs(i-j) > close_limit:
          d = dist(self.ca_x[i],self.ca_x[j])
          if d <= cutoff:
            self.CA_cmap[i,j] = 1
        self.CA_cmap[j,i] = self.CA_cmap[i,j]
# }}}

# {{{ avg
  def avb(self):
    self.make_backbone_id_list()
    self.protein_b = []
    self.backbone_b = []
    self.sidechain_b = []
    self.solvent_b = []
    self.other_b = []
    self.all_b = []
    res_counter = 0
    for i in range(self.atcounter):
      self.all_b.append(self.b[i])

    for i in self.backbone_id:
      self.backbone_b.append(self.b[i])
      self.protein_b.append(self.b[i])
    for i in self.sidechain_id:
      self.sidechain_b.append(self.b[i])
      self.protein_b.append(self.b[i])
    for i in self.solvent_id:
      self.solvent_b.append(self.b[i])
    for i in self.other_id:
      self.other_b.append(self.b[i])
# }}}

# {{{ make_backbone_id_list
  def make_backbone_id_list(self):
    self.CA_id = []
    self.backbone_id = []
    self.sidechain_id = []
    self.solvent_id = []
    self.other_id = []
    for i in range(self.atcounter):
      if self.resname[i].strip() in amino_acid_names_list:
        if self.atname[i] == ' CA ':
          self.CA_id.append(i)
        if self.atname[i] in [' N  ',' CA ',' C  ',' O  ']:
          self.backbone_id.append(i)
        else:
          self.sidechain_id.append(i)
      elif self.resname[i].strip() in solvent_names_list:
        self.solvent_id.append(i)
      else:
        self.other_id.append(i)
# }}}

# {{{ make_chain_backbone_id_list
  def make_chain_backbone_id_list(self):
    self.chain_CA_id = {}
    self.chain_backbone_id = {}
    self.chain_sidechain_id = {}
    self.chain_solvent_id = {}
    self.chain_other_id = {}
    for i in range(self.atcounter):
      if self.resname[i].strip() in amino_acid_names_list:
        if self.atname[i] == ' CA ':
          self.chain_CA_id.setdefault(self.chain[i], []).append(i)
        elif self.atname[i] in [' N  ',' CA ',' C  ',' O  ']:
          self.chain_backbone_id.setdefault(self.chain[i], []).append(i)
        else:
          self.chain_sidechain_id.setdefault(self.chain[i], []).append(i)
      elif self.resname[i].strip() in solvent_names_list:
        self.chain_solvent_id.setdefault(self.chain[i], []).append(i)
      else:
        self.chain_other_id.setdefault(self.chain[i], []).append(i)
# }}}

# {{{ phipsi
  def phipsi(self):
    self.backbone = []
    self.phi = []
    self.psi = []
    res_counter = 0
    for i in range(self.atcounter):
#debug#      print(i,self.atname[i],res_counter)
      if self.atname[i] == ' N  ':
        res_counter += 1
# this might not really be necessary as phipsi only uses N,CA,C
#        if self.to_atnum(self.chain[i],self.resnum[i],self.resname[i],' O  '):
#          coord_O = self.coord[self.to_atnum(self.chain[i],self.resnum[i],self.resname[i],' O  ')]
#        elif self.to_atnum(self.chain[i],self.resnum[i],self.resname[i],' O1 '):
#          coord_O = self.coord[self.to_atnum(self.chain[i],self.resnum[i],self.resname[i],' O1 ')]
#        elif self.to_atnum(self.chain[i],self.resnum[i],self.resname[i],' O2 '):
#          coord_O = self.coord[self.to_atnum(self.chain[i],self.resnum[i],self.resname[i],' O2 ')]
#        elif self.to_atnum(self.chain[i],self.resnum[i],self.resname[i],' OXT'):
#          coord_O = self.coord[self.to_atnum(self.chain[i],self.resnum[i],self.resname[i],' OXT')]
        self.backbone.append(
             (res_counter,
              self.chain[i],
              self.resnum[i],
              self.resname[i],
              self.coord[i],
              self.coord[self.to_atnum(self.chain[i],self.resnum[i],self.resname[i],' CA ')],
              self.coord[self.to_atnum(self.chain[i],self.resnum[i],self.resname[i],' C  ')],
#              coord_O
#              self.coord[self.to_atnum(self.chain[i],self.resnum[i],self.resname[i],' O  ')]
              )
             )

    b,c,d = map(N.asarray,self.backbone[0][4:7])
    e = N.asarray(self.backbone[1][4])
    tmppsi = torsion(b,c,d,e)
    if tmppsi > 180:
      tmppsi = tmppsi - 360.

    self.phi.append(180.0)
    self.psi.append(tmppsi)

    for i in range(1,len(self.backbone)-1):
      a = N.asarray(self.backbone[i-1][6])
      b,c,d = map(N.asarray,self.backbone[i][4:7])
      e = N.asarray(self.backbone[i+1][4])
# check that in same chain or not
      if self.backbone[i-1][1] == self.backbone[i][1] and self.backbone[i][1] == self.backbone[i+1][1]:
        #        print("a,b,c,d",a,b,c,d)
        tmpphi = torsion(a,b,c,d)
        #print("tmpphi",tmpphi)
        tmppsi = torsion(b,c,d,e)
      else:
        if self.backbone[i-1][1] == self.backbone[i][1]:
          tmpphi = torsion(a,b,c,d)
          tmppsi = 180.
        elif self.backbone[i][1] == self.backbone[i+1][1]:
          tmpphi = 180.
          tmppsi = torsion(b,c,d,e)
      if tmpphi > 180:
        tmpphi = tmpphi - 360.
      if tmppsi > 180:
        tmppsi = tmppsi - 360.

      self.phi.append(tmpphi)
      self.psi.append(tmppsi)

    i += 1
    #print(i,self.backbone[i-1],self.backbone[i])
    a = N.asarray(self.backbone[i-1][6])
    b,c,d = map(N.asarray,self.backbone[i][4:7])
    tmpphi = torsion(a,b,c,d)
    if tmpphi > 180:
      tmpphi = tmpphi - 360.
    self.phi.append(tmpphi)
    self.psi.append(180.0)
# }}}

# {{{ find_neighbours
  def find_neighbours(self,atnum,cutoff):
    """
    Find neighbours within 'cutoff' distance of the atom numbered 'atnum'
    Use the function 'to_atnum(chain,resnum,resname,atname)' to get an atom
    number to feed to this function.

    A dictionary is then created with keys of the central defined atom and values
    of the atom numbers within the cutoff distance.

    e.g.
      p=MyPDB.Protein()
      p.readPDB('test.pdb')
      p.find_neighbours(p.to_atnum('A',105,'cys','sg'),4.0)
      [[' ', 75, 'GLN ', ' CD '], [' ', 75, 'GLN ', ' OE1'], [' ', 75, 'GLN ', ' NE2'],
      [' ', 226, 'TYR ', ' CD2'], [' ', 286, 'HIS ', ' ND1'], [' ', 286, 'HIS ', ' CG '],
      [' ', 286, 'HIS ', ' CE1']]
      p.neighbours
      {668: [627, 628, 629, 1842, 2300, 2301, 2305]}

    """

    cutoff2 = cutoff**2
    self.neighbours[atnum] = []
    for i in range(self.atcounter):
      if ((self.chain[atnum] != self.chain[i] or self.resnum[atnum] != self.resnum[i])
          and distx2(self.coord[i],self.coord[atnum]) <= cutoff2):
        self.neighbours[atnum].append(i)
    atom_list = self.neighbours[atnum]
    full_atomlist = []
    for i in atom_list:
      full_atomlist.append([self.chain[i],self.resnum[i],self.resname[i],self.atname[i]])
    return full_atomlist
# }}}

# {{{ to_atnum
  def to_atnum(self,chain,resnum,resname,atname):
    """
    Convert an atom specified by chain,residue number, residue name, and atom name
    to an atom number in the "self" object -- useful for passing on to find_neighbours and
    other functions.

    """
#    print(chain,resnum,resname,elem,atname,self.atcounter)
    chain = chain.strip().upper()
    resname = resname.strip().upper()
#    elem = elem.strip().upper()
    atname = atname.strip().upper()
    for i in range(self.atcounter):
#      if self.resnum[i] == resnum:
#        print(self.chain[i],self.resnum[i],self.resname[i],self.elem[i],self.atname[i],i)
      if ((chain == '' or self.chain[i].strip() == chain) and
          self.resnum[i] == resnum and
          (resname == '' or self.resname[i].strip() == resname) and
#          (elem == '' or self.elem[i].strip() == elem) and
          (atname == '' or self.atname[i].strip() == atname)):
#        print("Found it", i)
        return i
    return None
# {{{

# {{{ renumber
  def renumber(self,start):
    lastresnum = -9999
    lastchain = ' '
    lastresname = ' '
    lastext = ' '
#    idatom = 0
    idres = start - 1
    icount = 0

    for i in range(len(self.resnum)):
      if i == 0:
        icount += 1
        lastchain = self.chain[i]
        lastresname = self.resname[i]
        lastresnum = self.resnum[i]
        lastext = self.resext[i]
        self.resnum[i] = start
        idres += 1
      else:
        if (self.chain[i] != lastchain or lastresnum != self.resnum[i] or
            lastresname != self.resname[i] or lastext != self.resext[i]):
          icount += 1
          idres += 1
          lastchain = self.chain[i]
          lastresname = self.resname[i]
          lastresnum = self.resnum[i]
          self.resnum[i] = idres
        else:
          self.resnum[i] = idres
# }}}

# {{{ get_sequence
  def get_sequence(self):
    """
    Returns a dictionary as well as assigning the sequence to the dictionary self.sequence.
    This function is called by the readPDBlines function, such that after reading a PDB file:

      P=MyPDB.Protein()
      P.readPDB('test1.pdb')

    P.sequence would contain the sequence as a dictionary whose keys
    are the chain IDs and the values are a list of the sequence in
    3-letter code.

    With the sequence in this form, one can then do

      import seq_convert
      P.one_letter = {}
      for chain in P.sequence.keys():
        P.one_letter[chain] = seq_convert.seq3_to_seq1(P.sequence[chain])

    """
    for i in range(self.atcounter):
      if self.name[i] == ' CA ' and (self.atalt[i] == ' ' or self.atalt[i] == 'A'):
        if self.chain[i] in self.sequence:
          self.sequence[self.chain[i]].append(self.resn[i])
        else:
          self.sequence[self.chain[i]] = [self.resn[i]]
    return self.sequence
# }}}

#class Trace:
#  def __init__(self):
#    self.trace = []

  def make_trace(self):
    self.trace = []
    for i in range(1,len(self.ca_x)):
      self.trace.append([self.ca_x[i-1], self.ca_x[i]])
      print(i, self.ca_x[i-1], self.ca_x[i])

if __name__ == '__main__':
  if len(sys.argv) > 1:
    file_in=sys.argv[1]
  else:
    file_in=raw_input('Enter file name: ')

  p=Protein()
  p.readPDB(file_in)
#  for i in range(len(a.ca_line)):
#    print(i, a.ca_line[i])

