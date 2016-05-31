from chempy import Atom,Bond
from chempy.models import Indexed
from copy import deepcopy
import time

def set_molecule_property(code, x,a):
    a.molecule_properties.append([code, x])

#
# dotmae (.mae) file reader

import re,string

# regular expressions for parsing

quoted_re = re.compile(r'"((?:\\.|[^"])*)"')
comment_re = re.compile(r'\s#[^#]*#')
_label_fmt_re = re.compile(r'%[A-Z]{2}')

_crlf_re = re.compile(r"\r\n|\r")

_p2sym = { 1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N',
    8: 'O', 9: 'F', 10: 'Ne', 11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si',
    15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca', 21: 'Sc',
    22: 'Ti', 23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28:
    'Ni', 29: 'Cu', 30: 'Zn', 31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se',
    35: 'Br', 36: 'Kr', 37: 'Rb', 38: 'Sr', 39: 'Y', 40: 'Zr', 41:
    'Nb', 42: 'Mo', 43: 'Tc', 44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag',
    48: 'Cd', 49: 'In', 50: 'Sn', 51: 'Sb', 52: 'Te', 53: 'I', 54:
    'Xe', 55: 'Cs', 56: 'Ba', 57: 'La', 58: 'Ce', 59: 'Pr', 60: 'Nd',
    61: 'Pm', 62: 'Sm', 63: 'Eu', 64: 'Gd', 65: 'Tb', 66: 'Dy', 67:
    'Ho', 68: 'Er', 69: 'Tm', 70: 'Yb', 71: 'Lu', 72: 'Hf', 73: 'Ta',
    74: 'W', 75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt', 79: 'Au', 80:
    'Hg', 81: 'Tl', 82: 'Pb', 83: 'Bi', 84: 'Po', 85: 'At', 86: 'Rn',
    87: 'Fr', 88: 'Ra', 89: 'Ac', 90: 'Th', 91: 'Pa', 92: 'U', 93:
    'Np', 94: 'Pu', 95: 'Am', 96: 'Cm', 97: 'Bk', 98: 'Cf', 99: 'Es',
    100: 'Fm', 101: 'Md', 102: 'No', 103: 'Lr', 104: 'Rf', 105: 'Db',
    106: 'Sg', 107: 'Bh', 108: 'Hs', 109: 'Mt', 110: 'Ds', 111: 'Rg',
    }

_ss2ss = { 0: 'L', 1: 'H', 2: 'S' }

# from schrodinger.structutils.color import get_rgb_from_color_index
# _colors = [(i, get_rgb_from_color_index(i)) for i in xrange(1, 160)]
_colors = [(1,(0,0,0)),(2,(160,160,160)),(3,(0,0,180)),
           (4,(30,30,225)),(5,(100,100,225)),(6,(112,219,147)),(7,(173,234,234)),
           (8,(0,255,127)),(9,(0,100,0)),(10,(30,225,30)),(11,(50,204,50)),
           (12,(153,204,30)),(13,(225,225,30)),(14,(234,130,50)),
           (15,(142,35,107)),(16,(225,30,30)),(17,(255,152,163)),
           (18,(234,173,234)),(19,(225,30,225)),(20,(159,95,159)),
           (21,(255,255,255)),(22,(165,42,42)),(23,(225,127,80)),
           (24,(105,105,105)),(25,(225,193,37)),(26,(225,105,180)),
           (27,(107,142,35)),(28,(205,133,63)),(29,(160,82,45)),(30,(70,130,180)),
           (31,(216,191,216)),(32,(245,222,179)),(33,(7,7,255)),(34,(15,15,255)),
           (35,(23,23,255)),(36,(31,31,255)),(37,(39,39,255)),(38,(47,47,255)),
           (39,(55,55,255)),(40,(63,63,255)),(41,(71,71,255)),(42,(79,79,255)),
           (43,(87,87,255)),(44,(95,95,255)),(45,(103,103,255)),
           (46,(111,111,255)),(47,(119,119,255)),(48,(127,127,255)),
           (49,(135,135,255)),(50,(143,143,255)),(51,(151,151,255)),
           (52,(159,159,255)),(53,(167,167,255)),(54,(175,175,255)),
           (55,(183,183,255)),(56,(191,191,255)),(57,(199,199,255)),
           (58,(207,207,255)),(59,(215,215,255)),(60,(223,223,255)),
           (61,(231,231,255)),(62,(239,239,255)),(63,(247,247,255)),
           (64,(255,255,255)),(65,(255,7,7)),(66,(255,15,15)),(67,(255,23,23)),
           (68,(255,31,31)),(69,(255,39,39)),(70,(255,47,47)),(71,(255,55,55)),
           (72,(255,63,63)),(73,(255,71,71)),(74,(255,79,79)),(75,(255,87,87)),
           (76,(255,95,95)),(77,(255,103,103)),(78,(255,111,111)),
           (79,(255,119,119)),(80,(255,127,127)),(81,(255,135,135)),
           (82,(255,143,143)),(83,(255,151,151)),(84,(255,159,159)),
           (85,(255,167,167)),(86,(255,175,175)),(87,(255,183,183)),
           (88,(255,191,191)),(89,(255,199,199)),(90,(255,207,207)),
           (91,(255,215,215)),(92,(255,223,223)),(93,(255,231,231)),
           (94,(255,239,239)),(95,(255,247,247)),(96,(255,255,255)),
           (97,(44,60,60)),(98,(0,0,128)),(99,(34,139,34)),(100,(0,255,255)),
           (101,(220,20,60)),(102,(210,127,36)),(103,(192,192,192)),
           (104,(224,224,224)),(105,(255,30,0)),(106,(255,60,0)),(107,(255,90,0)),
           (108,(255,120,0)),(109,(255,150,0)),(110,(255,180,0)),
           (111,(255,210,0)),(112,(255,240,0)),(113,(255,255,0)),
           (114,(204,255,0)),(115,(153,255,0)),(116,(102,255,0)),(117,(51,255,0)),
           (118,(0,204,51)),(119,(0,153,102)),(120,(0,102,153)),(121,(0,51,204)),
           (122,(0,0,255)),(123,(15,0,255)),(124,(30,0,255)),(125,(45,0,255)),
           (126,(60,0,255)),(127,(75,0,255)),(128,(90,0,255)),(129,(255,0,31)),
           (130,(223,63,0)),(131,(191,127,25)),(132,(255,191,51)),
           (133,(223,223,0)),(134,(191,255,51)),(135,(159,255,0)),
           (136,(127,191,25)),(137,(95,255,51)),(138,(63,191,76)),
           (139,(31,255,102)),(140,(0,191,127)),(141,(159,255,191)),
           (142,(159,210,255)),(143,(255,223,159)),(144,(191,159,255)),
           (145,(95,63,255)),(146,(0,31,191)),(147,(191,0,191)),(148,(255,0,159)),
           (149,(0,255,51)),
           (150,    (   0,   255,       102)),
           (151,    (   0,   255,       153)),
           (152,    (   0,   255,       204)),
           (153,    (   0,   204,       255)),
           (154,    (   0,   153,       255)),
           (155,    (   0,   102,       255)),
           (156,    (   0,    51,       255)),
           (157,    ( 120,     0,       255)),
           (158,    ( 180,     0,       255)),
           ]

_color_dict = {}

# representation mappings
# see pymol.constants.repres or layer1/Rep.h

REP_STICKS     = 0x000001
REP_SPHERES    = 0x000002
REP_SURFACE    = 0x000004
REP_LABELS     = 0x000008
REP_NB_SPHERES = 0x000010
REP_CARTOON    = 0x000020
REP_RIBBON     = 0x000040
REP_LINES      = 0x000080
REP_MESH       = 0x000100
REP_DOTS       = 0x000200
REP_DASHES     = 0x000400
REP_NONBONDED  = 0x000800
REP_CELL       = 0x001000
REP_CGO        = 0x002000
REP_CALLBACK   = 0x004000
REP_EXTENT     = 0x008000
REP_SLICE      = 0x010000
REP_ANGLES     = 0x020000
REP_DIHEDRALS  = 0x040000
REP_ELLIPSOIDS = 0x080000
REP_VOLUME     = 0x100000

REPS_ATOM      = REP_SPHERES | REP_NONBONDED
REPS_BOND      = REP_STICKS | REP_LINES
REPS_CARTOON   = REP_CARTOON | REP_RIBBON

# mmshare/mmlibs/mmct/mmctg.h:
# enum MM_CTAtomStyle
MMCT_ATOM_NOSTYLE    = 0
MMCT_ATOM_CIRCLE     = 1
MMCT_ATOM_CPK        = 2
MMCT_ATOM_BALLNSTICK = 3
# enum MM_CTBondStyle
MMCT_BOND_NOSTYLE    = 0
MMCT_BOND_WIRE       = 1
MMCT_BOND_TUBE       = 2
MMCT_BOND_BALLNSTICK = 3
# mmshare/mmlibs/mmct/mmct.h
# enum MM_CTRibbonStyle
MMCT_RIBBON_STYLE_NONE      = 0  #
MMCT_RIBBON_STYLE_CARTOON   = 1  # auto
MMCT_RIBBON_STYLE_RIBBON    = 2  # auto
MMCT_RIBBON_STYLE_TUBE      = 3  # tube, cartoon_tube_radius=0.30
MMCT_RIBBON_STYLE_THINTUBE  = 4  # tube, cartoon_tube_radius=0.15
MMCT_RIBBON_STYLE_CURVELINE = 5  # ribbon, ribbon_sampling=5
MMCT_RIBBON_STYLE_CALINE    = 6  # ribbon, ribbon_as_cylinders=0
MMCT_RIBBON_STYLE_CATUBE    = 7  # ribbon, ribbon_as_cylinders=1

_atom_rep = {
    MMCT_ATOM_NOSTYLE:    REP_NONBONDED,
    MMCT_ATOM_CIRCLE:     REP_NONBONDED,
    MMCT_ATOM_CPK:        REP_SPHERES,
    MMCT_ATOM_BALLNSTICK: REP_SPHERES,
}

_atom_set = {
    MMCT_ATOM_CPK:        [("sphere_scale", 0.85)],
    MMCT_ATOM_BALLNSTICK: [("sphere_scale", 0.18)],
}

_cartoon_rep = {
    1: REP_CARTOON,
    2: REP_CARTOON,
    3: REP_CARTOON,
    4: REP_CARTOON,
    5: REP_RIBBON,
    6: REP_RIBBON,
    7: REP_RIBBON,
}

_cartoon_set = {
    3: [("cartoon_tube_radius", 0.3)],
    4: [("cartoon_tube_radius", 0.15)],
    5: [("ribbon_sampling", 5)],
    6: [("ribbon_as_cylinders", 0)],
    7: [("ribbon_as_cylinders", 1)],
}

_bond_rep = {
    MMCT_BOND_NOSTYLE:    0,
    MMCT_BOND_WIRE :      REP_LINES,
    MMCT_BOND_TUBE:       REP_STICKS,
    MMCT_BOND_BALLNSTICK: REP_STICKS,
}

_bond_set = {
    MMCT_BOND_BALLNSTICK: [("stick_radius", 0.15)],
}

def prime_colors():
    for (key,rgb) in _colors:
        trgb = eval("0x%02x%02x%02x"%rgb)
        _color_dict[key] = trgb
        _color_dict[-key] = trgb        
        _color_dict[-1000-key] = trgb
        
class MAEReaderException(Exception):

    def __init__(self,value):
        self.value = value
    def __str__(self):
        return `self.value`

def default_atom_handler(c,x,a):
    a.atom_properties[c] = x
    return 0
    
def handler_i_m_representation(c, x, a):
    x = int(x)
    a.visible = _atom_rep.get(x, 0)
    try:
        a.__dict__.update(_atom_set[x])
    except KeyError:
        pass

def handler_i_m_ribbon_style(c, x, a):
    x = int(x)
    a.visible |= _cartoon_rep.get(x, 0)
    try:
        a.__dict__.update(_cartoon_set[x])
    except KeyError:
        pass

def handler_atom_name(c, x, a):
    if not a.name and x.strip():
        a.name = x

def handler_i_m_from_rep(x, b, key="from_rep"):
    x = int(x)
    setattr(b, key, x)
    try:
        b.__dict__.update(_bond_set[x])
    except KeyError:
        pass

def handler_i_m_to_rep(x, b):
    handler_i_m_from_rep(x, b, "to_rep")

def quotedsplit(handle):
    '''
    regex based string split. Handles double quotes and escapes inside
    double quotes.
    '''
    if not isinstance(handle, basestring):
        handle = handle.read()

    reduced = quoted_re.sub('"', handle)
    reduced = comment_re.sub(' ', reduced)
    pos = 0

    for tok in reduced.split():
        if tok == '"':
            match = quoted_re.search(handle, pos)
            pos = match.end()
            tok = match.group(1).replace('\\', '')
        yield tok

class MAEReader:

    def __init__(self,mimic=1, load_properties='',load_atom_properties=''):
        self.full_count = 0
        self.mimic = int(mimic)
        self.load_properties = load_properties
        self.load_all_atom_properties = load_atom_properties == '*'
        self.load_atom_properties = set(load_atom_properties.split(' '))
        prime_colors()
        
    def open_brace(self,tok):
        if self.stack[-1:] == ['f_m_ct']:
            self.full_bond_block = None
            self.full_atom_block = None
        self.meta_stack.append( (self.stack,self.spec) )
        self.stack = []
        self.spec = []
        
    def triple_colon(self,tok):
        if self.spec == []:
            self.spec = self.stack
            self.stack = []

    def parse_atom_block(self,count,spec,data,full):
        handler = {
            'i_m_mmod_type' : lambda c,x,a:setattr(a,'numeric_type',int(x)),
            'r_m_x_coord' : lambda c,x,a:a.coord.__setitem__(0,float(x)),
            'r_m_y_coord' : lambda c,x,a:a.coord.__setitem__(1,float(x)),
            'r_m_z_coord' : lambda c,x,a:a.coord.__setitem__(2,float(x)),
            'i_m_residue_number' : lambda c,x,a:setattr(a,'resi',x),
            's_m_insertion_code' : lambda c,x,a:setattr(a,'ins_code',x),
#            's_m_mmod_res' : lambda c,x,a:setattr(a,'mmod_res',x),
            's_m_chain_name' : lambda c,x,a:setattr(a,'chain',x),
            'i_m_color' : lambda c,x,a:setattr(a,'trgb',_color_dict.get(int(x),16777215)),
            'r_m_charge1' : lambda c,x,a:setattr(a,'partial_charge',float(x)),
#            'r_m_charge2' : lambda c,x,a:setattr(a,'',x),
            's_m_pdb_residue_name' : lambda c,x,a:setattr(a,'resn',x),
            's_m_pdb_atom_name' : handler_atom_name,
#            's_m_grow_name' : lambda c,x,a:setattr(a,'',x),
            'i_m_atomic_number' : lambda c,x,a:setattr(a,'symbol',_p2sym.get(int(x),'X')),
            'i_m_formal_charge' : lambda c,x,a:setattr(a,'formal_charge',int(x)),
            's_m_atom_name' : handler_atom_name,
            'i_m_secondary_structure' : lambda c,x,a:setattr(a,'ss',_ss2ss.get(int(x),'L')),
            's_m_label_format' : lambda c,x,a:setattr(a,'label_format',x),
            'i_m_label_color' : lambda c,x,a:(
            setattr(a,'label_trgb',_color_dict.get(int(x),16777215))),
            's_m_label_user_text' : lambda c,x,a:setattr(a,'label_user_text',x),
            'i_m_ribbon_style' : handler_i_m_ribbon_style,
            'i_m_representation' : handler_i_m_representation,
            'i_m_ribbon_color' : lambda c,x,a:(
            setattr(a,'ribbon_trgb',_color_dict.get(int(x),16777215)),
            setattr(a,'cartoon_trgb',_color_dict.get(int(x),16777215))),                      
#            setattr(a,'cartoon_trgb',_rib_color(int(x)))),
            'i_m_visibility' : lambda c,x,a:setattr(a,'visibility',int(x)),
            'r_m_pdb_tfactor' : lambda c,x,a:setattr(a,'b',float(x)),
#            'i_m_pdb_convert_problem' : lambda c,x,a:setattr(a,'',x),
            }
        if not self.mimic:
            handler['i_m_label_color']=lambda c,x,a:None
            handler['i_m_ribbon_color']=lambda c,x,a:None
            
        fn_list = [ ['id', lambda c,x,a:setattr(a,'id',int(x)) ] ]
        if self.load_all_atom_properties:
            for code in spec:
                fn_list.append([ code, handler.get(code, default_atom_handler) ])
        else:
            for code in spec:
                if code in self.load_atom_properties:
                    fn_list.append([ code, handler.get(code, default_atom_handler) ])
                else:
                    fn_list.append([ code, handler.get(code, lambda c,x,a: 0) ])
            
        data.reverse()
        if full != None:
            count = len(data)/len(fn_list)
            atom_list = deepcopy(full)
        else:
            atom_list = []
        for a in xrange(count):
            if full == None:
                at = Atom()
                at.coord = [0.0,0.0,0.0]
                at.visible = REP_NONBONDED
                at.hetatm = 0
                atom_list.append(at)
            else:
                at = atom_list[int(data[-1])-1] # trusting index
            if 'i_m_representation' in spec: # representation information provided
                at.visible = 0

            for code, fn in fn_list: 
                apply(fn,(code, data.pop(),at))

            if hasattr(at,'label_format'):
                if at.label_format != "%OF":
                    get = {
                        "%UT": at.label_user_text,
                        "%RT": at.resn,
                        "%RN": at.resi,
                        "%CH": at.chain,
                        "%EL": at.symbol,
                        "%NU": str(at.id),
                        "%FC": "%+d" % at.formal_charge if at.formal_charge else "",
                    }.get
                    label = _label_fmt_re.sub(lambda m: get(m.group(), ''), at.label_format)
                    if len(label.replace(" ","")):
                        at.visible |= REP_LABELS
                        at.label = label
                    elif hasattr(at,'label_trgb'):
                        del at.label_trgb
                del at.label_format
                
        for atom in atom_list:
            if hasattr(atom,'ins_code'):
                atom.resi=atom.resi+atom.ins_code
                del atom.ins_code
                
        return atom_list
    
    def parse_bond_block(self,count,spec,data,full):
        handler = {
            'i_m_from' : lambda x,a:a.index.__setitem__(0,int(x)-1),
            'i_m_to' : lambda x,a:a.index.__setitem__(1,int(x)-1),
            'i_m_order' : lambda x,a:setattr(a,'order',int(x)),
            'i_m_from_rep': handler_i_m_from_rep,
            'i_m_to_rep': handler_i_m_to_rep,
            'b_m_thin_bond': lambda x,a: int(x) and setattr(a,'stick_radius', 0.08),
            }
        fn_list = [ lambda x,a:setattr(a,'id',int(x)) ]
        for code in spec:
            fn_list.append(handler.get(code,lambda x,a:0))
        if full != None:
            count = len(data)/len(fn_list)
            bond_list = deepcopy(full)
        else:
            bond_list = []
        data.reverse()
        try:
            for a in xrange(count):
                if full == None:
                    bd = Bond()
                    bd.index = [0,0]
                    bond_list.append(bd)
                else:
                    bd = bond_list[int(data[-1])-1] # trusting index
                for fn in fn_list:
                    apply(fn,(data.pop(),bd))
        except IndexError:
            pass
        return bond_list

    def parse_ct_block(self,spec,data):
        handler = {
            's_m_title'   : lambda code,x,a:setattr(a,'title',x),
            'chempy_atoms' : lambda code,x,a:setattr(a,'atom',x),
            'chempy_bonds' : lambda code,x,a:setattr(a,'bond',x),
            's_m_entry_name' : lambda code,x,a:setattr(a,'name',x), 
            's_ffio_ct_type' : lambda code,x,a:setattr(a,'ct_type',x),
            }
        fn_list = [ ]
        if len(self.load_properties):
            if self.load_properties == '*':
                for code in spec:
                    fn_list.append([ code, handler.get(code, set_molecule_property)])
            else:
                load_properties_list = self.load_properties.split(' ')
                for code in spec:
                    if code in load_properties_list:
                        fn_list.append([ code, handler.get(code, set_molecule_property)])
                    else:
                        fn_list.append([ code, handler.get(code,lambda code, x,a:0)])
        else:
            for code in spec:
                fn_list.append([ code, handler.get(code,lambda code, x,a:0)])
        bond_list = []
        data.reverse()
        mdl = Indexed()
        for code, fn in fn_list:
            apply(fn,(code, data.pop(), mdl))
        return mdl

    def close_brace(self,tok):
        (spec,data) = (self.spec,self.stack)
        (self.stack,self.spec) = self.meta_stack.pop()
        if len(self.stack):
            tag = self.stack.pop()
            if isinstance(tag, basestring):
                if tag[0:7]=='m_atom[' and tag[-1:]==']':
                    self.spec.append('chempy_atoms')
                    cur_atom_block = self.parse_atom_block(int(tag[7:-1]),
                                                           spec,
                                                           data,
                                                           self.full_atom_block)
                    self.stack.append( cur_atom_block )
                    if self.full_atom_block == None:
                        self.full_atom_block = cur_atom_block
                elif tag[0:7]=='m_bond[' and tag[-1:]==']':
                    self.spec.append('chempy_bonds')
                    expected_length = int(tag[7:-1])
                    cur_bond_block = self.parse_bond_block( int(tag[7:-1]),spec,
                                                            data,
                                                            self.full_bond_block)
                    if self.full_bond_block == None:
                        self.full_bond_block = cur_bond_block
                    bond_list = filter(lambda x:x.index[0]<x.index[1],cur_bond_block)
                    self.stack.append( bond_list )
                elif tag[0:6]=='f_m_ct':
                    self.spec.append('f_m_ct')
                    self.full_count = self.full_count + 1
                    rec = self.parse_ct_block(spec,data)
                    if rec:
                        rec.mae_record = 'f_m_ct'
                    self.stack.append( rec )
                elif tag[0:6]=='p_m_ct':
                    self.spec.append('p_m_ct')
                    rec = self.parse_ct_block(spec,data)
                    if rec:
                        rec.mae_record = 'p_m_ct'
                    self.stack.append( rec )
            else:
                self.stack.append(tag)
                
    def listFromStr(self,mae_st):
        self.stack = []
        self.spec = []
        self.meta_stack = []
        self.full_atom_block = None
        self.full_bond_block = None
        
        if mae_st.find('\r'):
            mae_st = _crlf_re.sub('\n',mae_st)
        mae_st = string.replace(mae_st,"<>","0")
        _dispatch = {
            '{' : self.open_brace,
            ':::' : self.triple_colon,
            '}' : self.close_brace,
            }

        for tok in quotedsplit(mae_st):
            _dispatch.get(tok, self.stack.append)(tok)

        return self.stack
        
    def fix_bond_reps(self,mdl):
        atom = mdl.atom
        for bond in mdl.bond:
            try:
                from_rep = bond.from_rep
            except AttributeError:
                from_rep = MMCT_BOND_WIRE

            try:
                to_rep = bond.to_rep
            except AttributeError:
                to_rep = MMCT_BOND_WIRE

            from_at = atom[bond.index[0]]
            to_at = atom[bond.index[1]]

            from_at.visible |= _bond_rep.get(from_rep, 0)
            to_at.visible   |= _bond_rep.get(to_rep, 0)

        for at in mdl.atom:
            if not getattr(at, 'visibility', 1):
                at.visible &= REPS_CARTOON

