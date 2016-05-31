# PyMOL epymol.mae module

from dotmae import MAEReader, MAEReaderException
from pymol import cmd
from chempy import cpv
import math
import re
import itertools

def read_maestr(maestr,name,state=0,finish=1,discrete=-1,
                quiet=1,zoom=-1,multiplex=-1,mimic=1,
                object_props=None, atom_props=None, _self=cmd):
    cmd = _self
    mimic = int(mimic)
    if not quiet:
        print " Load: Processing Schrodinger .MAE file..."

    if object_props is None:
        object_props = _self.get('load_object_props_default')
    if atom_props is None:
        atom_props = _self.get('load_atom_props_default')

    mr = MAEReader(mimic, object_props, atom_props)
    list = mr.listFromStr(maestr)

    if list and getattr(list[0], 'ct_type', None) == 'full_system':
        if not quiet:
            print " Desmond CMS file detected, only loading first CT"
        list = list[:1]

    len_list = len(list);
    if len_list<2: # MAE contains no more than one molecule
        if discrete<0:
            discrete = 0
    else: # MAE contains multiple molecules
        if multiplex > 0:
            if len_list > mr.full_count:
                multiplex = 2 # multiplex at the record level (break out conformers, etc.)
        if multiplex < 0:
            if len_list == mr.full_count: # file only contains full records
                multiplex = 1 # then multiplex at the object level
        elif multiplex == 0:
            if mr.full_count>1: # multiple full records present
                if discrete<0:
                    discrete = 1 # use discrete object to prevent collisions
    if discrete < 0:
        discrete = 0
    preserve_chempy_ids=cmd.get('preserve_chempy_ids')
    cmd.set('preserve_chempy_ids')
    
    # nest p_m_ct behind f_m_ct
    nested_list = []
    for mdl in list:
        if (mdl.mae_record == 'f_m_ct') or not nested_list:
            nested_list.append( [mdl] )
        else:
            nested_list[-1].append(mdl)

    if multiplex < 0:
        multiplex = 1 if (len(nested_list) > 1) else 0
    if multiplex:
        cmd.group(name)

    unnamed_cnt = itertools.count(1)

    object_name_dict = {}
    nested_list.reverse()
    while nested_list:
        mdl_list = nested_list.pop()
        len_mdl_list = len(mdl_list)
        if len_mdl_list:
            mdl = mdl_list[0]
            mdl_name = name
            if multiplex:
                mdl_name += "." + (
                        getattr(mdl, 'title', None) or
                        getattr(mdl, 'name', None) or
                        str(unnamed_cnt.next()))
            mdl_name = cmd.get_legal_name(mdl_name)
            if object_name_dict.has_key(mdl_name) and multiplex == 1:
                mdl_name = cmd.get_unused_name(mdl_name + "_")
            if multiplex > 0:
                name_dict = {}
                if multiplex > 1: # break coordinate sets into separate objects
                    cnt = 0
                    for mdl in mdl_list:
                        cnt = cnt + 1
                        mr.fix_bond_reps(mdl)
                        state_name = mdl_name
                        if len_mdl_list>1:
                            state_name += "_" + str(cnt)
                        _self.load_model(mdl,state_name,state,zoom=zoom,
                                         discrete=discrete,quiet=quiet, _self=_self)
                        object_name_dict[state_name] = 1
                else: # only break objects
                    for mdl in mdl_list:
                        mr.fix_bond_reps(mdl)
                        _self.load_model(mdl,mdl_name,state,zoom=zoom,
                                         discrete=discrete,quiet=quiet, _self=_self)
                        object_name_dict[mdl_name] = 1                    
            else: # multiplex = 0 , so stuff everything into a single object 
                for mdl in mdl_list:
                    mr.fix_bond_reps(mdl)
                    _self.load_model(mdl,mdl_name,state,zoom=zoom,
                                     discrete=discrete,quiet=quiet, _self=_self)
                    object_name_dict[mdl_name] = 1
    cmd.set('preserve_chempy_ids',preserve_chempy_ids)

    if mimic:
        for name in object_name_dict.keys():
            # maestro-like settings
            _self.flag("ignore",name+" and not polymer","set")
            _self.set("stick_radius",0.18,name)
            _self.set("line_width",1.2,name)
            _self.set("valence_size",0.075,name)
            _self.set("valence", 1, name)

        _self.set("light",[0.4,-0.4,-1.0])
        _self.set("shininess",100)
        _self.set("light2",[-0.5,-0.7,-0.2])
        _self.set("light_count",3)

        _self.set("label_position",[1,0,2])
        _self.set("label_font_id",7)
#    _self.set("cartoon_highlight_color","gray50")    
    
def set_view(view, center, trans, volume, persp):
    '''
    view is their 4x4

    center is the centor of rotation in tranformed coordinates
    
    '''
    view = list(view)
    
    # print view[0:4]
    # print view[4:8]
    # print view[8:12]
    # print view[12:16]
    # print "\nmaestro_center: %8.3f %8.3f %8.3f"%tuple(center)
    # print "maestro_transl: %8.3f %8.3f %8.3f"%tuple(trans)

    mat = [ view[0:3], view[4:7], view[8:11] ]
    tmp = cpv.sub(center,view[12:15])
    pymol_center = cpv.transform(mat,tmp)
    
    # print "volumeX: %8.3f %8.3f"%(tuple(volume[0:2]))
    # print "volumeY: %8.3f %8.3f"%(tuple(volume[2:4]))
    # print "volumeZ: %8.3f %8.3f"%(tuple(volume[4:6]))

    # print "pymol_center: %8.3f %8.3f %8.3f"%tuple(pymol_center)

    # print persp
    if persp<0.0:
        pymol_ortho = 1.0
    else:
        pymol_ortho = 0.0

    if pymol_ortho==0.0:
        fov = 1 + 100*(persp - 1.0)/3.0
        cmd.set('field_of_view',fov)

    fov = (math.pi * float(cmd.get("field_of_view")) / 180.0)

    # unclear why we have to multiple by 1.17 to get the correct look
        
    camera_dist = 1.17 * (abs(volume[2] - volume[3])/2.0)/math.atan(fov/2.0)

    vol_center = [ (volume[0]+volume[1])/2.0,
                   (volume[2]+volume[3])/2.0,
                   (volume[4]+volume[5])/2.0 ]
    
    pymol_objective = [ trans[0] + center[0] - vol_center[0],
                        trans[1] + center[1] - vol_center[1],
                        -camera_dist ]
    
    # print "pymol_objective: %8.3f %8.3f %8.3f"%tuple(pymol_objective)

    pymol_front = center[2] + trans[2] - volume[4] - pymol_objective[2] 
    pymol_back  = center[2] + trans[2] - volume[5] - pymol_objective[2] 

    pymol_view = cmd.get_view()
    cur_view = ( view[0:3] + view[4:7] + view[8:11] +
                 pymol_objective + pymol_center +
                 [ pymol_front, pymol_back, pymol_ortho ] )

    # print cur_view

    cmd.set_view(tuple(cur_view))
    
