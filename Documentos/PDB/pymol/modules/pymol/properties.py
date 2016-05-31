#A* -------------------------------------------------------------------
#B* This file contains source code for the PyMOL computer program
#C* Copyright (c) Schrodinger, LLC. 
#D* -------------------------------------------------------------------
#E* It is unlawful to modify or remove this copyright notice.
#F* -------------------------------------------------------------------
#G* Please see the accompanying LICENSE file for further information. 
#H* -------------------------------------------------------------------
#I* Additional authors of this source file include:
#-* 
#-* 
#-*
#Z* -------------------------------------------------------------------

PROPERTY_AUTO = -1
PROPERTY_BOOLEAN = 1
PROPERTY_INT = 2
PROPERTY_FLOAT = 3
PROPERTY_FLOAT3 = 4
PROPERTY_COLOR = 5
PROPERTY_STRING = 6

DEFAULT_SUCCESS = None

if __name__=='pymol.properties':
    import cmd
    from cmd import _cmd, DEFAULT_ERROR, DEFAULT_SUCCESS, QuietException
    from querying import get_color_index_from_string_or_list
    self = cmd

    def get_property(propname, name, state=0, quiet=1, _self=cmd):
        '''
DESCRIPTION

    Get an object-level property

ARGUMENTS

    propname = string: Name of the property

    name = string: Name of a single object

    state = int: Object state, 0 for all states, -1 for current state
    {default: 0}

    '''
        state, quiet = int(state)-1, int(quiet)
        r = DEFAULT_ERROR
        proptype = None
        propval = None
        if not len(str(propname)):
            return None
        try:
            _self.lock(_self)   
            r = _self._cmd.get_property(_self._COb, propname, name, state, quiet)
            if isinstance(r, list):
                proptype = r[0]
                propval = r[1]
        except:
            return None
        finally:
            _self.unlock(None,_self)
        if not quiet:
            if propval and proptype == PROPERTY_COLOR:
                try:
                    propval = dict((k,v) for (v,k) in _self.get_color_indices())[propval]
                except:
                    pass
            if propval:
                print " get_property: '%s' in object '%s' : %s"%(propname, name, propval)
            else:
                print " get_property: '%s' in object '%s' not found"%(propname, name)

        return propval

    def get_property_list(object, state=0, quiet=1, _self=cmd):
        '''
DESCRIPTION

    Get all properties for an object (for a particular state) as a list

ARGUMENTS

    object = string: Name of a single object

    state = int: Object state, 0 for all states, -1 for current state
    {default: 0}
    '''
        state, quiet = int(state)-1, int(quiet)
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)   
            r = _self._cmd.get_property(_self._COb, '', object, state, quiet)
            if not quiet:
                print " get_property_list: %s : %s"%(object, r)
        finally:
            _self.unlock(None,_self)
        return r
        
    def set_property(name, value, object='*', state=0, proptype=PROPERTY_AUTO, quiet=1, _self=cmd):
        '''
DESCRIPTION

    Set an object-level property

ARGUMENTS

    name = string: Name of the property

    value = string (or arbitrary type!?): Value to be set

    object = string: Space separated list of objects or * for all objects
    {default: *}

    proptype = int: The type of the property (PROPERTY_AUTO detects the type by default)
                    if the proptype is PROPERTY_COLOR (5), then the value can
                    be a color string, a list of 3 floats, or an integer/color index 

    state = int: Object state, 0 for all states, -1 for current state
    {default: 0}
    '''
        state, quiet = int(state)-1, int(quiet)
        proptype = int(proptype)
        if proptype == PROPERTY_COLOR:
            value = get_color_index_from_string_or_list(value, _self=_self)
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)   
            r = _self._cmd.set_property(_self._COb, name, value, object, proptype, state, quiet)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise QuietException
        return r

    def set_atom_property(name, value, selection='all', state=0, proptype=PROPERTY_AUTO, quiet=1, _self=cmd):
        '''
DESCRIPTION

    Set an atom-level property

ARGUMENTS

    name = string: Name of the property

    value = string (or arbitrary type!?): Value to be set

    selection = string: a selection-expression
    {default: all}

    proptype = int: The type of the property (PROPERTY_AUTO detects the type by default)
                    if the proptype is PROPERTY_COLOR (5), then the value can
                    be a color string, a list of 3 floats, or an integer/color index 

    state = int: Object state, 0 for all states, -1 for current state
    {default: 0}
    '''
        state, quiet = int(state)-1, int(quiet)
        proptype = int(proptype)
        if proptype == PROPERTY_COLOR:
            value = get_color_index_from_string_or_list(value, _self=_self)
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)   
            r = _self._cmd.set_atom_property(_self._COb, name, value, selection, proptype, state, quiet)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise QuietException
        return r
