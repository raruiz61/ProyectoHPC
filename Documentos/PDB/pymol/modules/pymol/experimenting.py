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

if __name__=='pymol.experimenting':
    
    import selector
    from cmd import _cmd,lock,unlock,Shortcut,QuietException, \
          DEFAULT_ERROR, DEFAULT_SUCCESS, _raising, is_ok, is_error        
    import cmd
    import threading
    import pymol
    import string
    
    def get_bond_print(obj,max_bond,max_type,_self=cmd):
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.get_bond_print(_self._COb,str(obj),int(max_bond),int(max_type))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException                  
        return r

    def spheroid(object="",average=0,_self=cmd):  # EXPERIMENTAL
        '''
DESCRIPTION

    "spheroid" averages trajectory frames together to create
    an ellipsoid-like approximation of the actual anisotropic
    motion exhibited by the atom over a series of trajectory frames.

USAGE

    spheroid object,average

    average = number of states to average for each resulting spheroid state

    '''
        r = DEFAULT_ERROR
        try:
            print "Warning: 'spheroid' is experimental, incomplete, and unstable."
            _self.lock(_self)
            r = _cmd.spheroid(_self._COb,str(object),int(average))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException                  
        return r

    def mem(_self=cmd):
        '''
DESCRIPTION

    "mem" Dumps current memory state to standard output. This is a
    debugging feature, not an official part of the API.

    '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.mem(_self._COb)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException                  
        return r


    def check(selection=None, preserve=0):
        '''
DESCRIPTION

    "check" is unsupported command that may eventually have something
    to do with assigning forcefield parameters to a selection of
    atoms.
    
'''        
        # This function relies on code that is not currently part of PyMOL/ChemPy
        # NOTE: the realtime module relies on code that is not yet part of PyMOL/ChemPy
        from chempy.tinker import realtime
        if selection==None:
            arg = cmd.get_names("objects")
            arg = arg[0:1]
            if arg:
                if len(arg):
                    selection = arg
        if selection!=None:
            selection = selector.process(selection)
            realtime.assign("("+selection+")",int(preserve))
            realtime.setup("("+selection+")")

    def fast_minimize(*args, **kwargs):
        '''
DESCRIPTION

    "fast_minimize" is an unsupported nonfunctional command that may
    eventually have something to do with doing a quick clean up of the
    molecular structure.
    
'''        
        kwargs['_setup'] = 0
        return minimize(*args, **kwargs)

    def minimize(sele='', iter=500, grad=0.01, interval=50, _setup=1, _self=cmd):
        '''
DESCRIPTION

    "fast_minimize" is an unsupported nonfunctional command that may
    eventually have something to do with minimization.
    
'''        
        from chempy.tinker import realtime  

        if not sele:
            names = _self.get_names("objects")
            if not names:
                return
            sele = names[0]
        sele = '(' + sele + ')'

        if not int(_setup) or realtime.setup(sele):
            _self.async(realtime.mini, int(iter), float(grad), int(interval), sele)
        else:
            print " minimize: missing parameters, can't continue"


    def dump(fnam,obj,_self=cmd):
        '''
DESCRIPTION

    "dump" is an unsupported command which may have something to do
    with outputing isosurface meshes and surface objects to a file.

    '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.dump(_self._COb,str(fnam),obj)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException                  
        return r


    def dummy(*arg):
        return None

    def test(group=0,index=0,_self=cmd): # generic test routine for development
        '''
DESCRIPTION

    "dump" is an unsupported internal command.

    '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)   
            r=_cmd.test(_self._COb,int(group),int(index))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException                  
        return r

    def import_coords(coords,name,state,_self=cmd): # experimental
        r = DEFAULT_ERROR      
        try:
            _self.lock(_self)   
            r = _cmd.import_coords(_self._COb,str(name),int(state)-1,coords)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException                  
        return r

    def load_coords(model, oname, state=1): # UNSUPPORTED
        '''
        WARNING: buggy argument list, state get's decremented twice!
        '''
        return pymol.importing.load_coordset(model, oname, int(state)-1)

    class MatchMaker(object):
        '''
DESCRIPTION

    Experimental and subject to change!

    API only. Matches two atom selections and provides two matched
    subselections with equal atom count. May involve temporary objects
    or named selections which will be automatically deleted.

ARGUMENTS

    mobile = string: first atom selection

    target = string: second atom selection

    match = string: method how to match atoms
        * none: (dummy)
        * in: match atoms by "in" operator
        * like: match atoms by "like" operator
        * align: match atoms by cmd.align (without refinement)
        * super: match atoms by cmd.super (without refinement)
        * <name of alignment object>: use given alignment

RESULT

    Properties "mobile" and "target" hold the matched subselections as
    selection strings.
        '''
        def __init__(self, mobile, target, match, mobile_state=0, target_state=0):
            self.autodelete = True
            self.temporary = []

            self.mobile_state = mobile_state
            self.target_state = target_state

            if match == 'none':
                self.mobile = mobile
                self.target = target
            elif match in ['in', 'like']:
                self.mobile = '(%s) %s (%s) & state %i' % (mobile, match, target, target_state)
                self.target = '(%s) %s (%s) & state %i' % (target, match, mobile, mobile_state)
            elif match in ['align', 'super']:
                self.align(mobile, target, match, mobile_state, target_state)
            elif match in cmd.get_names('all') and \
                    cmd.get_type(match) in ('object:', 'object:alignment'):
                self.from_alignment(mobile, target, match)
            else:
                print ' Error: unkown match method', match
                raise pymol.CmdException

        def check(self):
            return cmd.count_atoms(self.mobile) == cmd.count_atoms(self.target)

        def align(self, mobile, target, match, mobile_state=0, target_state=0):
            '''
            Align mobile to target using the alignment method given by "match"
            '''
            aln_obj = cmd.get_unused_name('_')
            self.temporary.append(aln_obj)

            align = cmd.keyword[match][0]
            align(mobile, target, cycles=0, transform=0, object=aln_obj,
                    mobile_state=mobile_state, target_state=target_state)
            cmd.disable(aln_obj)

            self.from_alignment(mobile, target, aln_obj)

        def from_alignment(self, mobile, target, aln_obj):
            '''
            Use alignment given by "aln_obj" (name of alignment object)
            '''
            cmd.refresh()

            self.mobile = '(%s) and %s' % (mobile, aln_obj)
            self.target = '(%s) and %s' % (target, aln_obj)
            if self.check():
                return

            # difficult: if selections spans only part of the alignment or
            # if alignment object covers more than the two objects, then we
            # need to pick those columns that have no gap in any of the two
            # given selections

            mobileidx = set(cmd.index(mobile))
            targetidx = set(cmd.index(target))
            mobileidxsel = []
            targetidxsel = []

            for column in cmd.get_raw_alignment(aln_obj):
                mobiles = mobileidx.intersection(column)
                if len(mobiles) == 1:
                    targets = targetidx.intersection(column)
                    if len(targets) == 1:
                        mobileidxsel.extend(mobiles)
                        targetidxsel.extend(targets)

            self.mobile = cmd.get_unused_name('_mobile')
            self.target = cmd.get_unused_name('_target')
            self.temporary.append(self.mobile)
            self.temporary.append(self.target)

            mobile_objects = set(idx[0] for idx in mobileidxsel)
            target_objects = set(idx[0] for idx in targetidxsel)

            if len(mobile_objects) == len(target_objects) == 1:
                mobile_index_list = [idx[1] for idx in mobileidxsel]
                target_index_list = [idx[1] for idx in targetidxsel]
                cmd.select_list(self.mobile, mobile_objects.pop(), mobile_index_list, mode='index',
                        state=self.mobile_state)
                cmd.select_list(self.target, target_objects.pop(), target_index_list, mode='index',
                        state=self.target_state)
            else:
                cmd.select(self.mobile, ' '.join('%s`%d' % idx for idx in mobileidxsel),
                        state=self.mobile_state)
                cmd.select(self.target, ' '.join('%s`%d' % idx for idx in targetidxsel),
                        state=self.target_state)

        def __del__(self):
            if not self.autodelete:
                return
            for name in self.temporary:
                cmd.delete(name)

    def focal_blur(aperture=2.0, samples=10, ray=0, filename='', quiet=1, _self=cmd):
        '''
DESCRIPTION

    Creates fancy figures by introducing a focal blur to the image.
    The object at the origin will be in focus.

USAGE

    focal_blur [ aperture [, samples [, ray [, filename ]]]]

ARGUMENTS

    aperture = float: aperture angle in degrees {default: 2.0}

    samples = int: number of images for averaging {default: 10}

    ray = 0/1: {default: 0}

    filename = str: write image to file {default: temporary}

AUTHORS

    Jarl Underhaug, Jason Vertrees and Thomas Holder

EXAMPLES

    focal_blur 3.0, 50
        '''
        from tempfile import mkdtemp
        from shutil import rmtree
        from math import sin, cos, pi, sqrt

        try:
            import Image
        except ImportError:
            try:
                from PIL import Image
            except ImportError:
                raise pymol.CmdException("Python Image Library (PIL) not available")

        aperture, samples = float(aperture), int(samples)
        ray, quiet = int(ray), int(quiet)
        avg = None
        handle_list = []

        # make sure we have a still image
        _self.mstop()
        _self.unset('rock')
        _self.unset('text')

        # we need at least two samples
        if samples < 2:
            raise pymol.CmdException("need samples > 1")

        # Get the orientation of the camera and the light
        view = _self.get_view()
        lights_names = ['light', 'light2', 'light3', 'light4', 'light5',
                'light6', 'light7', 'light8', 'light9']
        lights = [_self.get_setting_tuple(lights_names[i])[1]
                for i in range(_self.get_setting_int('light_count') - 1)]

        # Create a temporary directory
        tmpdir = mkdtemp()

        try:
            # Rotate the camera and the light in order to create the blur
            for frame in range(samples):
                # Populate angles as Fermat's spiral
                theta = frame * pi * 110.0/144.0
                radius = 0.5 * aperture * sqrt(frame/float(samples-1))
                x = cos(theta) * radius
                y = sin(theta) * radius
                xr = x/180.0*pi
                yr = y/180.0*pi

                # Rotate the camera
                _self.turn('x', x)
                _self.turn('y', y)

                # Rotate the light
                for light, light_name in zip(lights, lights_names):
                    ly = light[1] * cos(xr) - light[2] * sin(xr)
                    lz = light[2] * cos(xr) + light[1] * sin(xr)
                    lx = light[0] * cos(yr) + lz * sin(yr)
                    lz = lz * cos(yr) - lx * sin(yr)
                    _self.set(light_name, [lx, ly, lz])

                # Save the image to temporary directory
                curFile = "%s/frame-%04d.png" % (tmpdir, frame)
                _self.png(curFile, ray=ray, quiet=1)

                if not quiet:
                    print " Created frame %i/%i (%0.0f%%)" % (frame + 1, samples,
                            100 * (frame + 1) / samples)

                # Create the average/blured image
                handle = open(curFile, "rb")
                handle_list.append(handle)
                img = Image.open(handle)
                avg = Image.blend(avg, img, 1.0 / (frame + 1)) if avg else img
                del img

                # Return the camera and the light to the original orientation
                _self.set_view(view)
                for light, light_name in zip(lights, lights_names):
                    _self.set(light_name, light)

            if not filename:
                filename = '%s/avg.png' % (tmpdir)

            # Load the blured image
            avg.save(filename)
            _self.load(filename)

        finally:
            del avg

            # needed for Windows
            for handle in handle_list:
                handle.close()

            # Delete the temporary files
            rmtree(tmpdir)

    def callout(name, label, pos='', screen='auto', state=-1, color='front',
            quiet=1, _self=cmd):
        '''
DESCRIPTION

    Create a new screen-stabilized callout object.

ARGUMENTS

    name = str: object name

    label = str: label text

    pos = str or list: anchor in model space as 3-float coord list or atom
    selection. If empty, don't draw an arrow. {default: }

    screen = str or list: position on screen as 2-float list between [-1,-1]
    (lower left) and [1,1] (upper right) or "auto" for smart placement.
    {default: auto}
        '''
        state, quiet = int(state), int(quiet)
        sele = ''

        if state == -1:
            state = _self.get_setting_int('state')

        if isinstance(pos, basestring):
            if ',' in pos:
                pos = _self.safe_list_eval(pos)
            else:
                sele, pos = pos, None

        if isinstance(screen, basestring):
            if screen == 'auto':
                screen = (0.0, 0.0)
            else:
                screen = _self.safe_list_eval(screen)

        if len(screen) == 2:
            screen += (0.0,)
        elif len(screen) != 3:
            raise pymol.CmdException('invalid screen argument')

        _self.pseudoatom(name, sele, 'CALL', 'CAL',
                str(_self.count_atoms('?' + name) + 1),
                'C', 'CALL', 'PS', 0.5, color=color, label=label, pos=pos,
                state=state, quiet=quiet)

        expr = ('('
                's["label_relative_mode"],'
                's["label_connector"],'
                's["label_connector_color"],'
                's["label_color"],'
                ') = (1, %s, "default", "default",)' % (
                    bool(pos or sele),
                ))

        expr_state = ('('
                's["label_screen_point"],'
                ') = (%s,)' % (
                    screen,
                ))

        _self.alter("last " + name, expr)
        _self.alter_state(state, "last " + name, expr_state)

    def desaturate(selection="all", a=0.5, quiet=1, _self=cmd):
        '''
DESCRIPTION

    Desaturate the colors in the given selection.

ARGUMENTS

    selection = str: atom selection {default: all}

    a = float [0..1]: desaturation factor {default: 0.5}
        '''
        a, b = float(a), 1.0 - float(a)
        def desat(color):
            rgb = [int(255 * (i * b + a)) for i in _self.get_color_tuple(color)]
            return 0x40000000 + (rgb[0] << 16) + (rgb[1] << 8) + rgb[2]
        _self.alter(selection, "color = desat(color)", space={'desat': desat})
        _self.recolor()
