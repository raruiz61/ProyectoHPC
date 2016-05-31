'''
Various morphing workflows

(c) 2013 Schrodinger, Inc.
(c) 2011 Thomas Holder
(c) 2009 DeLano Scientific LLC.
'''

import cmd
from pymol import CmdException

def morph(name, sele1, sele2=None, state1=-1, state2=-1, refinement=3,
        steps=30, method='rigimol', match='align', quiet=1, _self=cmd):
    '''
DESCRIPTION

    Creates an interpolated trajectory between two or multiple conformations.
    If the two input objects are not the same, match them based on sequence
    alignment.

    This command supports two methods: rigimol and linear. RigiMOL is an
    incentive feature and only available to official PyMOL sponsors. Linear
    morphing is quick and robust but likely to produce distorted intermediates.

ARGUMENTS

    name = string: name of object to create

    sele1 = string: atom selection of first conformation

    sele2 = string: atom selection of second conformation {default: <sele1>}

    state1 = int: sele1 state {default: 1}. If state1=0 and sele1 has N
    states, create N morphings between all consecutive states and back from
    state N to 1 (so the morph will have N*steps states). If state2=0, create
    N-1 morphings and stop at last state.

    state2 = int: sele2 state {default: 2 if sele1=sele2, else 1}

    refinement = int: number of sculpting refinement cycles to clean
    distorted intermediates {default: 3}

    steps = int: number of states for sele2 object {default: 30}

    method = string: rigimol or linear {default: rigimol}

EXAMPLE

    fetch 1akeA 4akeA, async=0
    align 1akeA, 4akeA
    morph mout, 1akeA, 4akeA
    '''
    state1, state2 = int(state1), int(state2)
    refinement , steps, quiet = int(refinement), int(steps), int(quiet)
    method = method.lower()

    n_sele1 = _self.count_states(sele1)
    if not sele2 and n_sele1 < 2:
        raise CmdException('sele1 has only one state and no sele2 given')

    try:
        from epymol.rigimol import run_morph
    except ImportError:
        if method == 'rigimol':
            print ' Warning: rigimol not available (incentive-only feature),' \
                    ' fallback to linear morphing'
            method = 'linear'
        if refinement > 0:
            print ' Warning: refine not available (incentive-only feature)'
            refinement = 0

    if not name:
        name = _self.get_unused_name('morph')

    if state1 == 0:
        # multi-state morph

        if sele2 and sele2 != sele1:
            raise CmdException('sele1 and sele2 must be equal if state1=0')

        tmp = _self.get_unused_name('_tmpmorph')
        source1 = _self.get_unused_name('_tmpsource1')
        source2 = _self.get_unused_name('_tmpsource2') if cmd.count_discrete(sele1) else None

        if state2 < 1: state2 = n_sele1 - state2

        for i in range(1, state2):
            j = i % n_sele1 + 1
            _self.select(source1, sele1, 0, state=i)
            if source2:
                _self.select(source2, sele1, 0, state=j)
            morph(tmp, source1, source2, i, j, refinement, steps, method, 'none', quiet)
            _self.create(name, tmp, 0, -1)
            _self.delete(tmp)

        _self.delete(source1)
        _self.delete(source2)
        return

    if not sele2 or sele2 == sele1:
        if state1 < 1: state1 = 1
        if state2 < 1: state2 = 2
        source = sele1
    else:
        from .querying import get_selection_state
        from .experimenting import MatchMaker

        if state1 < 1: state1 = get_selection_state(sele1)
        if state2 < 1: state2 = get_selection_state(sele2)

        mm = MatchMaker(sele1, sele2, match, state1, state2)

        # 2-state morph-in object
        source = _self.get_unused_name('_sele1')
        _self.create(source, mm.mobile, state1, 1)
        _self.create(source, source, 1, 2)
        _self.update(source, mm.target, 2, state2, matchmaker=0)

        state1, state2 = 1, 2
        mm.temporary.append(source)

    if method != 'rigimol':
        run_morph = run_morph_linear

    run_morph(source, name, refinement, state1, state2, steps, 0, quiet)

def run_morph_linear(source, target, refinement, first, last, steps, n_state,
        quiet, _self=cmd):
    # This must take exactly the same arguments like rigimol.run_morph

    steps_h = steps / 2

    for x in range(steps_h):
        _self.create(target, source, first, x + 1)
    for x in range(steps_h, steps):
        _self.create(target, source, last, x + 1)

    _self.smooth(target, 1, steps, 1, steps)

    if refinement > 0:
        try:
            from epymol import rigimol
            rigimol.refine(refinement, target, quiet=quiet)
        except ImportError:
            print " Warning: refine not available (incentive-only feature)"

# vi: ts=4:sw=4:smarttab:expandtab
