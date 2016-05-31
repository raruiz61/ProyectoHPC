#!/bin/bash

# create wrapper, if filename provided
if [[ -n "$1" ]]; then
    wrapper="$1"
    echo '======================================================================'
    echo "Creating wrapper script: '$wrapper'"
    echo '======================================================================'
    echo ''
    echo '#!/bin/sh' > $wrapper
    echo '#' >> $wrapper
    echo '# PyMOL startup script' >> $wrapper
    echo '#' >> $wrapper
    echo '# Set PYMOL_PATH to point to this directory' >> $wrapper
    echo '#' >> $wrapper
    echo '# ================================================' >> $wrapper
    echo "PYMOL_PATH=`pwd`" >> $wrapper
    echo '# ================================================' >> $wrapper
    echo '#' >> $wrapper
    echo 'export PYMOL_PATH' >> $wrapper
    echo 'PYTHONHOME=$PYMOL_PATH/ext' >> $wrapper
    echo 'PYTHONPATH=$PYTHONHOME/lib/python2.7:$PYTHONPATH' >> $wrapper
    echo 'PYTHONPATH=$PYTHONHOME/lib/python2.7/lib-tk:$PYTHONPATH' >> $wrapper
    echo 'LD_LIBRARY_PATH=$PYMOL_PATH/ext/lib:$LD_LIBRARY_PATH' >> $wrapper
    echo 'LANG=C' >> $wrapper
    echo 'export LD_LIBRARY_PATH' >> $wrapper
    echo 'export PYTHONHOME' >> $wrapper
    echo 'export PYTHONPATH' >> $wrapper
    echo 'export LANG' >> $wrapper
    echo 'exec $PYMOL_PATH/pymol.exe "$@"' >> $wrapper
    chmod 755 $wrapper
else
    echo '======================================================================'
    echo 'Running setup.sh no longer creates a wrapper script unless you provide'
    echo 'an output filename (./setup.sh /path/to/pymolwrapper). We recommend to'
    echo 'just use the existing "pymol" launch script from this directory. You'
    echo 'may want to create a symlink in $PATH, for example:'
    echo 'ln -s $PWD/pymol /usr/local/bin/'
    echo '======================================================================'
    echo ''
    wrapper=pymol
fi

# explicit path for execution
case "$wrapper" in
    /*|./*|../*)
        ;;
    *)
        wrapper="./$wrapper"
        ;;
esac

echo '======================================================================'
echo 'Checking if any library from "ext/libextra" is required on this system'
echo '======================================================================'
echo ''

# remove extra libs from ext/lib
(cd ext/libextra && names=`ls` && cd ../lib && rm -f $names)

# add extra libs to ext/lib if there seems to be a problem like
# missing library or incompatible library
while true; do
    name=(`$wrapper -cq 2>&1 | grep 'lib[^ /]*\.so' | sed 's!.*\(lib[^ /]*\.so[.0-9]*\).*!\1!'`)
    test -z "$name" && break
    if [ -e "ext/lib/$name" -o ! -e "ext/libextra/$name" ]; then
        echo "warning: problem with $name"
        break
    fi
    (cd ext/lib && ln -s ../libextra/$name .)
done

