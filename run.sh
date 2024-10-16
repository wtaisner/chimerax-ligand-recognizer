# command without displaying a graphics window (--nogui)
# and exits immediately after installation (exit true).
# The initial --exit flag guarantees that ChimeraX will exit even
# if installation fails for some reason.

# set ChimeraX alias
alias chimerax=/usr/bin/chimerax-daily  # change this to your ChimeraX path
# clean previous builds
chimerax --nogui --cmd 'devel clean . exit true'
rm -rf build dist ChimeraX_HelloWorld.egg-info
# build the package
chimerax --nogui --cmd 'devel build . exit true'
# install the package
chimerax --nogui --cmd 'devel install . exit true'
# run the command
chimerax --nogui --cmd 'blobus validatus 8SOR/map_model_difference_1.ccp4'