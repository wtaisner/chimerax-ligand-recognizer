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

# test the command: first open the structure and the map, then run the command:
chimerax --nogui --cmd 'open 8SOR/8sor.cif; open 8SOR/map_model_difference_1.ccp4; blobus validatus /R:101'
