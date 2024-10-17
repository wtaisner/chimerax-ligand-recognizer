# command without displaying a graphics window (--nogui)
# and exits immediately after installation (exit true).
# The initial --exit flag guarantees that ChimeraX will exit even
# if installation fails for some reason.

# set ChimeraX alias -- change this to your ChimeraX path
#alias chimerax=/usr/bin/chimerax
alias chimerax=/usr/bin/chimerax-daily
# clean previous builds
chimerax --nogui --cmd 'devel clean . exit true'
rm -rf build dist ChimeraX_HelloWorld.egg-info
# build the package
chimerax --nogui --cmd 'devel build . exit true'
# install the package
chimerax --nogui --cmd 'devel install . exit true'

# test the command: first open the structure and the map, then run the command and exit
chimerax --nogui --exit --cmd 'open 8SOR/8sor.cif; open 8SOR/map_model_difference_1.ccp4; blobus validatus /A:1401'
chimerax --nogui --exit --cmd  'open 8SOR/8sor.cif; open 8SOR/map_model_difference_1.ccp4; blob validate /A:1401'
