#!/bin/sh
cd ${0%/*} || exit 1    # Run from this directory

# Parse arguments for library compilation
. $WM_PROJECT_DIR/wmake/scripts/AllwmakeParseArguments

#-

echo -n "Patching the hexRef8.H file... " 
patchDir=polyTopoChange/polyTopoChange/hexRef8
cp $patchDir/hexRef8.H $FOAM_SRC/dynamicMesh/$patchDir/hexRef8.H
cd $WM_PROJECT_DIR
./Allwmake -update -j 8 
echo "done".
cd -                

wmake all;

#------------------------------------------------------------------------------
