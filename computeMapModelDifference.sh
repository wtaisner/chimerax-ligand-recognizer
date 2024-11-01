#!/bin/bash
# Author: John Heumann

if [ $# -ne 3 ]
then
  echo "Usage: $0 <dir> <map> <model>"
  echo "  Construct a cryoEM map - model difference map using Phenix"
  echo "    <dir>  : the directory containing <map> and <model>"
  echo "    <map>  : a cryoEM map of the <model> plus any ligands"
  echo "    <model>: a ligand-free .pdb or .cif model of the <map>"
  exit 1
fi

DIR=$1
MAP=$2
MODEL=$3
TYPE="${MODEL##*.}"

cd "$DIR" >/dev/null

if [ "$TYPE" == "pdb" ]
then
  RES=$(grep -i 'RESOLUTION.*ANGSTROMS' "$MODEL" | head -1 | \
    sed -e's/  */ /g' | cut -f 4 -d ' ')
elif [ "$TYPE" == "cif" ]
then
  RES=$(grep -i reconstruction.resolution "$MODEL" | grep -v method | \
    sed -e's/  */ /g' | cut -f 2 -d ' ')
else
  echo "ERROR: unrecognized model type!"
  exit 2
fi

echo "phenix.real_space_diff_map $MODEL $MAP resolution=$RES"
phenix.real_space_diff_map "$MODEL" "$MAP" "resolution=$RES"
rm -rf map_model_difference_2.ccp4
echo "Difference map is in $(pwd)/map_model_difference_1.ccp4"

exit 0
