#!/bin/bash
if [ $# -ne 3 ]
then
  echo "Usage: $0 <dir> <map> <model>"
  exit 1
fi

DIR=$1
MAP=$2
MODEL=$3
TYPE="${MODEL##*.}"

cd "$DIR" >/dev/null

if [ "$TYPE" == "pdb" ]
then
  RES=$(grep -i 'RESOLUTION.*ANGSTROMS' *.pdb | head -1 | \
    sed -e's/  */ /g' | cut -f 4 -d ' ')
elif [ "$TYPE" == "cif" ]
then
  RES=$(grep -i reconstruction.resolution *stripped.cif | grep -v method | \
    sed -e's/  */ /g' | cut -f 2 -d ' ')
else
  echo "ERROR: unrecognized model type!"
  exit 2
fi

echo "phenix.real_space_diff_map $MODEL $MAP resolution=$RES"
phenix.real_space_diff_map "$MODEL" "$MAP" "resolution=$RES"

exit 0
