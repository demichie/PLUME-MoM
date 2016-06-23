#!/bin/sh

echo "### $0 ###"

#-------------------------------------------------------------
result=$(grep -i 'hysplit_dir' input_file.py | cut -c 15-)

temp="${result%\"}"
result="${temp#\"}"
temp="${result%\'}"
result="${temp#\'}"

MDL="$result"

result=$(grep -i 'runname' input_file.py | cut -c 11-)

temp="${result%\"}"
result="${temp#\"}"

temp="${result%\'}"
result="${temp#\'}"

DUMP="cdump_$result"

#----------------------------------------------------------

python run_plumemom.py 

python create_hysplit_emitimes_control.py

python create_hysplit_setup_ascdata.py
 

${MDL}/exec/hycs_std  

echo "'TITLE&','### $0 ### &'" >LABELS.CFG
${MDL}/exec/parxplot -iPARDUMP -k1 -z80 -j${MDL}/graphics/arlmap
evince parxplot.ps
  
${MDL}/exec/concplot -i$DUMP -j${MDL}/graphics/arlmap -s0 -z80 -d1

evince concplot.ps

rm -f LABELS.CFG


  