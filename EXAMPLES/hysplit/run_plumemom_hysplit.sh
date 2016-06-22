#!/bin/sh

echo "### $0 ###"

#-------------------------------------------------------------
result=$(grep -i 'hysplit_dir' input_file.py | cut -c 15-)

temp="${result%\"}"
result="${temp#\"}"

MDL="$result/trunk"

#----------------------------------------------------------

python run_plumemom.py 

python create_hysplit_emitimes_control.py

python create_hysplit_setup_ascdata.py
 

${MDL}/exec/hycs_std  

echo "'TITLE&','### $0 ### &'" >LABELS.CFG
${MDL}/exec/parxplot -iPARDUMP -k1 -z80 -j${MDL}/graphics/arlmap
evince parxplot.ps
  
${MDL}/exec/concplot -icdump_Eya2010 -j${MDL}/graphics/arlmap -s0 -z80

evince concplot.ps

rm -f LABELS.CFG


  
