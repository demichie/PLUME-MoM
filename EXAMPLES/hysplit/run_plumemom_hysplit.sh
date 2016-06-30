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

DUMP_ACC="cdumpcum_$result"

#----------------------------------------------------------

python run_plumemom.py 

python create_hysplit_emitimes_control.py

python create_hysplit_setup_ascdata.py
 

${MDL}/exec/hycs_std  

echo "'TITLE&','### $0 ### &'" >LABELS.CFG
${MDL}/exec/parxplot -iPARDUMP -k1 -z80 -j${MDL}/graphics/arlmap
evince parxplot.ps
  
${MDL}/exec/concplot -i$DUMP -j${MDL}/graphics/arlmap -s0 -z80 -d1 -ukg -oconcplot.ps

evince concplot.ps

${MDL}/exec/concacc -i$DUMP -o$DUMP_ACC

${MDL}/exec/concplot -i$DUMP_ACC -j${MDL}/graphics/arlmap -s0 -t0 -z80 -d1 -ukg -oconcplot_cum.ps

evince concplot_cum.ps


rm -f LABELS.CFG


grep -A100000 POINTS input_file.py|grep -v "POINTS" > con2stn.inp

sed -i 's/P/0/' con2stn.inp
sed 's/=/ /' con2stn.inp > con2stn.tmp
sed 's/\[/ /' con2stn.tmp > con2stn.inp
sed 's/,/ /' con2stn.inp > con2stn.tmp
sed 's/\]//' con2stn.tmp > con2stn.inp

${MDL}/exec/con2stn -i$DUMP_ACC -scon2stn.inp -d0 -p0 -xi -z1 -ocon2stn.txt

python extract_samples.py

rm con2stn.inp
rm con2stn.tmp

evince gsd.pdf

