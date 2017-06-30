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

DUMP_SUM="cdumpsum_$result"

PDUMP="pdump_$result"


#----------------------------------------------------------

python run_plumemom.py 

python create_hysplit_emitimes_control.py

python create_hysplit_setup_ascdata.py
 

${MDL}/exec/hycs_std  

echo "-------------- start postprocessing ---------------"

python create_maptext.py 

echo "'TITLE&','### $0 ### &'" >LABELS.CFG
${MDL}/exec/parxplot -i$PDUMP -k1 -z20 -j${MDL}/graphics/arlmap

${MDL}/exec/par2asc -i$PDUMP -oPARDUMP.txt 
    
${MDL}/exec/concplot -i$DUMP -j${MDL}/graphics/arlmap -s0 -z20 -d1 -ukg -oconcplot.ps

${MDL}/exec/concacc -i$DUMP -o$DUMP_ACC

${MDL}/exec/concsum -i$DUMP_ACC -o$DUMP_SUM

${MDL}/exec/concplot -i$DUMP_ACC -j${MDL}/graphics/arlmap -s0 -t0 -z20 -d1 -ukg -oconcplot_cum.ps

rm -f LABELS.CFG

# echo "-------------- extract loading and GSD at locs ---------------"

grep -A100000 POINTS input_file.py|grep -v "POINTS" > con2stn.tmp0

sed 's/P/0/' con2stn.tmp0 > con2stn.tmp1
sed 's/=/ /' con2stn.tmp1 > con2stn.tmp2
sed 's/\[/ /' con2stn.tmp2 > con2stn.tmp3
sed 's/,/ /' con2stn.tmp3 > con2stn.tmp4
sed 's/\]//' con2stn.tmp4 > con2stn.tmp5
sed '/^$/d' con2stn.tmp5 > con2stn.inp # elimina le ultime righe bianche dal file con2stn0.inp

${MDL}/exec/con2stn -i$DUMP_ACC -scon2stn.inp -d0 -p0 -xi -z1 -r0 -ocon2stn.txt

python extract_samples.py

rm con2stn.tmp*

echo "-------------- convert ps to pdf ---------------"

filelist=`(find . -name \*.ps)`

for i in $filelist; do

        ps2pdf $i

        rm $i
done

#python check_deposit.py




