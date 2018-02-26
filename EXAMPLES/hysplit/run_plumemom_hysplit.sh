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

DUMP_PART="cdump_part_$result"

DUMP_ACC_PART="cdumpcum_part_$result"

DUMP_SUM_PART="cdumpsum_part_$result"

PDUMP_PART="pdump_part_$result"

DUMP_GAS="cdump_gas_$result"

DUMP_ACC_GAS="cdumpcum_gas_$result"

DUMP_SUM_GAS="cdumpsum_gas_$result"

PDUMP_GAS="pdump_gas_$result"

#----------------------------------------------------------

python run_plumemom.py 

python create_hysplit_emittimes_control.py

python create_hysplit_setup_ascdata.py
 
echo "-------------- particles dispersion simulation ---------------"

${MDL}/exec/hycs_std part  

echo "-------------- gas dispersion simulation ---------------"

${MDL}/exec/hycs_std gas  

echo "-------------- start postprocessing ---------------"

#python create_maptext.py !check!

echo "'PARTICLES &','### $0 ### &'" >LABELS.CFG

${MDL}/exec/parxplot -i$PDUMP_PART -k1 -z20 -j${MDL}/graphics/arlmap

${MDL}/exec/par2asc -i$PDUMP_PART -oPARDUMP_PART.txt 
    
${MDL}/exec/concplot -i$DUMP_PART -j${MDL}/graphics/arlmap -s0 -z20 -d1 -ukg -oconcplot_part.ps

${MDL}/exec/concacc -i$DUMP_PART -o$DUMP_ACC_PART

${MDL}/exec/concsum -i$DUMP_ACC_PART -o$DUMP_SUM_PART

${MDL}/exec/concplot -i$DUMP_ACC_PART -j${MDL}/graphics/arlmap -s0 -t0 -z20 -d1 -ukg -oconcplot_cum_part.ps

echo "'GAS &','### $0 ### &'" >LABELS.CFG

${MDL}/exec/parxplot -i$PDUMP_GAS -k1 -z20 -j${MDL}/graphics/arlmap

${MDL}/exec/par2asc -i$PDUMP_GAS -oPARDUMP_GAS.txt 
    
${MDL}/exec/concplot -i$DUMP_GAS -j${MDL}/graphics/arlmap -s0 -z20 -d1 -ukg -oconcplot_gas.ps

${MDL}/exec/concacc -i$DUMP_GAS -o$DUMP_ACC_GAS

${MDL}/exec/concsum -i$DUMP_ACC_GAS -o$DUMP_SUM_GAS

${MDL}/exec/concplot -i$DUMP_ACC_GAS -j${MDL}/graphics/arlmap -s0 -t0 -z20 -d1 -ukg -oconcplot_cum_gas.ps


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




