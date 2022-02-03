#!/bin/bash

# Converting 'xyz' format file into SLAM input 'geo.txt'
#
TAR_XYZ=$1		# As 1st input argument, target xyz file name
SP_CATION_NAME=$2	# As 2nd input argument, user must specify sp-ion name
IF_SHELL=$3		# As 3rd input argument, if MM ions are with shell (optional)
# ARG_1 : TARGET_XYZ FILE  / ARG_2 : ATOM TYPE OF LONE PAIR / ARG_3 : IF THERE IS SHELL FOR MM
#################################################################################
if [ -z $TAR_XYZ ]; then							# CHECK IF TARGET XYZ FILE
	printf "%s\n" "Error, there is no xyz file found for input ..."
	exit 1
else
	if [ ! -f $TAR_XYZ ]; then						# DOUBLE CHECK IF FILE DOES EXIST
		printf "%s\n" "Error, the input xyz file does not exist ..."
		exit 1
	fi
fi

if [ -z $SP_CATION_NAME ]; then							# CHECK LP CATION NAME
	printf "%s\n" "Error, the first argument place is empty ..."
	exit 1
fi

#if [ -z $IF_SHELL ]; then							# CHECK IF SHELL IS USED
#	printf "%s\n" "MM atoms/ions are treated as rigid ion  model"
#else 
#	printf "%s\n" "MM atoms/ions are treated as core-shell model"
#fi										# END OF READ INPUT PARAMETERS
#################################################################################

ATOM_NUMBER=$( sed -n 1p $TAR_XYZ )						# READ NUMBER OF ATOMS
MODEL_EN=$( sed -n 2p $TAR_XYZ)							# READ COMMENT LINE ... USUALLY ENERGY (EV) BUT SOMETIMES NOTHING WRITTEN

declare -a MM	# MM ION ARR
declare -a SP	# SP ION ARR
MM_CNT=0	# MM ION CNT
SP_CNT=0	# SP ION CNT

for (( i=0; i<$ATOM_NUMBER; i++)); do
	
	line_flag=$( echo "$i+3" | bc )				# SET LINE NUMBER
	read_line=$(sed -n "$line_flag"p $TAR_XYZ)		# READ LINE
	spliter=( $read_line )					# SPLIT LINE

	if [ $SP_CATION_NAME == ${spliter[0]} ]; then		# CHECK IF SP CATION
		for (( j=0; j<4; j++ )); do
			SP[$SP_CNT*4+$j]=${spliter[$j]}
		done
		SP_CNT=$( echo "$SP_CNT + 1" | bc )
	else							# CHECK IF MM ION
		for (( j=0; j<4; j++ )); do
			MM[$MM_CNT*4+$j]=${spliter[$j]}
		done
		MM_CNT=$( echo "$MM_CNT + 1" | bc )
	fi
done								# READ CONFIG DONE

if [ -f tmp ]; then
	rm tmp
fi
touch tmp		# 'tmp' temporal workspace


# WRITE SLAM 'geo.txt' OUTPUT WILL BE SAVED IN 
OUT="geo.txt.out"
#
for (( i=0; i<$MM_CNT; i++ )); do
	if [ -z $IF_SHELL ]; then 				# $IF_SHELL is empty then use RIM 
		printf "%3.2s%3.2s%12.6lf%12.6lf%12.6lf\n" ${MM[$i*4+0]} "c" ${MM[$i*4+1]} ${MM[$i*4+2]} ${MM[$i*4+3]}	>> tmp
	else							# $IF_SHELL is not empty then use core-shell
		printf "%3.2s%3.2s%12.6lf%12.6lf%12.6lf\n" ${MM[$i*4+0]} "c" ${MM[$i*4+1]} ${MM[$i*4+2]} ${MM[$i*4+3]}	>> tmp
		printf "%3.2s%3.2s%12.6lf%12.6lf%12.6lf\n" ${MM[$i*4+0]} "s" ${MM[$i*4+1]} ${MM[$i*4+2]} ${MM[$i*4+3]}	>> tmp
	fi	
done		# MM WRT DONE

for (( i=0; i<$SP_CNT; i++ )); do
	printf "%3.2s%15.6lf%12.6lf%12.6lf\n" ${SP[$i*4+0]} ${SP[$i*4+1]} ${SP[$i*4+2]} ${SP[$i*4+3]} >> tmp
done		# SP WRT DONE

lm_dummy=$( wc "tmp" | awk '{print $1}' )		# CNT TOTAL SPECIES
SLAM_MM_CNT=$( echo "$lm_dummy - $SP_CNT" | bc )	# CNT MM SPECIES

echo "#src_path:"$PWD"/"$TAR_XYZ > $OUT
echo $SLAM_MM_CNT $SP_CNT >> $OUT
cat tmp >> $OUT
rm tmp


