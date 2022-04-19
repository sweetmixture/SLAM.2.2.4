#!/bin/bash

MAIN_INPUT=$1		# SLAM STANDARD OUTPUT FILE
DEBUG=$2
NO_RELAX=$2
CUR_PATH=$PWD		# 

if [ -z '$1' ]; then
    exit 1
fi

#################################################################################
# ENV VARIABLES
#################################################################################
SRC_PATH="/Users/woongkyujee/SLAM/SLAM.2.2.4/utils_slam/slam_dipole_analyser"
PYTHON_SRC_PATH=""
#SRC_PATH="/Users/woongkyujee/Desktop/Thesis_Sub_Chapters/SLAM_SRC/SLAM_TURING_01112021/SLAM.2.2.4/utils_slam/slam_dipole_analyser"
PYTHON_SRC_PATH="$SRC_PATH/calc_slam_dipole.py"
#################################################################################

# TMP FILES
CONFIG_TMP="config"	# CONFIG INFO
MO_TMP="mo"		# MO INFO
SPECIES_TMP="type"	# SPECIES INFO
GEN_INFO="info"		# GENERAL INFO INCLUDING POTENTIAL + CONFIG INFO

if [ -z $MAIN_INPUT ];	then
    echo "Error ... please check if the first argument is SLAM standard output!"
    echo "May ask help by feeding '-help'"
    exit 1
elif [ "$MAIN_INPUT" == "-help" ]; then


    echo "1> MainOutputs"
    echo "      dip_slam.out    (atomic dipole moment [Angs*e] unit, defined by the core positions of atoms + cluster total dipole)"
    echo "      dip_slam.vesta  (compatible with 'vesta' package for atomic dipole visualization as vectors)"
    echo ""
    echo " Before run this shell utility script, please do check the related environment variables"
    echo "      path/to/the/python/script ... see 'ENV VARIABLE' section in the script"
    echo "              the python code is for doing all the computation related with giving dipole momenet result of the cluster"
    echo "2> HowToRun"
    echo "Usage: bash get_slam_dipole.sh [standard slam output (energy relaxation used - recommended)] ... [option]"
    echo "Options and arguments"
    echo " -d           : keep all the temporal files"
    echo "              info    (slam input information - cluster configuration + used potential)"
    echo "              config  (xyz configuration including shell position if necessary)"
    echo "              mo      (molecular orbitals of lone pair ions/atoms for their ground states)"
    echo "              type    (object types including charge information)"
    echo " -nr          : no relaxation, ignore if the system is fully relaxed"
    exit 1
fi

if [ ! -f "$CUR_PATH"/"$MAIN_INPUT" ]; then
    echo "Error ... main input file does not exist in the specified path!"
    echo "May ask help by feeding '-help'"
    exit 1
fi

#################################################################################
# INTERNAL VARIABLES ... DO NOT MAKE ANY CHANGES
#################################################################################
POS_INT="Position Integral Reference"

#OPTI_SUCCESS_FLAG="Optimisation Meets Termination Condition"
OPTI_SUCCESS_FLAG="Termination Condition"
OPTI_FAIL_FLAG="Optimiser Failed"

CONFIG_SC_FLAG="CONFIGURATION_XYZ_SC_INFO"
LP_MO_FLAG="Lone Pair Molecular Orbital"

MM_INFO_FLAG="MM atoms/ions"
SP_INFO_FLAG="SP atoms/ions"

#################################################################################



# CHECK STANDARD OUTPUT RESULTS
SUCCESS_FLAG=$(grep "$OPTI_SUCCESS_FLAG" $MAIN_INPUT)		# IF SUCCESS $SUCCESS_FLAG IS NOT EMPTY
if [ ! -z "$SUCCESS_FLAG" ]; then

	IF_SUCCESS="1"

else
	IF_SUCCESS="0"

fi

POS_INTEGRAL=$( grep  "$POS_INT" $MAIN_INPUT | awk '{print $5}' )

# EXTRACT NUMBER OF SPECIES
READ_LN=$( grep -n "$CONFIG_SC_FLAG" $MAIN_INPUT | tail -1 | awk '{print $1}' | sed "s/://" )
NL=$(echo "$READ_LN + 1" | bc )

MM_CNT=$( sed -n "$NL"p $MAIN_INPUT | awk '{print $1}' )
QM_CNT=$( sed -n "$NL"p $MAIN_INPUT | awk '{print $2}' )
CNT=$( echo "$MM_CNT + $QM_CNT" | bc )												# GET TOTAL MM + TOTAL LP


# EXTRACT CONFIGURATION INFO ... WHICH USES THE LATEST (FROM THE BOTTOM OF $MAIN_INPUT)
READ_LN=$( grep -n "$CONFIG_SC_FLAG" $MAIN_INPUT | tail -1 | awk '{print $1}' | sed "s/://" )
L_STA=$( echo "$READ_LN + 2"  | bc )
L_END=$( echo "$L_STA + $CNT - 1" | bc )
echo $MM_CNT $QM_CNT > $CONFIG_TMP
sed -n "$L_STA","$L_END"p $MAIN_INPUT >> $CONFIG_TMP

# EXTRACT SLAM MOLECULAR ORBITAL INFO
READ_LN=$( grep -n "$LP_MO_FLAG" $MAIN_INPUT | tail -1 | awk '{print $1}' | sed "s/://" )
L_STA=$( echo "$READ_LN + 5"  | bc )
L_END=$( echo "$L_STA + $QM_CNT - 1" | bc )
echo $QM_CNT > $MO_TMP
sed -n "$L_STA","$L_END"p $MAIN_INPUT >> $MO_TMP

# EXTRACE SPECIES INFO
echo $MM_CNT $QM_CNT > $SPECIES_TMP
READ_LN=$( echo  "$MM_CNT + 4" | bc )
grep -A $READ_LN "$MM_INFO_FLAG" $MAIN_INPUT | tail -"$MM_CNT" | awk '{ printf "%2s%2s%12.4f\n", $1, $2, $6}' >> $SPECIES_TMP
READ_LN=$( echo  "$QM_CNT + 4" | bc )
grep -A $READ_LN "$SP_INFO_FLAG" $MAIN_INPUT | tail -"$QM_CNT" | awk '{ printf "%2s%20s\n", $1, $5}' >> $SPECIES_TMP

############################################################################### MAY NEED MODIFICATION DEPENDS ON THE SHELL TYPE
SHELL_TYPE=$(echo "$SHELL" | cut -c 6-)
if [ $SHELL_TYPE == "bash" ]; then
    sed -i "s/\//   /g" $SPECIES_TMP	# FOR GNU ... does not need 
elif [ $SHELL_TYPE == "zsh" ]; then
    sed -i "" "s/\//   /g" $SPECIES_TMP # OSX	  ... needs argument after -i 
fi
###############################################################################

# EXTRACT GENERAL INPUT INFO
L_STA=$( grep -n "####" $MAIN_INPUT | sed -n 3p | sed 's/[^0-9]*//g')
L_END=$( grep -n "####" $MAIN_INPUT | sed -n 4p | sed 's/[^0-9]*//g')

sed -n "$L_STA","$L_END"p $MAIN_INPUT > $GEN_INFO



# RUN PYTHON SCRIPT ------ INITIATE MAIN COMPUTATION 
python $PYTHON_SRC_PATH $POS_INTEGRAL $CONFIG_TMP $MO_TMP $SPECIES_TMP $MM_CNT $QM_CNT $IF_SUCCESS> dip_slam.out

if [ "$DEBUG" == "-d" ]; then
	echo "debugging mode requested ... leaving all temporal files"
else
	rm -rf $CONFIG_TMP $MO_TMP $SPECIES_TMP
fi
