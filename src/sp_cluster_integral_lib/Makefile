# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# (c) Author:    Woongkyu Jee, woong.jee.16@ucl.ac.uk, wldndrb1@gmail.com
# Created:   02.06.2019 ~
# 	
# University College London, Department of Chemistry, KLMC
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# ENVS - MAY REQUIRE CHANGE BELOW
#
# IF USING UCL / ARCHER HPC / MakeSure load appropriate mods e.g., intel / intel-mpi / gnu-mpi / gsl 
# MCC=mpiicc -O3 -std=c11 -mcmodel=large
# CC=icc -O3 -std=c11 -mcmodel=large
# IC=icc -O3 -std=c11 -mcmodel=large
#
IC=icc -O3 -std=c11 -I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include
# for local osx dev .. intel compiler use
# CC=gcc -O3 -std=c11

# Compile Integral Library
#
# * CXH:
#       basic Coulomb interaction + Analytical Derivatives
#       e.g., CH: Scalar (Energy) / CDH: 1st-derivative / CDDH: 2nd-derivative ...
#
# * BM_SXH: BM (Born-Mayer type exponetial repulsive)
#       similar convention used
#       e.g., *_SH: Scalar (Energy) / *_SDH: 1st-derivative ...
#
# 'sp_cluster_integrator*': integral solver, as function of model parameters
#       Making pre-calculated look-up-table
#

dummy.o: dummy.c CH_Integrals.o CDH_Integrals.o CDDH_Integrals.o  CDDDH_Integrals.o \
        BM_SH_Integrals.o BM_SDH_Integrals.o BM_SDDH_Integrals.o BM_SDDDH_Integrals.o
	$(IC) -c dummy.c

CH_Integrals.o: CH_Integrals.c CH_Integrals.h 
	$(IC) -o $@ -c CH_Integrals.c

BM_SH_Integrals.o: BM_SH_Integrals.c BM_SH_Integrals.h 
	$(IC) -o $@ -c BM_SH_Integrals.c

CDH_Integrals.o: CDH_Integrals.c CDH_Integrals.h 
	$(IC) -o $@ -c CDH_Integrals.c

BM_SDH_Integrals.o: BM_SDH_Integrals.c BM_SDH_Integrals.h 
	$(IC) -o $@ -c BM_SDH_Integrals.c

CDDH_Integrals.o: CDDH_Integrals.c CDDH_Integrals.h 
	$(IC) -o $@ -c CDDH_Integrals.c

BM_SDDH_Integrals.o: BM_SDDH_Integrals.c BM_SDDH_Integrals.h 
	$(IC) -o $@ -c BM_SDDH_Integrals.c

CDDDH_Integrals.o: CDDDH_Integrals.c CDDDH_Integrals.h 
	$(IC) -o $@ -c CDDDH_Integrals.c

#BM_SDDDH_Integrals.o: BM_SDDDH_Integrals.c BM_SDDDH_Integrals.h 
#	$(IC) -o $@ -c BM_SDDDH_Integrals.c

# DEV - YUKAWA SHORT_RANGE POT TEST

#SH_Integrals.o: SH_Integrals.c SH_Integrals.h
#	$(IC) -o $@ -c SH_Integrals.c
#SDH_Integrals.o: SDH_Integrals.c SDH_Integrals.h
#	$(IC) -o $@ -c SDH_Integrals.c

clean:
	rm *.o
