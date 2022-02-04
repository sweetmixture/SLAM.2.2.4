#!/bin/bash

# dev - local use only
ICC="icc -O3 -std=c11 -I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include"

$ICC -c CH_Integrals.c
$ICC -c CDH_Integrals.c
$ICC -c CDDH_Integrals.c
$ICC -c CDDDH_Integrals.c

$ICC -c BM_SH_Integrals.c
$ICC -c BM_SDH_Integrals.c
$ICC -c BM_SDDH_Integrals.c
$ICC -c BM_SDDDH_Integrals.c

