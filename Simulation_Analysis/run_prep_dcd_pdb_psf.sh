#!/bin/bash

cp ../*.dcd .
cp ../all_stress*.txt .

cp ../spider_HA1B1A1B1x1440_silkworm_A1B1x0_water_only.psf ref.psf
cp ../spider_HA1B1A1B1x1440_silkworm_A1B1x0.data ref.data

module load vmd
vmd -dispdev text -e make_refpdb.tcl
