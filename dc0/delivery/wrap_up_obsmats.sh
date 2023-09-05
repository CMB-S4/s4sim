#!/bin/bash

indir="/global/cfs/cdirs/cmbs4/dc/dc0/staging/obsmat/obsmats"
outdir="/global/cfs/cdirs/cmbs4/dc/dc0/mission/spsat"

for band in 025 040 085 095 145 155 230 280; do
    case $band in
        025) band_in=f030;;
        230) band_in=f0220;;
        *) band_in=f${band};;
    esac
    fname_in=${indir}/full/${band_in}/obsmat_${band_in}.npz
    fname_out=${outdir}/split01/${band}/dc0_chlat_t01.01_${band}_mat03.npz
    rsync -avrP $fname_in $fname_out
done
