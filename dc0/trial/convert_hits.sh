#!/bin/bash

mkdir -p hitmaps

toast_healpix_convert hitmap_outputs_chlat/mapmaker_hits.h5 --outdir hitmaps --outmap chlat_hits.fits
toast_healpix_convert hitmap_outputs_splat/mapmaker_hits.h5 --outdir hitmaps --outmap splat_hits.fits
toast_healpix_convert hitmap_outputs_spsat_2/mapmaker_hits.h5 --outdir hitmaps --outmap spsat_hits.fits
