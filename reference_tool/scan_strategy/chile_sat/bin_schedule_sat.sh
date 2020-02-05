##!/bin/bash

OMP_NUM_THREADS=2 mpirun -np 4 \
	       python3 $PREFIX/bin/toast_s4_sim.py \
	       ../bin_schedule.par \
	       @bin_schedule_sat.par \
	       >& bin_schedule_sat.log
