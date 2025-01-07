# CMB-S4 transient simulation template

Make sure to (edit and) run
- `./make_schedule.sh`
- `./make_focalplane.sh`
- `python make_catalog.py`

before submitting run_sim.slurm

## make_schedule.py

The script generates a TOAST observing schedule for one calendar day. It should
be obvious how to change that. The single schedule is then broken up into
single observation schedules that can be mapped separately.

## make_focalplane.sh

Uses `s4sim` to create TOAST focalplane files for the CMB-S4 CHLAT instrument

## make_catalog.py

Example python script to generate three catalog entries.  Just a placeholder.

## run_sim.slurm

SLURM script that loops over the split observing schedules and maps the first
one that does not have a matching log file created yet. To complete the
simulation, you need to submit the same job repeatedly.  It is fine to queue
them all together.
