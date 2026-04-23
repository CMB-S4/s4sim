#/bin/bash

fname_targets=max-depth-targets.txt
if [[ ! -e $fname_targets ]]; then
    python make_targets.py
fi

fname_schedule=schedule.sat.txt
fname_log=${fname_schedule/.txt/.log}

echo "Writing schedule to $fname_schedule"

toast_ground_schedule \
    @schedule.sat.par \
    @max-depth-targets.txt \
    --out $fname_schedule \
    >& $fname_log

fname_analysis=${fname_schedule/.txt/.analysis}
echo "Writing analysis to $fname_analysis"
toast_analyze_schedule $fname_schedule >& $fname_analysis

fname_plot=${fname_schedule/.txt/.png}
echo "Projecting schedule to $fname_plot"
toast_project_schedule $fname_schedule --fov 35 --nside 64 --bg smooth_npipe6v20_353_map.fits --bg-percentile 75 --bg-pol --no-step-plots --azstep 5
mv hits_tot.png $fname_plot
