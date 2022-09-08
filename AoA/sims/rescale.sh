#!/bin/bash

rm -rf scaled_outputs

python3 rescale.lat.py &
python3 rescale.sat.py &

wait

chmod g+rX -R scaled_outputs &
chgrp cmbs4 -R scaled_outputs &

python3 plot_depth.py &

wait
