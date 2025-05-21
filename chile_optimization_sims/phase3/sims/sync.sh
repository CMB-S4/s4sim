#!/bin/bash

rsync -avrP scaled_outputs/* /global/cfs/cdirs/cmbs4/chile_optimization/simulations/phase3/noise_depth/

rsync -avrP scaled_daily_outputs/* /global/cfs/cdirs/cmbs4/chile_optimization/simulations/phase3/noise_depth/daily/

chgrp cmbs4 -R /global/cfs/cdirs/cmbs4/chile_optimization/simulations/phase3/noise_depth
chmod g+rX-w -R /global/cfs/cdirs/cmbs4/chile_optimization/simulations/phase3/noise_depth
