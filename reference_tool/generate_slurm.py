#!/usr/bin/env python

""" This script generates SLURM scripts to simulate input maps for the simulation tool.
"""


import os
import sys


"""
  --bands BANDS         Comma-separated list of bands: ULFL1 (20 GHz, LAT),
                        LFL1 (27 GHz LAT), LFL2 (39 GHz, LAT), LFS1 (30 GHz,
                        SAT), LFS2 (40 GHz, SAT), MFL1 (93 GHz, LAT), MFL2
                        (145 GHz, LAT), MFLS1 (85 GHz, SAT), MFLS2 (145.1 GHz,
                        SAT), MFHS1 (95 GHz, SAT), MFHS2 (155.1 GHz, SAT),
                        HFL1(225 GHz, LAT), HFL2 (278 GHz, LAT), HFS1 (220
                        GHz, SAT), HFS2 (270 GHz, SAT).Length of list must
                        equal --tubes
  --tubes TUBES         Comma-separated list of optics tubes: LT0 (HFL), LT1
                        (HFL), LT2 (HFL), LT3 (HFL), LT4 (HFL), LT5 (MFL), LT6
                        (MFL), LT7 (MFL), LT8 (MFL), LT9 (MFL), LT10 (MFL),
                        LT11 (MFL), LT12 (MFL), LT13 (MFL), LT14 (MFL), LT15
                        (MFL), LT16 (MFL), LT17 (LFL), LT18 (LFL), LT19 (HFL),
                        LT20 (HFL), LT21 (HFL), LT22 (HFL), LT23 (HFL), LT24
                        (MFL), LT25 (MFL), LT26 (MFL), LT27 (MFL), LT28 (MFL),
                        LT29 (MFL), LT30 (MFL), LT31 (MFL), LT32 (MFL), LT33
                        (MFL), LT34 (MFL), LT35 (MFL), LT36 (LFL), LT37 (LFL),
                        LT38 (HFL), LT39 (HFL), LT40 (HFL), LT41 (HFL), LT42
                        (MFL), LT43 (MFL), LT44 (MFL), LT45 (MFL), LT46 (MFL),
                        LT47 (MFL), LT48 (MFL), LT49 (MFL), LT50 (MFL), LT51
                        (MFL), LT52 (MFL), LT53 (MFL), LT54 (LFL), LT55 (LFL),
                        LT56 (ULFL), ST0 (MFLS), ST1 (MFLS), ST2 (MFLS), ST3
                        (MFLS), ST4 (MFLS), ST5 (MFLS), ST6 (MFHS), ST7
                        (MFHS), ST8 (MFHS), ST9 (MFHS), ST10 (MFHS), ST11
                        (MFHS), ST12 (HFS), ST13 (HFS), ST14 (HFS), ST15 (HFS),
                        ST16 (LFS), ST17 (LFS).Length of list must equal
                        --bands
"""

input_map_dir = "/global/cscratch1/sd/zonca/cmbs4/map_based_simulations/202006_foregrounds_extragalactic_cmb_tophat"

flavors = (
    "noise",
    "atmosphere",
    "cmb-unlensed", # cmb_unlensed_solardipole_nest
    "cmb-lensing", # cmb_lensing_signal
    "cmb-tensors", # cmb_tensor_nest
    "foreground", # combined_foregrounds
)

telescopes = {
    "LAT": {
        "LT0": ["HFL1", "HFL2"],  #  225 & 278 GHz
        "LT5": ["MFL1", "MFL2"],  #   93 & 145 GHz
        "LT17": ["LFL1", "LFL2"],  #  27 &  39 GHz
        "LT56": ["ULFL1"],  #               20 GHz
    },
    "SAT": {
        "ST0": ["MFLS1", "MFLS2"],  #  85 & 145.1 GHz - SAT0 - FOV 14.5 deg
        "ST6": ["MFHS1", "MFHS2"],  #  95 & 155.1 GHz - SAT2 - FOV 14.5 deg
        "ST12": ["HFS1", "HFS2"],  #  220 & 270 GHz - SAT4 - FOV 17.5 deg
        "ST16": ["LFS1", "LFS2"],  #   30 &  40 GHz - SAT5 - FOV 17.5 deg
    },
}

for telescope, tubes in telescopes.items():
    if telescope == "LAT":
        nside = 4096
        fsample = 200
        hwprpm = None
        scan_rate = 1
        scan_accel = 1
        poly_order = 15
        ground_order = 25
        fpradius = 4.3
        nnode = 64
        nthread = 16
        nnode_group = 1
        madampars = {
            "madam-concatenate-messages": None,
            "madam-allreduce": None,
            "madam-precond-width": 30,
            "nside-submap": 16,
            "madam-baseline-length": 1,
            "madam-noisefilter": None,
        }
        cosecant_scan = True
        thinfp = 8
    elif telescope == "SAT":
        nside = 512
        fsample = 20
        hwprpm = 120
        scan_rate = 1
        scan_accel = 1
        poly_order = 5
        ground_order = 10
        nnode = 16
        fpradius = 20.5
        nthread = 4
        nnode_group = 1
        madampars = {
            "no-madam-allreduce": None,
            "madam-precond-width": 30,
            "nside-submap": 16,
            "madam-baseline-length": 1,
            "madam-noisefilter": None,
        }
        cosecant_scan = False
        thinfp = 4
    else:
        raise RuntimeError("Unknown telescope: {}".format(telescope))

    # For now, we disable destriping and only output filtered maps
    madampars = {
        # Comment out skip-madam for a pure binned map (2 separate Madam calls)
        "skip-madam": None,
        "no-destripe": None,
    }

    for site in "chile", "pole":
        if site == "chile":
            weather = "weather_Atacama.fits"
        elif site == "pole":
            weather = "weather_South_Pole.fits"
            hwprpm = None
            cosecant_scan = False
        else:
            raise RuntimeError("Unknown site: {}".format(site))

        schedule = "scan_strategy/{}_{}/schedules/{}_schedule_{}.txt".format(
            site, telescope.lower(), site, telescope.lower()
        )

        for tube, bands in tubes.items():
            for band in bands:
                thinfp_temp = thinfp
                if band.startswith("ULF") or band.startswith("LF"):
                    thinfp_temp = 1
                elif band.startswith("HFS"):
                    thinfp_temp = 8
                hardware = "hardware_{}_{}.toml.gz".format(telescope, band[:-1])
                for flavor in flavors:
                    rootname = "{}_{}_{}_{}".format(site, flavor, telescope, band)
                    os.makedirs("slurm", exist_ok=True)
                    os.makedirs("logs", exist_ok=True)
                    
                    params = {
                        "bands": band,
                        "tubes": tube,
                        "sample-rate": fsample,
                        "scan-rate": scan_rate,
                        "scan-accel": scan_accel,
                        "nside": nside,
                        "schedule": schedule,
                        "weather": weather,
                        "site": site,
                        "madam-concatenate-messages": None,
                        "no-madam-allreduce": None,
                        "focalplane-radius": fpradius,
                        "madam-prefix": rootname,
                        "thinfp": thinfp_temp,
                        "hardware": hardware,
                    }
                    if cosecant_scan:
                        params["scan-cosecant-modulate"] = None
                    if poly_order is not None:
                        params["polyfilter"] = None
                        params["poly-order"] = poly_order
                    if ground_order is not None:
                        params["groundfilter"] = None
                        params["ground-order"] = ground_order

                    if flavor == "noise":
                        params["simulate-noise"] = None
                        params["hits"] = None
                        params["wcov"] = None
                        params["wcov-inv"] = None
                        params["MC-count"] = 8
                        walltime = "02:00:00"
                    elif flavor == "atmosphere":
                        params["simulate-atmosphere"] = None
                        params["no-hits"] = None
                        params["no-wcov"] = None
                        params["no-wcov-inv"] = None
                        params["MC-count"] = 8
                        walltime = "24:00:00"
                    elif flavor in [
                        "cmb-unlensed",
                        "cmb-lensing",
                        "cmb-tensors",
                        "foreground",
                    ]:
                        # params["input-map"] = input_map
                        params["no-hits"] = None
                        params["no-wcov"] = None
                        params["no-wcov-inv"] = None
                        params["skip-madam"] = None
                        signal_name = {
                            "cmb-unlensed" : "cmb_unlensed_solardipole_nest",
                            "cmb-lensing" : "cmb_lensing_signal",
                            "cmb-tensors" : "cmb_tensor_nest",
                            "foreground" : "combined_foregrounds",
                        }[flavor]
                        num = "0000"
                        params["input-map"] = os.path.join(
                            input_map_dir,
                            str(nside),
                            signal_name,
                            num,
                            "cmbs4_{}_uKCMB_{}-{}_nside{}_{}.fits".format(
                                signal_name, telescope, band, nside, num
                            )
                        )
                        walltime = "00:30:00"
                    else:
                        raise RuntimeError(
                            "Unknown simulation flavor: '{}'".format(flavor)
                        )

                    params.update(madampars)
                        
                    fname_slurm = os.path.join("slurm", "{}.slrm".format(rootname))
                    with open(fname_slurm, "w") as slurm:
                        for line in [
                            "#!/bin/bash",
                            "#SBATCH --partition=regular",
                            "#SBATCH --time={}".format(walltime),
                            "#SBATCH --nodes={}".format(nnode),
                            "#SBATCH --job-name={}".format(rootname),
                            "#SBATCH --licenses=SCRATCH",
                            "#SBATCH --constraint=knl",
                            "#SBATCH --core-spec=4",
                            "#SBATCH --account=mp107",
                            "\nulimit -c unlimited",
                            "export MALLOC_MMAP_THRESHOLD_=131072",
                            'export PYTHONSTARTUP=""',
                            "export PYTHONNOUSERSITE=1",
                            "export HOME=$SCRATCH",
                            "export OMP_NUM_THREADS={}".format(nthread),
                            "export OMP_PLACES=threads",
                            "export OMP_PROC_BIND=spread",
                            "\nlet nnode={}".format(nnode),
                            "let ntask_node=64/$OMP_NUM_THREADS",
                            "let ntask=$nnode*$ntask_node",
                            "let ncore=4*$OMP_NUM_THREADS",
                            "# Make sure nnode is divisible by nnode_group",
                            "let nnode_group={}".format(nnode_group),
                            "let groupsize=nnode_group*ntask_node",
                            '\necho "Running with"',
                            'echo "            nnode = ${nnode}"',
                            'echo "  OMP_NUM_THREADS = ${OMP_NUM_THREADS}"',
                            'echo "       ntask_node = ${ntask_node}"',
                            'echo "            ntask = ${ntask}"',
                            'echo "      nnode_group = ${nnode_group}"',
                            'echo "        groupsize = ${groupsize}"',
                            '\nexport PYTHONSTARTUP=""',
                            "export PYTHONNOUSERSITE=1",
                            "\nlogfile=logs/{}.log\n".format(rootname),
                            "if [[ ! -e $logfile ]]; then",
                            '    echo "Writing $logfile"',
                            "    srun -n $ntask -c $ncore --cpu_bind=cores \\",
                            "        toast_s4_sim.py @general.par \\",
                            "        --group-size $groupsize \\",
                        ]:
                            slurm.write(line + "\n")

                        for key in sorted(params.keys()):
                            if params[key] is None:
                                slurm.write("        --{} \\\n".format(key))
                            else:
                                slurm.write(
                                    "        --{} {} \\\n".format(key, params[key])
                                )

                        for line in [
                            "    >& ${logfile}",
                            "else",
                            '    echo "$logfile exists"',
                            "fi",
                        ]:
                            slurm.write(line + "\n")
