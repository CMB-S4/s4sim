#!/usr/bin/env python

""" This script generates SLURM scripts to simulate input maps for the simulation tool.
"""


from glob import glob
import os
import sys

"""

  --bands BANDS         Comma-separated list of bands: ULFPL1 (20 GHz, Pole
                        LAT), LFL1 (27 GHz LAT), LFL2 (39 GHz, LAT), LFPL1 (27
                        GHz Pole LAT), LFPL2 (39 GHz, Pole LAT), LFS1 (30 GHz,
                        SAT), LFS2 (40 GHz, SAT), MFL1 (93 GHz, LAT), MFL2
                        (145 GHz, LAT), MFPL1 (93 GHz, Pole LAT), MFPL2 (145
                        GHz, Pole LAT), MFLS1 (85 GHz, SAT), MFLS2 (145.1 GHz,
                        SAT), MFHS1 (95 GHz, SAT), MFHS2 (155.1 GHz, SAT),
                        HFL1(225 GHz, LAT), HFL2 (278 GHz, LAT), HFPL1 (225
                        GHz, Pole LAT), HFPL2 (278 GHz, Pole LAT), HFS1 (220
                        GHz, SAT), HFS2 (270 GHz, SAT).
  --tubes TUBES         Comma-separated list of optics tubes: LAT0-LFL : LT63,
                        LT66, LT67, LT70, LT75, LT78, LT79, LT82. LAT0-MFL :
                        LT19..LT22, LT25..LT31, LT34..LT62, LT64, LT65, LT68,
                        LT69, LT71..LT74, LT76, LT77, LT80, LT81, LT83, LT84.
                        LAT0-HFL : LT0..LT18, LT23, LT24, LT32, LT33. LAT1-LFL
                        : LT148, LT151, LT152, LT155, LT160, LT163, LT164,
                        LT167. LAT1-MFL : LT104..LT107, LT110..LT116,
                        LT119..LT150, LT153, LT154, LT156..LT159, LT161,
                        LT162, LT165, LT166, LT168, LT169. LAT1-HFL :
                        LT85..LT103, LT108, LT109, LT117, LT118. LAT2-ULFPL :
                        LT178, LT182, LT184, LT188. LAT2-LFPL : LT232, LT234,
                        LT236, LT239, LT242, LT245, LT248, LT251, LT254.
                        LAT2-MFPL : LT170..LT176, LT180, LT186, LT189..LT206,
                        LT208, LT210, LT212, LT214, LT216, LT218, LT220,
                        LT222, LT224, LT226, LT228, LT230, LT231, LT233,
                        LT235, LT237, LT238, LT240, LT241, LT243, LT244,
                        LT246, LT247, LT249, LT250, LT252, LT253. LAT2-HFPL :
                        LT177, LT179, LT181, LT183, LT185, LT187, LT207,
                        LT209, LT211, LT213, LT215, LT217, LT219, LT221,
                        LT223, LT225, LT227, LT229. SAT0-MFLS : ST0. SAT0-MFHS
                        : ST1. SAT0-HFS : ST2. SAT1-MFLS : ST3. SAT1-MFHS :
                        ST4. SAT1-HFS : ST5. SAT2-MFLS : ST6. SAT2-MFHS : ST7.
                        SAT2-HFS : ST8. SAT3-MFLS : ST9. SAT3-MFHS : ST10.
                        SAT3-HFS : ST11. SAT4-LFS : ST14. SAT4-MFLS : ST12.
                        SAT4-MFHS : ST13. SAT5-LFS : ST17. SAT5-MFLS : ST15.
                        SAT6-MFHS : ST16.

"""

input_map_dir = "/global/cscratch1/sd/zonca/cmbs4/map_based_simulations/202006_foregrounds_extragalactic_cmb_tophat"

flavors = (
    "noise",
    "atmosphere",
    #"cmb-unlensed", # cmb_unlensed_solardipole_nest
    #"cmb-lensing", # cmb_lensing_signal
    #"cmb-tensors", # cmb_tensor_nest
    #"foreground", # combined_foregrounds
)

telescopes = {
    "LAT": {
        "LT0": ["HFL1", "HFL2"], #           225 & 278 GHz
        "LT5": ["MFL1", "MFL2"],  #           93 & 145 GHz
        "LT17": ["LFL1", "LFL2"],  #          27 &  39 GHz
        "LT38": ["HFPL1", "HFPL2"],  #  Pole 225 & 278 GHz
        "LT42": ["MFPL1", "MFPL2"],  #  Pole  93 & 145 GHz
        "LT54": ["LFPL1", "LFPL2"],  #  Pole  27 &  39 GHz
        "LT56": ["ULFPL1"],  #          Pole 20 GHz
    },
    "SAT": {
        "ST0": ["MFLS1", "MFLS2"],  #  85 & 145.1 GHz - SAT0 - FOV 14.5 deg
        "ST6": ["MFHS1", "MFHS2"],  #  95 & 155.1 GHz - SAT2 - FOV 14.5 deg
        "ST12": ["HFS1", "HFS2"],  #  220 & 270 GHz - SAT4 - FOV 17.5 deg
        "ST16": ["LFS1", "LFS2"],  #   30 &  40 GHz - SAT5 - FOV 17.5 deg
    },
}


def get_n_obs(site, telescope):
    fnames = glob(f"scan_strategy/{site}_{telescope}/split_schedules/*txt".lower())
    return len(fnames)


for telescope in "SAT", "LAT":
    # For now, we disable destriping and only output filtered maps
    madampars = {
        # Comment out skip-madam for a pure binned map (2 separate Madam calls)
        "skip-madam": None,
        "no-destripe": None,
    }
    atm_cache = "atm_cache"
    for site in "pole", "chile":
        schedule = "scan_strategy/{}_{}/schedules/{}_schedule_{}.txt".format(
            site, telescope.lower(), site, telescope.lower()
        )
        hwprpm = None
        cosecant_scan = False
        poly_order_2d = None
        if site == "chile":
            weather = "weather_Atacama.fits"
        elif site == "pole":
            weather = "weather_South_Pole.fits"
        if telescope == "LAT":
            nside = 4096
            fsample = 200
            fpradius = 5.0
            nthread = 4
            """
            madampars = {
                "madam-concatenate-messages": None,
                "madam-allreduce": None,
                "madam-precond-width": 30,
                "nside-submap": 16,
                "madam-baseline-length": 1,
                "madam-noisefilter": None,
            }
            """
            thinfp = {
                "ULFPL" : 4,
                "LFPL" : 4,
                "MFPL" : 16,
                "HFPL" : 16,
                "LFL" : 4,
                "MFL" : 16,
                "HFL" : 16,
            }
            if site == "pole":
                nnode_group = 2
                nnode = 181 * nnode_group
                scan_rate = 1.0
                scan_accel = 1.0
                telescope_name = "LAT2"
                pixel_types = {"ULFPL" : None, "LFPL" : None, "MFPL" : None, "HFPL" : None}
                poly_order = 10
                ground_order = 100
            elif site == "chile":
                nnode_group = 2
                nnode = 500
                scan_rate = 0.5
                scan_accel = 3.0
                cosecant_scan = True
                telescope_name = "LAT0"
                pixel_types = {"LFL" : None, "MFL" : None, "HFL" : None}
                poly_order = 25
                ground_order = 15
                #poly_order_2d = 1
                atm_cache = "atm_cache_test"
            else:
                raise RuntimeError(f"Unknown site: {site}")
        elif telescope == "SAT":
            nnode = 16
            fpradius = 16
            nthread = 4
            nnode_group = 1
            """
            madampars = {
                "no-madam-allreduce": None,
                "madam-precond-width": 30,
                "nside-submap": 16,
                "madam-baseline-length": 1,
                "madam-noisefilter": None,
            }
            """
            thinfp = {
                "LFS" : 4,
                "MFLS" : 4,
                "MFHS" : 4,
                "HFS" : 8,
            }
            nside = 512
            fsample = 20
            poly_order = 3
            telescope_name = None
            pixel_types = {
                "LFS" : "--tubes ST14",
                "MFLS" : "--tubes ST0",
                "MFHS" : "--tubes ST1",
                "HFS" : "--tubes ST2",
            }
            if site == "pole":
                scan_rate = 1.5
                scan_accel = 0.97
                ground_order = 50
            elif site == "chile":
                scan_rate = 1.0
                scan_accel = 1.0
                hwprpm = 120
                ground_order = 10
            else:
                raise RuntimeError(f"Unknown site: {site}")
        else:
            raise RuntimeError(f"Unknown telescope: {telescope}")

        
        for pixel_type, tubes in pixel_types.items():
            for i_band in [1, 2]:
                band = f"{pixel_type}{i_band}"
                if band == "ULFPL2":
                    continue
                thinfp_temp = thinfp[pixel_type]
                #hardware = f"hardware_{telescope}_{pixel_type}.toml.gz"
                hardware = f"hardware_{telescope}_{pixel_type}.pkl"
                for flavor in flavors:
                    parfiles = "@general.par"
                    rootname = "{}_{}_{}_{}".format(site, flavor, telescope, band.replace("P", ""))
                    os.makedirs("slurm", exist_ok=True)
                    os.makedirs("logs", exist_ok=True)

                    params = {
                        "bands": band,
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
                        "out" : "out",
                    }

                    # HACK for fixed CHLAT noise level
                    #if telescope == "LAT" and site == "chile":
                    #    params["no-elevation-noise"] = None

                    if telescope_name is not None:
                        params["telescope"] = telescope_name
                    if tubes is not None:
                        params["tubes"] = tubes
                    if hwprpm is not None:
                        params["hwp-rpm"] = hwprpm
                    if cosecant_scan:
                        params["scan-cosecant-modulate"] = None
                    if poly_order is not None:
                        params["polyfilter"] = None
                        params["poly-order"] = poly_order
                    if poly_order_2d is not None:
                        params["polyfilter2D"] = None
                        params["poly-order2D"] = poly_order_2d
                    if ground_order is not None:
                        params["groundfilter"] = None
                        params["ground-order"] = ground_order

                    if flavor == "noise":
                        params["simulate-noise"] = None
                        params["hits"] = None
                        params["wcov"] = None
                        params["wcov-inv"] = None
                        params["MC-count"] = 8
                        walltime = "08:00:00"
                    elif flavor == "atmosphere":
                        parfiles += f" @atmosphere_{site}.par"
                        params["simulate-atmosphere"] = None
                        params["no-hits"] = None
                        params["no-wcov"] = None
                        params["no-wcov-inv"] = None
                        params["MC-count"] = 8
                        params["atm-cache"] = atm_cache
                        walltime = "08:00:00"
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
                            #"cmb-tensors" : "cmb_tensor_nest",
                            "cmb-tensors" : "cmb_tensor",
                            "foreground" : "combined_foregrounds",
                        }[flavor]
                        num = "0000"
                        if flavor == "cmb-tensors":
                            # cmb_tensor/4096/cmbs4_cmb_tensor_uKCMB_LAT-HFL1_nside4096_0000.fits
                            params["input-map"] = os.path.join(
                                "cmb_tensor",
                                str(nside),
                                "cmbs4_cmb_tensor_uKCMB_{}-{}_nside{}_{}.fits".format(
                                    telescope, band.replace("P", ""), nside, num
                                )
                            )
                        else:
                            params["input-map"] = os.path.join(
                                input_map_dir,
                                str(nside),
                                signal_name,
                                num,
                                "cmbs4_{}_uKCMB_{}-{}_nside{}_{}.fits".format(
                                    signal_name, telescope, band.replace("P", ""), nside, num
                                )
                            )
                        walltime = "02:00:00"
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
                            "\nlogfile=logs/{}.log\n".format(rootname),
                            "if [[ ! -e $logfile ]]; then",
                            '    echo "Writing $logfile at" `date`',
                            "    date > ${logfile}",    
                            "    srun -n $ntask -c $ncore --cpu_bind=cores \\",
                            "        toast_s4_sim.py {} \\".format(parfiles),
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
                            "    >> ${logfile} 2>&1",
                            "    date >> ${logfile}",
                            '    echo "Done with $logfile at" `date`',
                            "else",
                            '    echo "$logfile exists"',
                            "fi",
                        ]:
                            slurm.write(line + "\n")
