"""
Main script for running the simulation
"""
import os
import shutil
import numpy as np
import tellurium as te
from numpy import random
from modelmaker import generate_gillespie_model
import methods

_DATA_ = os.path.__file__

if __name__ == "__main__":

    yaml_path = "params.yaml"
    output = "model.txt"
    intermediates = [8] + [i for i in range(9, 384, 3)]
    my_model, index_to_species, uid = generate_gillespie_model(yaml_path, intermediates, output)

    outdir = os.path.join(os.path.dirname(__file__), "output", str(uid))
    os.makedirs(outdir, exist_ok=True)
    shutil.copy2(yaml_path, os.path.join(outdir, "params.yaml"))
    shutil.copy2(output, os.path.join(outdir, "model.txt"))
    

    species_to_index = {v: k for k, v in index_to_species.items()}

    # Load the model from the imported module
    seed = 1999
    r = te.loada(my_model)
    r.integrator = "gillespie"
    r.integrator.seed = seed
    random.seed(seed)
    S = r.S

    # Simulation parameters
    n_timepoints = 100000  #  N timepoints
    n_simu = 500 # N replciates
    every = 10 # capture every
    n_points = n_timepoints // every
    mu = n_timepoints // 16
    sigma = n_timepoints // 32

    n_species = len(species_to_index )

    results_raw = []
    results_group = []
    results_dlc_homologous = np.zeros((n_simu, n_points))

    for s in range(n_simu):
        delay = abs(int(random.normal(mu, sigma))) # resection delay (gaussian)
        n_points_delay = int(delay / every)
        n_points_dyn = int((n_timepoints - delay) / every)
        if n_points_delay + n_points_dyn != n_points:
            n_points_delay += 1

        timepoints_delay = np.linspace(0, delay, n_points_delay)
        timepoints_dyn = np.linspace(delay, n_timepoints, n_points_dyn)


        r.reset()

        # modeling the resection step, delaying the homology search
        r.S = 0
        start, end = timepoints_delay[0], timepoints_delay[-1]
        s_delay = r.simulate(start, end, n_points_delay)

        # simualtion dynamics
        r.S = S
        start, end = timepoints_dyn[0], timepoints_dyn[-1]
        s_dyn = r.simulate(start, end, n_points_dyn)

        # concatenating the results
        s_total = np.zeros((n_species, n_points))
        s_total[0, :] = np.arange(0, n_timepoints, every)
        s_total[1, :] = np.concatenate(([200] * n_points_delay, s_dyn[:, 1]))

        for i in range(2, n_species):
            s_total[i, :] = np.concatenate((s_delay[:, i], s_dyn[:, i]))

        results_raw.append(s_total)

        # DLC analysis
        dmh_idx = species_to_index["DHM"]
        results_dlc_homologous[s, :] = s_total[dmh_idx, :]


        recomb_idx = species_to_index["R"]
        if s_total[recomb_idx, :][-1] != 0:
            # recombination occured
            recomb_time = np.where(s_total[recomb_idx, :] > 0.)[0][0]
            results_dlc_homologous[s][recomb_time+1:] = 0

        results_group.append(methods.make_group(s_total, species_to_index))
        # methods.plot_trajectories(results_group[-1], os.path.join(outdir, f"S{s}_trajectories.png"))
        print(f"Simulation {s} completed.")

    agg_result_groups = methods.aggregate_groups(results_group)
    methods.plot_trajectories(agg_result_groups, os.path.join(outdir, "Aggregated_trajectories.png"))

    convo_dlc = [0, 100, 200, 500, 1000]
    for c in convo_dlc:
        methods.plot_dlc(results_dlc_homologous, c, os.path.join(outdir, f"Aggregated_homologous_DLC_convo{c}.png"))