"""
Main script for running the simulation
"""
import numpy as np
import tellurium as te
from numpy import random
from model import gillespie_model

import methods

if __name__ == "__main__":
    # Load the model from the imported module
    seed = 1999
    my_model = gillespie_model
    r = te.loada(my_model)
    r.integrator = "gillespie"
    r.integrator.seed = seed
    random.seed(seed)

    # Simulation parameters
    n_timepoints = 100000  #  N timepoints
    n_simu = 100 # N replciates
    every = 10 # capture every
    n_points = n_timepoints // every
    mu = n_timepoints // 16
    sigma = n_timepoints // 32


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
        r.S = 200
        start, end = timepoints_dyn[0], timepoints_dyn[-1]
        s_dyn = r.simulate(start, end, n_points_dyn)

        # concatenating the results
        s_total = np.zeros((23, n_points))
        s_total[0, :] = np.arange(0, n_timepoints, every)
        s_total[1, :] = np.concatenate(([200] * n_points_delay, s_dyn[:, 1]))

        for i in range(2, 23):
            s_total[i, :] = np.concatenate((s_delay[:, i], s_dyn[:, i]))

        results_raw.append(s_total)

        # DLC analysis
        results_dlc_homologous[s, :] = s_total[20, :]
        if s_total[22, :][-1] != 0:
            # recombination occured
            recomb_time = np.where(s_total[22, :] > 0.)[0][0]
            results_dlc_homologous[s][recomb_time+1:] = 0

        results_group.append(methods.make_group(s_total))

        print(f"Simulation {s} completed.")

    agg_result_groups = methods.aggregate_groups(results_group)
    # methods.plot_trajectories(agg_result_groups)
    methods.plot_dlc(results_dlc_homologous)