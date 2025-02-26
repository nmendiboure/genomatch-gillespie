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
    n_tp = 100000  # Final simulation time
    n_output_tp = n_tp // 10  # Number of output time points
    n_simu = 1000

    results_raw = []
    results_group = []
    results_dlc = []

    for i in range(n_simu):
        r.reset()

        # t_start = abs(int(random.normal(N // 16, N // 32)))
        # if t_start > N:
        #     t_start = N // 16

        t_start = 0
        t_end = n_tp

        s = r.simulate(t_start, n_tp, n_output_tp)
        results_raw.append(s)
        results_dlc.append([s[:, 20], s[:, 21]])
        if s[:, 22][-1] != 0:
            recombined_time = np.where(s[:, 22] > 0.)[0][0]
            results_dlc[i][0][recombined_time+1:] = 0
            results_dlc[i][1][recombined_time+1:] = 0

        results_group.append(methods.make_group(s))

        print(f"Simulation {i} completed.")

    agg_result_groups = methods.aggregate_groups(results_group)
    # methods.plot_trajectories(agg_result_groups)
    methods.plot_dlc(results_dlc)