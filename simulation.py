"""
This script runs the simulation of the model.
"""

import os
import shutil
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import tellurium as te
from numpy import random

import methods
from modelmaker import generate_gillespie_model


def run_simulation(
    s, species_to_index, my_model, S, n_timepoints, every, mu, sigma, outdir
):
    """
    Run a single simulation.
    """
    seed = 1999 + s  # Ensure different seeds
    random.seed(seed)
    r = te.loada(my_model)
    r.integrator = "gillespie"
    r.integrator.seed = seed

    n_points = n_timepoints // every
    n_species = len(species_to_index)

    delay = abs(int(random.normal(mu, sigma)))
    n_points_delay = int(delay / every)
    n_points_dyn = int((n_timepoints - delay) / every)
    if n_points_delay + n_points_dyn != n_points:
        n_points_delay += 1

    timepoints_delay = np.linspace(0, delay, n_points_delay)
    timepoints_dyn = np.linspace(delay, n_timepoints, n_points_dyn)

    r.reset()
    r.S = 0
    start, end = timepoints_delay[0], timepoints_delay[-1]
    s_delay = r.simulate(start, end, n_points_delay)

    r.S = S
    start, end = timepoints_dyn[0], timepoints_dyn[-1]
    s_dyn = r.simulate(start, end, n_points_dyn)

    s_total = np.zeros((n_species, n_points))
    s_total[0, :] = np.arange(0, n_timepoints, every)
    s_total[1, :] = np.concatenate(([200] * n_points_delay, s_dyn[:, 1]))

    for i in range(2, n_species):
        s_total[i, :] = np.concatenate((s_delay[:, i], s_dyn[:, i]))

    dmh_idx = species_to_index["DHM"]
    recomb_idx = species_to_index["R"]
    results_dlc_homologous = s_total[dmh_idx, :]

    if s_total[recomb_idx, :][-1] != 0:
        recomb_time = np.where(s_total[recomb_idx, :] > 0.0)[0][0]
        results_dlc_homologous[recomb_time + 1 :] = 0

    result_group = methods.make_group(s_total, species_to_index)

    print(f"Simulation {s} done.")

    return s, result_group, results_dlc_homologous


def main():
    """
    Main function to run the simulation.
    """
    yaml_path = "params.yaml"
    output = "model.txt"
    intermediates = [8] + [i for i in range(9, 96, 3)] + [96, 144, 192, 240, 384]
    my_model, index_to_species, uid = generate_gillespie_model(
        yaml_path, intermediates, output
    )

    outdir = os.path.join(os.path.dirname(__file__), "output", str(uid))
    os.makedirs(outdir, exist_ok=True)
    shutil.copy2(yaml_path, os.path.join(outdir, "params.yaml"))
    shutil.copy2(output, os.path.join(outdir, "model.txt"))

    species_to_index = {v: k for k, v in index_to_species.items()}
    S = te.loada(my_model).S

    n_timepoints = 500000
    n_simu = 1000
    n_threads = 24
    every = 10
    mu = n_timepoints // 16
    sigma = n_timepoints // 32

    results_group = []
    results_dlc_homologous = np.zeros((n_simu, n_timepoints // every))

    with ProcessPoolExecutor(max_workers=n_threads) as executor:
        futures = [
            executor.submit(
                run_simulation,
                s,
                species_to_index,
                my_model,
                S,
                n_timepoints,
                every,
                mu,
                sigma,
                outdir,
            )
            for s in range(n_simu)
        ]

        for future in futures:
            s, result_group, dlc_homologous = future.result()
            results_group.append(result_group)
            results_dlc_homologous[s, :] = dlc_homologous

    agg_result_groups = methods.aggregate_groups(results_group)
    methods.plot_trajectories(
        agg_result_groups, os.path.join(outdir, "Aggregated_trajectories.png")
    )

    convo_dlc = [0, 100, 200, 500, 1000]
    for c in convo_dlc:
        methods.plot_dlc(
            results_dlc_homologous,
            c,
            os.path.join(outdir, f"Aggregated_homologous_DLC_convo{c}.png"),
        )


if __name__ == "__main__":
    main()
