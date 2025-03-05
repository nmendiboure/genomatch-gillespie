import numpy as np
import tellurium as te
from mpi4py import MPI
from numpy import random

from genomatchgp import methods

COMM = MPI.COMM_WORLD
RANK = COMM.Get_rank()
SIZE = COMM.Get_size()


def run(
    uid: int,
    model_id: int,
    species_to_index: dict,
    my_model: str,
    params: dict,
    outdir: str,
):
    """
    Run a single simulation.
    """

    print(f"[Process {RANK}] :: SIMULATION {uid}", flush=True)

    # ensure different seeds
    seed0 = params["seed_zero"]
    seed = seed0 + uid + model_id

    n_timepoints = params["timepoints"]
    every = params["every"]
    k = params["gamma_k"]
    theta = params["gamma_theta"]

    random.seed(seed)
    r = te.loada(my_model)
    r.integrator = "gillespie"
    r.integrator.seed = seed

    n_points = n_timepoints // every
    n_species = len(species_to_index)

    delay = int(random.gamma(k, theta))
    delay = min(delay, n_timepoints - 1)

    n_points_delay = int(delay / every)
    n_points_dyn = int((n_timepoints - delay) / every)
    if n_points_delay + n_points_dyn != n_points:
        n_points_delay += 1

    timepoints_delay = np.linspace(0, delay, n_points_delay)
    timepoints_dyn = np.linspace(delay, n_timepoints, n_points_dyn)

    r.reset()
    # resection step
    r.S = 0
    start, end = timepoints_delay[0], timepoints_delay[-1]
    s_delay = r.simulate(start, end, n_points_delay)

    # synthesis step
    n = params["N"]
    r.S = n
    start, end = timepoints_dyn[0], timepoints_dyn[-1]
    s_dyn = r.simulate(start, end, n_points_dyn)

    s_total = np.zeros((n_species, n_points))
    s_total[0, :] = np.arange(0, n_timepoints, every)
    s_total[1, :] = np.concatenate(([n] * n_points_delay, s_dyn[:, 1]))

    for i in range(2, n_species):
        s_total[i, :] = np.concatenate((s_delay[:, i], s_dyn[:, i]))

    result_group = methods.make_group(s_total, species_to_index)

    np.savez_compressed(
        f"{outdir}/simulation_{uid}.npz",
        **result_group,
    )
