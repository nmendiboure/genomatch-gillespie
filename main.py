import tellurium as te
import numpy as np
from numpy import random
import matplotlib.pyplot as plt
from model import my_model # This module contains the 'model' string with your Antimony model


def make_group(result: np.ndarray) -> dict:
    """
    Group the results of the simulation into different categories.
    """
    
    homo_8_9 = result[:, 2] + result[:, 3]
    homo_12_15 = result[:, 4] + result[:, 5]
    homo_24_48 = result[:, 6] + result[:, 7]
    homo_96_384 = result[:, 8] + result[:, 9] + result[:, 10],
    hetero_8_9 = result[:, 11] + result[:, 12]
    hetero_12_15 = result[:, 13] + result[:, 14]
    hetero_24_48 = result[:, 15] + result[:, 16]
    hetero_96_384 = result[:, 17] + result[:, 18] + result[:, 19]
    all_9_8 = hetero_8_9 + homo_8_9
    all_12_15 = hetero_12_15 + homo_12_15
    all_24_48 = hetero_24_48 + homo_24_48
    all_96_384 = hetero_96_384 + homo_96_384
    hetero_dloop = result[:, 20]
    homo_dloop = result[:, 21]
    recombined = result[:, 22]
    occupied = all_9_8 + all_12_15 + all_24_48 + all_96_384 + hetero_dloop + homo_dloop + recombined


    res = {
        "time": result[:, 0],
        "free sites": result[:, 1],
        "occupied sites": np.sum(result[:, 2:23], axis=1),
        "homologies 8-9 nts": homo_8_9,
        "homologies 12-15 nts": homo_24_48,
        "homologies 24-48 nts": homo_24_48,
        "homologies 96-384 nts": homo_96_384, 
        "heterologies 8-9 nts": hetero_8_9,
        "heterologies 12-15 nts": hetero_12_15,
        "heterologies 24-48 nts": hetero_24_48,
        "heterologies 96-384 nts": hetero_96_384,
        "all": occupied,
        "all associations 8-9 nts": all_9_8,
        "all associations 12-15 nts": all_12_15,
        "all associations 24-48 nts": all_24_48,
        "all associations 96-384 nts": all_96_384,
        "D-loop homologies": homo_dloop,
        "D-loop heterologies": hetero_dloop,
        "Recombined": recombined
    }

    return res


def plot_trajecotries(res):
    """
    Plot the trajectories of the different molecules in the model.
    """

    t = res["time"]
    S = res["S"]
    HM_8_9 = res["HM_8_9"]
    HM_12_24 = res["HM_12_24"]
    HM_48_384 = res["HM_48_384"]
    HT_8_9 = res["HT_8_9"]
    HT_12_24 = res["HT_12_24"]
    HT_48_384 = res["HT_48_384"]
    F_8_9 = res["F_8_9"]
    F_12_24 = res["F_12_24"]
    F_48_384 = res["F_48_384"]
    DHM = res["DHM"]
    DHT = res["DHT"]
    R = res["R"]

    # Create a 3x2 grid of subplots
    _, axs = plt.subplots(nrows=3, ncols=2, figsize=(24, 22))
    plt.subplots_adjust(hspace=0.35, wspace=0.3)

    # Plot 1: Free Binding Sites
    ax = axs[0, 0]
    ax.plot(t, S, color="blue", label="S (Free Sites)")
    ax.set_title("Free Binding Sites")
    ax.set_xlabel("Time")
    ax.set_ylabel("Number of Sites")
    ax.legend()

    # Plot 2: Homologous Complexes
    ax = axs[0, 1]
    ax.plot(t, HM_8_9, color="red", label="HM 8-9 nt")
    ax.plot(t, HM_12_24, color="orange", label="HM 12-24 nt")
    ax.plot(t, HM_48_384, color="magenta", label="HM 48-384 nt")
    ax.set_title("Homologous Complexes")
    ax.set_xlabel("Time")
    ax.set_ylabel("Number of Complexes")
    ax.legend()

    # Plot 3: Heterologous Complexes
    ax = axs[1, 0]
    ax.plot(t, HT_8_9, color="red", linestyle="--", label="HT 8-9 nt")
    ax.plot(t, HT_12_24, color="orange", linestyle="--", label="HT 12-24 nt")
    ax.plot(t, HT_48_384, color="magenta", linestyle="--", label="HT 48-384 nt")
    ax.set_title("Heterologous Complexes")
    ax.set_xlabel("Time")
    ax.set_ylabel("Number of Complexes")
    ax.legend()

    # Plot 4: Total Associations by Binding Length
    ax = axs[1, 1]
    ax.plot(t, F_8_9, color="green", label="F 8-9 nt total")
    ax.plot(t, F_12_24, color="blue", label="F 12-24 nt total")
    ax.plot(t, F_48_384, color="black", label="F 48-384 nt total")
    ax.set_title("Total Associations by Binding Length")
    ax.set_xlabel("Time")
    ax.set_ylabel("Number of Associations")
    ax.legend()

    # Plot 5: D-loop Species
    ax = axs[2, 0]
    ax.plot(t, DHM, color="cyan", label="DHM (Homologous D-loop)")
    ax.plot(t, DHT, color="grey", label="DHT (Heterologous D-loop)")
    ax.set_title("D-loop Species")
    ax.set_xlabel("Time")
    ax.set_ylabel("Number of D-loops")
    ax.legend()

    # Plot 6: Recombined State
    ax = axs[2, 1]
    ax.plot(t, R, color="black", label="R (Recombined)")
    ax.set_title("Recombined State")
    ax.set_xlabel("Time")
    ax.set_ylabel("Number of Recombinations")
    ax.legend()

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    # Load the model from the imported module
    SEED = 1999
    MODEL = my_model
    r = te.loada(MODEL)
    r.integrator = "gillespie"
    r.integrator.seed = SEED
    random.seed(SEED)

    # Simulation parameters
    N = 100000  # Final simulation time
    N_OUTPUTS = N // 10  # Number of output time points
    N_SIMULATIONS = 100

    
    results = []

    for i in range(N_SIMULATIONS):
        r.reset()

        t_start = abs(int(random.normal(N // 16, N // 32)))
        if t_start > N:
            t_start = N // 16

        s = r.simulate(t_start, N, N_OUTPUTS)
        results.append(s)
        print(f"Simulation {i} completed.")




    # group = make_group(s)
    # plot_trajecotries(group)

