import tellurium as te
import numpy as np
import matplotlib.pyplot as plt
from model import my_model # This module contains the 'model' string with your Antimony model


def resection(self, timepoints: int, tot_time_h: int, exo1_speed_bph: int) -> int:
    """
    Simulate the resection process of the filament using a normal distribution.
    """
    # 1. Calculate 1 hour in timepoints
    # 2. Calculate the speed of the Exo1 enzyme in timepoints
    # 3. Calculate the timepoints needed for the resection process of the filament
    # 4. Simulate the end of the resection process using a normal distribution, with
    #   sigma being 1/4 h for the DSB induction.
    one_hour_timepoints = timepoints / tot_time_h
    exo1_speed_timepoints = exo1_speed_bph / one_hour_timepoints
    resection_timepoints = int(self.size / exo1_speed_timepoints)

    mu = resection_timepoints
    sigma = 1 / 4 * one_hour_timepoints
    resection_end = abs(int(np.random.normal(mu, sigma)))
    if resection_end > timepoints:
        resection_end = timepoints
    return resection_end


def make_group(result_array: np.ndarray) -> dict:
    """
    Group the results of the simulation into different categories.
    """
    # Extract time
    time = result_array[:, 0]

    # Free binding sites
    S = result_array[:, 1]

    # Homologous complexes
    HM8 = result_array[:, 2]
    HM9 = result_array[:, 3]
    HM12 = result_array[:, 4]
    HM15 = result_array[:, 5]
    HM24 = result_array[:, 6]
    HM48 = result_array[:, 7]
    HM96 = result_array[:, 8]
    HM192 = result_array[:, 9]
    HM384 = result_array[:, 10]

    # Heterologous complexes
    HT8 = result_array[:, 11]
    HT9 = result_array[:, 12]
    HT12 = result_array[:, 13]
    HT15 = result_array[:, 14]
    HT24 = result_array[:, 15]
    HT48 = result_array[:, 16]
    HT96 = result_array[:, 17]
    HT192 = result_array[:, 18]
    HT384 = result_array[:, 19]

    # D-loops and recombined state
    DHM = result_array[:, 20]
    DHT = result_array[:, 21]
    R = result_array[:, 22]

    # Grouped complexes
    HM_8_9 = HM8 + HM9
    HM_12_24 = HM12 + HM15 + HM24
    HM_48_384 = HM48 + HM96 + HM192 + HM384

    HT_8_9 = HT8 + HT9
    HT_12_24 = HT12 + HT15 + HT24
    HT_48_384 = HT48 + HT96 + HT192 + HT384

    # Total associations by binding length (F = HM + HT)
    F_8_9 = HM_8_9 + HT_8_9
    F_12_24 = HM_12_24 + HT_12_24
    F_48_384 = HM_48_384 + HT_48_384

    res = {
        "time": time,
        "S": S,
        "HM_8_9": HM_8_9,
        "HM_12_24": HM_12_24,
        "HM_48_384": HM_48_384,
        "HT_8_9": HT_8_9,
        "HT_12_24": HT_12_24,
        "HT_48_384": HT_48_384,
        "F_8_9": F_8_9,
        "F_12_24": F_12_24,
        "F_48_384": F_48_384,
        "DHM": DHM,
        "DHT": DHT,
        "R": R
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
    MODEL = my_model
    r = te.loada(MODEL)
    r.integrator = "gillespie"
    r.integrator.seed = 1999

    # Simulation parameters
    N = 100000  # Final simulation time
    N_outputs = N // 10  # Number of output time points

    # Run the Gillespie stochastic simulation
    s = r.simulate(0, N, N_outputs)

    group = make_group(s)
    plot_trajecotries(group)

