"""
This module contains the core functions for the simulation of the DNA binding model.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

def make_group(result: np.ndarray, species_to_index) -> dict:
    """
    Group the results of the simulation into different categories.
    """
    def sum_species(species_list):
        res = np.zeros(result.shape[1])
        for species in species_list:
            res += result[species_to_index[species], :]
        return res
    
    all_homo = [k for k in species_to_index if k.startswith("HM")]
    all_hetero = [k for k in species_to_index if k.startswith("HT")]

    homo_8_9 = sum_species(["HM8", "HM9"])
    homo_12_21 = sum_species([h for h in all_homo if 12 <= int(h[2:]) <= 21])
    homo_24_48 = sum_species([h for h in all_homo if 24 <= int(h[2:]) <= 48])
    homo_51_144 = sum_species([h for h in all_homo if 51 <= int(h[2:]) <= 144])
    homo_144_plus = sum_species([h for h in all_homo if  144 < int(h[2:])])

    hetero_8_9 = sum_species(["HT8", "HT9"])
    hetero_12_21 = sum_species([h for h in all_hetero if 12 <= int(h[2:]) <= 21])
    hetero_24_48 = sum_species([h for h in all_hetero if 24 <= int(h[2:]) <= 48])
    hetero_51_144 = sum_species([h for h in all_hetero if 51 <= int(h[2:]) <= 144])
    hetero_144_plus = sum_species([h for h in all_hetero if  144 < int(h[2:])])

    all_8_9 = homo_8_9 + hetero_8_9
    all_12_21 = homo_12_21 + hetero_12_21
    all_24_48 = homo_24_48 + hetero_24_48
    all_51_144 = homo_51_144 + hetero_51_144
    all_144_plus = homo_144_plus + hetero_144_plus

    homo_dloop = result[species_to_index["DHM"], :]
    hetero_dloop = result[species_to_index["DHT"], :]
    recombined = result[species_to_index["R"], :]

    occupied = all_8_9 + all_12_21 + all_24_48 + all_51_144 + all_144_plus + homo_dloop + hetero_dloop + recombined

    res: dict = {
        "time": result[0, :],
        "free sites": result[1, :],
        "occupied sites": np.sum(result[2:23, :], axis=1),
        "homo 8-9 nts": homo_8_9,
        "homo 12-21 nts": homo_12_21,
        "homo 24-48 nts": homo_24_48,
        "homo 51-144 nts": homo_51_144,
        "homo 144+ nts": homo_144_plus, 
        "hetero 8-9 nts": hetero_8_9,
        "hetero 12-21 nts": hetero_12_21,
        "hetero 24-48 nts": hetero_24_48,
        "hetero 51-144 nts": hetero_51_144,
        "hetero 144+ nts": hetero_144_plus,
        "all": occupied,
        "all associations 8-9 nts": all_8_9,
        "all associations 12-21 nts": all_12_21,
        "all associations 24-48 nts": all_24_48,
        "all associations 51-144 nts": all_51_144,
        "all associations 144+ nts": all_144_plus,
        "D-loop homologies": homo_dloop,
        "D-loop heterologies": hetero_dloop,
        "Recombined": recombined
    }

    return res


def aggregate_groups(groups: list[dict]):
    n = len(groups)
    sum_groups: dict = {}

    for i in range(n):
        group = groups[i]
        for k, v in group.items():
            if k not in sum_groups:
                sum_groups[k] = v
            else:
                sum_groups[k] += v
    
    agg_groups = {k: v / n for k, v in sum_groups.items()}
    return agg_groups



def plot_trajectories(one_res, outpath="", show=False):
    """
    Plot the trajectories of the different molecules in the model.
    """

    t = one_res["time"]

    #Create a 3x2 grid of subplots with more space
    fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(14, 12), constrained_layout=True)

    # Define plot settings
    title_font = {"fontsize": 12}
    label_font = {"fontsize": 12}
    legend_font = fm.FontProperties(size=8)

    # Plot 0: Free Binding Sites
    ax0 = axs[0, 0]
    ax0.plot(t, one_res["free sites"], color="blue", label="S (Free Sites)")
    ax0.plot(t, one_res["all"], color="red", label="Total Occupancy")
    ax0.set_title("Filament Binding Dynamics", **title_font)
    ax0.set_xlabel("T", **label_font)
    ax0.set_ylabel("N", **label_font)
    ax0.legend(prop=legend_font)
    ax0.grid(True, linestyle="--", alpha=0.6)

    # Plot 1: Homologous Complexes
    ax1 = axs[0, 1]
    ax1.set_title("Homologous pre-synaptic complexes", **title_font)
    ax1.set_xlabel("T", **label_font)
    ax1.set_ylabel("N", **label_font)
    ax1.grid(True, linestyle="--", alpha=0.6)

    # Plot 2: Heterologous Complexes
    ax2 = axs[1, 0]
    ax2.set_title("Heterologous pre-synaptic complexes", **title_font)
    ax2.set_xlabel("T", **label_font)
    ax2.set_ylabel("N", **label_font)
    ax2.grid(True, linestyle="--", alpha=0.6)

    # Plot 3: Total Associations by Binding Length
    ax3 = axs[1, 1]
    ax3.set_title("All pre-synaptic complexes", **title_font)
    ax3.set_xlabel("T", **label_font)
    ax3.set_ylabel("N", **label_font)
    ax3.grid(True, linestyle="--", alpha=0.6)

    # Plot 4: D-loop Species
    ax4 = axs[2, 0]
    ax4.set_title("D-loop", **title_font)
    ax4.set_xlabel("Time", **label_font)
    ax4.set_ylabel("N", **label_font)
    ax4.grid(True, linestyle="--", alpha=0.6)

    # Plot 5: Recombined State
    ax5 = axs[2, 1]
    ax5.set_title("Recombined State", **title_font)
    ax5.set_xlabel("Time", **label_font)
    ax5.set_ylabel("N", **label_font)
    ax5.grid(True, linestyle="--", alpha=0.6)

    # Loop through the data and add to appropriate subplot
    for k, v in one_res.items():
        if k.startswith("homo"):
            ax1.plot(t, v, label=k)
        elif k.startswith("hetero"):
            ax2.plot(t, v, label=k)
        elif k.startswith("all associations"):
            ax3.plot(t, v, label=k)
        elif k.startswith("D-loop"):
            ax4.plot(t, v, label=k)
        elif k.startswith("Recombined"):
            ax5.plot(t, v, label=k)

    # Add legends
    ax1.legend(prop=legend_font)
    ax2.legend(prop=legend_font)
    ax3.legend(prop=legend_font)
    ax4.legend(prop=legend_font)
    ax5.legend(prop=legend_font)

    fig.tight_layout()

    if show:
        plt.show()
    
    if outpath != "":
        plt.savefig(outpath, format="png")

    plt.close(fig)


def plot_dlc(dlc_results: np.ndarray, convo: int = 100, outpath="", show=False):
    n, t = dlc_results.shape

    aggr_dlc = np.zeros(t)
    for i in range(n):
        aggr_dlc += dlc_results[i]
 
    aggr_dlc /= n


    if convo > 0:
        aggr_dlc = np.convolve(aggr_dlc, np.ones(convo) / convo, mode='same')

    xt = np.arange(0, t)
    fig, ax = plt.subplots()
    xt = np.arange(t)
    ax.plot(xt, aggr_dlc)

    if show:
        plt.show()

    if outpath:
        fig.savefig(outpath, format="png")

    plt.close(fig)