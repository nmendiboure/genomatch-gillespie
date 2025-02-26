"""
This module contains the core functions for the simulation of the DNA binding model.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

def make_group(result: np.ndarray) -> dict:
    """
    Group the results of the simulation into different categories.
    """
    
    homo_8_9 = result[:, 2] + result[:, 3]
    homo_12_15 = result[:, 4] + result[:, 5]
    homo_24_48 = result[:, 6] + result[:, 7]
    homo_96_384 = result[:, 8] + result[:, 9] + result[:, 10]
    hetero_8_9 = result[:, 11] + result[:, 12]
    hetero_12_15 = result[:, 13] + result[:, 14]
    hetero_24_48 = result[:, 15] + result[:, 16]
    hetero_96_384 = result[:, 17] + result[:, 18] + result[:, 19]
    all_9_8 = hetero_8_9 + homo_8_9
    all_12_15 = hetero_12_15 + homo_12_15
    all_24_48 = hetero_24_48 + homo_24_48
    all_96_384 = hetero_96_384 + homo_96_384
    homo_dloop = result[:, 20]
    hetero_dloop = result[:, 21]
    recombined = result[:, 22]
    occupied = all_9_8 + all_12_15 + all_24_48 + all_96_384 + hetero_dloop + homo_dloop + recombined

    res: dict = {
        "time": result[:, 0],
        "free sites": result[:, 1],
        "occupied sites": np.sum(result[:, 2:23], axis=1),
        "homo 8-9 nts": homo_8_9,
        "homo 12-15 nts": homo_24_48,
        "homo 24-48 nts": homo_24_48,
        "homo 96-384 nts": homo_96_384, 
        "hetero 8-9 nts": hetero_8_9,
        "hetero 12-15 nts": hetero_12_15,
        "hetero 24-48 nts": hetero_24_48,
        "hetero 96-384 nts": hetero_96_384,
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



def plot_trajectories(one_res):
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
    plt.show()


def plot_dlc(dlc_results: list):

    n = len(dlc_results)
    t = len(dlc_results[0][0])
    xt = np.arange(0, t, 1)


    #Create a 3x2 grid of subplots with more space
    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(14, 12))

    # Define plot settings
    title_font = {"fontsize": 12}
    label_font = {"fontsize": 12}

    aggr_dlc = [np.zeros(t), np.zeros(t)]
    for i in range(n):
        aggr_dlc[0] += dlc_results[i][0]
        aggr_dlc[1] += dlc_results[i][1]

    aggr_dlc[0] /= n
    aggr_dlc[1] /= n

    if aggr_dlc[0].max() > 0:
        #  make a convolution of the data
        aggr_dlc[0] = np.convolve(aggr_dlc[0], np.ones(200) / 200, mode='same')


    # Plot 0: Free Binding Sites
    ax0 = axs[0]
    ax0.plot(xt, aggr_dlc[0], color="blue")
    ax0.set_title("Homologous D-loops %", **title_font)
    ax0.set_xlabel("T", **label_font)
    ax0.set_ylabel("N", **label_font)
    ax0.grid(True, linestyle="--", alpha=0.6)

    # Plot 1: Homologous Complexes
    ax1 = axs[1]
    ax1.plot(xt, aggr_dlc[1], color="red")
    ax1.set_title("Heterologous D-loops %", **title_font)
    ax1.set_xlabel("T", **label_font)
    ax1.set_ylabel("N", **label_font)
    ax1.grid(True, linestyle="--", alpha=0.6)

    fig.tight_layout()
    plt.show()



