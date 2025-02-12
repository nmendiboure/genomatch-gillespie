import tellurium as te
import matplotlib.pyplot as plt
import model  # This module contains the 'model' string with your Antimony model

if __name__ == "__main__":
    # Load the model from the imported module
    model_str = model.model
    r = te.loada(model_str)

    # Simulation parameters
    final_time = 100000  # Final simulation time
    num_points = 10000  # Number of output time points

    # Run the Gillespie stochastic simulation
    result = r.gillespie(0, final_time, num_points)

    # Extract time
    time = result[:, 0]

    # Free binding sites
    S = result[:, 1]

    # Homologous complexes
    HM8 = result[:, 2]
    HM9 = result[:, 3]
    HM12 = result[:, 4]
    HM15 = result[:, 5]
    HM24 = result[:, 6]
    HM48 = result[:, 7]
    HM96 = result[:, 8]
    HM192 = result[:, 9]
    HM384 = result[:, 10]

    # Heterologous complexes
    HT8 = result[:, 11]
    HT9 = result[:, 12]
    HT12 = result[:, 13]
    HT15 = result[:, 14]
    HT24 = result[:, 15]
    HT48 = result[:, 16]
    HT96 = result[:, 17]
    HT192 = result[:, 18]
    HT384 = result[:, 19]

    # D-loops and recombined state
    DHM = result[:, 20]
    DHT = result[:, 21]
    R = result[:, 22]

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

    # Create a 3x2 grid of subplots
    fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(24, 22))
    plt.subplots_adjust(hspace=0.35, wspace=0.3)

    # Plot 1: Free Binding Sites
    ax = axs[0, 0]
    ax.plot(time, S, color='blue', label='S (Free Sites)')
    ax.set_title('Free Binding Sites')
    ax.set_xlabel('Time')
    ax.set_ylabel('Number of Sites')
    ax.legend()

    # Plot 2: Homologous Complexes
    ax = axs[0, 1]
    ax.plot(time, HM_8_9, color='red', label='HM 8-9 nt')
    ax.plot(time, HM_12_24, color='orange', label='HM 12-24 nt')
    ax.plot(time, HM_48_384, color='magenta', label='HM 48-384 nt')
    ax.set_title('Homologous Complexes')
    ax.set_xlabel('Time')
    ax.set_ylabel('Number of Complexes')
    ax.legend()

    # Plot 3: Heterologous Complexes
    ax = axs[1, 0]
    ax.plot(time, HT_8_9, color='red', linestyle='--', label='HT 8-9 nt')
    ax.plot(time, HT_12_24, color='orange', linestyle='--', label='HT 12-24 nt')
    ax.plot(time, HT_48_384, color='magenta', linestyle='--', label='HT 48-384 nt')
    ax.set_title('Heterologous Complexes')
    ax.set_xlabel('Time')
    ax.set_ylabel('Number of Complexes')
    ax.legend()

    # Plot 4: Total Associations by Binding Length
    ax = axs[1, 1]
    ax.plot(time, F_8_9, color='green', label='F 8-9 nt total')
    ax.plot(time, F_12_24, color='blue', label='F 12-24 nt total')
    ax.plot(time, F_48_384, color='black', label='F 48-384 nt total')
    ax.set_title('Total Associations by Binding Length')
    ax.set_xlabel('Time')
    ax.set_ylabel('Number of Associations')
    ax.legend()

    # Plot 5: D-loop Species
    ax = axs[2, 0]
    ax.plot(time, DHM, color='cyan', label='DHM (Homologous D-loop)')
    ax.plot(time, DHT, color='grey', label='DHT (Heterologous D-loop)')
    ax.set_title('D-loop Species')
    ax.set_xlabel('Time')
    ax.set_ylabel('Number of D-loops')
    ax.legend()

    # Plot 6: Recombined State
    ax = axs[2, 1]
    ax.plot(time, R, color='black', label='R (Recombined)')
    ax.set_title('Recombined State')
    ax.set_xlabel('Time')
    ax.set_ylabel('Number of Recombinations')
    ax.legend()

    plt.tight_layout()
    plt.show()
