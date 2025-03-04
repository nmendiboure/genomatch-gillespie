"""Model generator for the Gillespie algorithm."""

import hashlib


def generate_gillespie_model(config):
    """
    Generates a Gillespie model based on a YAML configuration file.
    """
    
    data_string = "".join([str(value) for value in config.values()])
    uid = int(hashlib.sha256(data_string.encode()).hexdigest(), 16) % 2**32

    # Load parameters from the YAML file
    population_size = config["N"]
    homologous_fraction = config["f"]
    kon = config["kon"]
    koff1 = config["koff1"]
    kext = config["kext"]
    ktol = config["ktol"]
    kdloop = config["kdloop"]
    koff2 = config["koff2"]
    kre = config["kre"]

    model = f"""
    model rad51_recombination()
        // Parameters
        N = {population_size};
        f = {homologous_fraction};
        kon = {kon};
        koff1 = {koff1};
        kext = {kext};
        ktol = {ktol};
        kdloop = {kdloop};
        koff2 = {koff2};
        kre = {kre};
    
        // Species
        S = N;
    """

    species_index = {}  # Dictionary to map species names to indices
    species_index[0] = "Time"
    species_index[1] = "S"
    index = 2


    intermediates = config["intermediates"]

    # Declare homologous and heterologous complexes
    for label in ["HM", "HT"]:
        model += f"\n        // {label} Complexes\n"
        for i in intermediates:
            species_name = f"{label}{i}"
            model += f"        {species_name} = 0;\n"
            species_index[index] = species_name
            index += 1

    # D-loop species
    model += "\n        // D-loop species\n"
    for label in ["DHM", "DHT"]:
        model += f"        {label} = 0;\n"
        species_index[index] = label
        index += 1

    # Recombined state
    model += "\n        // Recombined state\n        R = 0;\n"
    species_index[index] = "R"
    index += 1

    reaction_id = 1

    # Association reactions
    model += "\n        // Association reactions\n"
    model += f"        R{reaction_id}: S -> HM8; kon * f * S;\n"
    reaction_id += 1
    model += f"        R{reaction_id}: S -> HT8; kon * (1 - f) * S;\n"
    reaction_id += 1

    # Dissociation reactions
    for label in ["HM", "HT"]:
        model += (
            "\n        // Dissociation reactions (homologous)\n"
            if label == "HM"
            else "\n        // Dissociation reactions (heterologous)\n"
        )
        for i in intermediates:
            c = min(4, sum(i > x for x in [8, 9, 12, 15]))
            model += f"        R{reaction_id}: {label}{i} -> S; koff1 / (1.4^{c}) * {label}{i};\n"
            reaction_id += 1

    # Extension reactions
    for label in ["HM", "HT"]:
        model += (
            "\n        // Extension reactions (homologous)\n"
            if label == "HM"
            else "\n        // Extension reactions (heterologous)\n"
        )
        for idx in range(len(intermediates) - 1):
            rate = "kext * ktol" if label == "HT" else "kext"
            model += f"        R{reaction_id}: {label}{intermediates[idx]} -> {label}{intermediates[idx+1]}; {rate} * {label}{intermediates[idx]};\n"
            reaction_id += 1

    # D-loop formation reactions
    dloop_powers = [
        (15, 24, 5),
        (24, 48, 4),
        (48, 96, 3),
        (96, 144, 2),
        (144, 192, 1),
    ]

    for label in ["HM", "HT"]:
        model += f"\n        // D-loop formation reactions ({'homologous' if label == 'HM' else 'heterologous'})\n"

        for i in reversed(intermediates):
            p = (
                0
                if i < 15
                else 1
                / (
                    10
                    ** next(
                        (x for min_i, max_i, x in dloop_powers if min_i <= i < max_i), 0
                    )
                )
            )

            model += f"        R{reaction_id}: {label}{i} -> D{label}; kdloop * {p} * {label}{i};\n"
            reaction_id += 1

    # Fate of D-loops
    model += "\n        // D-loop fate reactions\n"
    model += f"        R{reaction_id}: DHM -> S; koff2 * DHM;\n"
    reaction_id += 1
    model += f"        R{reaction_id}: DHT -> S; koff2 * DHT;\n"
    reaction_id += 1
    model += f"        R{reaction_id}: DHM -> R; kre * DHM;\n"
    reaction_id += 1
    model += f"        R{reaction_id}: DHT -> R; kre * DHT;\n    end\n"

    return model, species_index, uid