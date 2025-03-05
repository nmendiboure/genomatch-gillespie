""" Commands for the genomatchgp package """

import os
import shutil
from os.path import join
from docopt import docopt

from matplotlib.pyplot import plot
import yaml
import numpy as np
from mpi4py import MPI

from genomatchgp import methods
from genomatchgp.modelmaker import generate_gillespie_model
from genomatchgp.simulation import run


COMM = MPI.COMM_WORLD
RANK = COMM.Get_rank()
SIZE = COMM.Get_size()


class AbstractCommand:
    """Base class for the commands"""

    def __init__(self, command_args, global_args):
        """
        Initialize the commands.

        :param command_args: arguments of the command
        :param global_args: arguments of the program
        """
        self.args = docopt(self.__doc__, argv=command_args)
        self.global_args = global_args

    def execute(self):
        """Execute the commands"""
        raise NotImplementedError


class Run(AbstractCommand):

    """
    Run a single simulation.

    Usage:
        run <parameters> [--cells CELLS] [--output OUTPUT]

    Arguments:
        <parameters>  Path to the parameters .yaml file

    Options:
        -c CELLS , --cells CELLS            Number of cells (replicates)
        -o CELLS , --output OUTPUT          Output directory
    """

    def execute(self):

        yaml_path = self.args["<parameters>"]
        outdir = self.args["--output"]
        n_cells = int(self.args["--cells"])

        if outdir is None:
            outdir = os.path.join(os.path.dirname(__file__), "data", "output")

        with open(yaml_path, "r", encoding="utf-8") as file:
            data_yaml = yaml.safe_load(file)

        n_timepoints = data_yaml["timepoints"]
        n_points = n_timepoints // data_yaml["every"]
        
        if RANK == 0:
            my_model, index_to_species, model_uid = generate_gillespie_model(data_yaml)
            outdir = os.path.join(outdir, str(model_uid))
            os.makedirs(outdir, exist_ok=True)
            shutil.copy2(yaml_path, os.path.join(outdir, "params.yaml"))
            species_to_index = {v: k for k, v in index_to_species.items()}
        else:
            my_model = None
            model_uid = None
            species_to_index = None
            index_to_species = None
            outdir = None
        
        my_model = COMM.bcast(my_model, root=0)
        model_uid = COMM.bcast(model_uid, root=0)
        species_to_index = COMM.bcast(species_to_index, root=0)
        index_to_species = COMM.bcast(index_to_species, root=0)
        outdir = COMM.bcast(outdir, root=0)


        sims_per_rank = n_cells // SIZE
        start_idx = RANK * sims_per_rank
        end_idx = (RANK + 1) * sims_per_rank if RANK != SIZE - 1 else n_cells
        
        records_dir = join(outdir, "records")
        plots_dir = join(outdir, "plots")
        os.makedirs(records_dir, exist_ok=True)
        os.makedirs(plots_dir, exist_ok=True)
        for s in range(start_idx, end_idx):
            run(s, model_uid, species_to_index, my_model, data_yaml, records_dir)

        COMM.Barrier()

        if RANK == 0:
            all_results = []
            for s in range(n_cells):
                npz = np.load(join(records_dir, f"simulation_{s}.npz"))
                all_results.append({f : npz[f] for f in npz.files})
 

            agg_result_groups = methods.aggregate_groups(all_results)
            methods.plot_trajectories(agg_result_groups, os.path.join(plots_dir, "aggregated_trajectories.png"))
            
            convo_dlc = [0, 100, 200, 500, 1000]
            for c in convo_dlc:
                methods.plot_dlc(
                    agg_result_groups["DLC homologous"],
                    c,
                    os.path.join(plots_dir, f"aggregated_homologous_DLC_convo{c}.png"),
                )
