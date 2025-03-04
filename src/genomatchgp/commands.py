""" Commands for the genomatchgp package """

import os
import shutil
from os.path import join
from docopt import docopt

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
        
        local_results_group = []
        local_results_dlc_homologous = np.zeros((end_idx - start_idx, n_points))
        
        for s in range(start_idx, end_idx):
            result_group, dlc_homologous = run(s, model_uid, species_to_index, my_model, data_yaml)
            local_results_group.append(result_group)
            local_results_dlc_homologous[s - start_idx, :] = dlc_homologous
        
        all_results_group = COMM.gather(local_results_group, root=0)
        all_results_dlc_homologous = COMM.gather(local_results_dlc_homologous, root=0)
        
        if RANK == 0:
            results_group = [item for sublist in all_results_group for item in sublist]
            results_dlc_homologous = np.vstack(all_results_dlc_homologous)
            
            agg_result_groups = methods.aggregate_groups(results_group)
            methods.plot_trajectories(agg_result_groups, os.path.join(outdir, "Aggregated_trajectories.png"))
            
            convo_dlc = [0, 100, 200, 500, 1000]
            for c in convo_dlc:
                methods.plot_dlc(
                    results_dlc_homologous,
                    c,
                    os.path.join(outdir, f"Aggregated_homologous_DLC_convo{c}.png"),
                )
