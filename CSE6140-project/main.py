#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
File: main.py
Description:
    This script runs algorithms and sets parameters given 4 command line arguments:
        i) the filename of a dataset,
        ii) the method to use (BF: brute force, Approx: approximation algorithm, LS: local search),
        iii) the cut-off time (in seconds),
        iv) a random seed (only used for local search algorithm).

Author: Kaichen Shen
Email: kshen77@gatech.edu
"""


import argparse
import os
import multiprocessing
import time
from exact_algorithm import exact_algorithm
from approximate_algorithm import approximate_algorithm
from simulated_annealing import simulated_annealing
from evolutionary_algorithm import evolutionary_algorithm


def main():
    # set command line arguments
    parser = argparse.ArgumentParser(description="TSP Solver Executable")
    parser.add_argument("-inst", required=True, help="Input TSP file name")
    parser.add_argument("-alg", required=True, choices=["BF", "Approx", "LS", "EA"], help="Algorithm to use: BF (brute force), Approx (approximation), LS (local search), EA (evolutionary algorithm)")
    parser.add_argument("-time", required=True, type=int, help="Cut-off time in seconds")
    parser.add_argument("-seed", type=int, help="Random seed (only for local search anf evolutionary algorithm)")
    args = parser.parse_args()

    # extract command line arguments
    filename = args.inst
    algorithm = args.alg
    cutoff = args.time
    seed = args.seed

    # check the existence of TSP file
    if not os.path.exists(filename):
        print(f"Error: File {filename} not found.")
        return

    # run corresponding algorithm
    if algorithm == "BF":
        start_time = time.perf_counter()
        tour, best_quality = exact_algorithm(filename, cutoff)
        end_time = time.perf_counter()
    elif algorithm == "Approx":
        start_time = time.perf_counter()
        tour, best_quality = approximate_algorithm(filename) # for approximation algorithm, "cutoff" can be omitted
        end_time = time.perf_counter()
    elif algorithm == "LS":
        if seed is None:
            print("Error: Seed is required for local search algorithm.")
            return
        start_time = time.perf_counter()
        tour, best_quality = simulated_annealing(filename, cutoff, seed)
        end_time = time.perf_counter()
    elif algorithm == "EA":
        if seed is None:
            print("Error: Seed is required for evolutionary algorithm.")
            return
        start_time = time.perf_counter()
        tour, best_quality = evolutionary_algorithm(filename, cutoff, seed)
        end_time = time.perf_counter()
    else:
        print(f"Error: Unknown algorithm {algorithm}.")
        return
    
    # calculate running time
    elapsed_time = end_time - start_time
    print(algorithm + " on " + filename + f" running time: {elapsed_time:.6f} seconds")

    # generate names of output files
    base_name = os.path.basename(filename).split('.')[0]
    if algorithm == "Approx":
        output_filename = f"{base_name}_{algorithm}"
    else:
        output_filename = f"{base_name}_{algorithm}_{cutoff}"
    if seed is not None and algorithm == "LS":
        output_filename += f"_{seed}"
    if seed is not None and algorithm == "EA":
        output_filename += f"_{seed}"
    output_filename += ".sol"

    # write output files
    with open(output_filename, "w") as f:
        f.write(f"{best_quality}\n")
        f.write(",".join(map(str, tour)) + "\n")
    print(f"Output written to {output_filename}")

if __name__ == "__main__":
    multiprocessing.freeze_support()
    main()