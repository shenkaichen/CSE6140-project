#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
File: simulated_annealing.py
Description:
    This script implements the simulated annealing solution to TSP.

Author: Kaichen Shen
Email: kshen77@gatech.edu
"""


from utility import read_tsp, construct_adj_matrix, np, copy, random, multiprocessing


def simulated_annealing_worker(file_path, return_dict):
    """
    Use simulated annealing to solve TSP.

    Parameters:
        file_path (str): The path to the TSP file.
        return_dict (dict): An empty dictionary.

    Updates:
        (tour, total_distance): tour is a list of node id (the first and last node are different). total_distance is the distance travelled in the tour.
        return_dict["tour"] = tour, return_dict["total_distance"] = total_distance
    """

    coordinates = read_tsp(file_path)
    adj_matrix = construct_adj_matrix(coordinates)

    # constants
    K = 10                # Boltzmann constant
    T = 100               # initial temperature
    nsteps = 200          # number of iterations
    M = 10                # the temperature will be lowered every M steps
    coolingFraction = 0.7 # fraction by which temperature is lowered

    # initialize tour and total_distance
    num_nodes = len(adj_matrix)
    tour = random.sample(list(range(num_nodes)), num_nodes)
    tour = tour
    total_distance = 0
    for i in range(len(tour) - 1):
        total_distance += adj_matrix[tour[i]][tour[i + 1]]
    total_distance += adj_matrix[tour[len(tour) - 1]][tour[0]]

    # start iteration
    for i in range(1, nsteps):
        # swap node pair at random to get temp tour
        tour_temp = copy.deepcopy(tour)
        idx1, idx2 = random.sample(range(len(tour)), 2)
        tour_temp[idx2], tour_temp[idx1] = tour_temp[idx1], tour_temp[idx2]

        # calculate total distance of temp tour
        total_distance_temp = 0
        for i in range(len(tour_temp) - 1):
            total_distance_temp += adj_matrix[tour_temp[i]][tour_temp[i + 1]]
        total_distance_temp += adj_matrix[tour_temp[len(tour_temp) - 1]][tour_temp[0]]

        # update tour
        if (total_distance_temp < total_distance):
            # update tour directly if tour_temp has lower total distance
            tour = copy.deepcopy(tour_temp)
            return_dict['tour'] = tour
            return_dict['total_distance'] = total_distance
        else:
            # update tour with probability if tour_temp has higher total distance
            delta_total_distance = total_distance_temp - total_distance
            p = np.e ** ((-1 * delta_total_distance) / (K * T))
            if (random.random() < p):
                tour = copy.deepcopy(tour_temp)
                return_dict['tour'] = tour
                return_dict['total_distance'] = total_distance

        # lower the temperature every M steps
        if (i % M == 0):
            T *= coolingFraction


def simulated_annealing(file_path, cutoff, seed):
    """
    Use simulated annealing to solve TSP with time cutoff. If the algorithm's running time exceeds cutoff, the algorithm will be terminated and return the solution found so far.

    Parameters:
        file_path (str): The path to the TSP file.
        cutoff (int): Maximum seconds the algorithm is allowed to execute.
        seed (int): random seed

    Updates:
        (tour, total_distance): tour is a list of node id (the first and last node are different). total_distance is the distance travelled in the tour.
    """

    manager = multiprocessing.Manager()
    return_dict = manager.dict()

    # set random seed
    random.seed(seed)

    process = multiprocessing.Process(target=simulated_annealing_worker, args=(file_path, return_dict))
    process.start()
    process.join(timeout=cutoff)

    if process.is_alive():
        print(f"Execution exceeded cutoff time of {cutoff} seconds.")
        process.terminate()
    
    # Ensure return_dict has valid data
    tour = return_dict.get("tour", [])
    total_distance = return_dict.get("total_distance", float('inf'))

    if tour:
        # Increment node indices if a valid tour is found
        tour = [node + 1 for node in tour]
    else:
        print("Warning: No valid tour found. Returning an empty tour.")

    return tour, total_distance