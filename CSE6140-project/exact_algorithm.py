#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
File: exact_algorithm.py
Description:
    This script implements the brute-force solution to TSP.

Author: Kaichen Shen
Email: kshen77@gatech.edu
"""


from utility import read_tsp, construct_adj_matrix, itertools, np, multiprocessing


def exact_algorithm_worker(file_path, return_dict):
    """
    Brute-force solution to TSP.

    Parameters:
        file_path (str): The path to the TSP file.
        return_dict (dict): An empty dictionary.

    Updates:
        (tour, total_distance): tour is a list of node id (the first and last node are different). total_distance is the distance travelled in the tour.
        return_dict["tour"] = tour, return_dict["total_distance"] = total_distance
    """

    coordinates = read_tsp(file_path)
    adj_matrix = construct_adj_matrix(coordinates)

    # initialize tour and total_distance
    tour = []
    total_distance = np.inf

    # use itertools.permutations to find all possible tours
    nodes = [coord[0] for coord in coordinates]
    all_tours = itertools.permutations(nodes)

    for tour_temp in all_tours:
        tour_temp = list(tour_temp)
        # calculate distance of current tour
        total_distance_temp = 0
        for i in range(len(tour_temp) - 1):
            total_distance_temp += adj_matrix[tour_temp[i] - 1][tour_temp[i + 1] - 1] # "-1" is necessary because adj_matrix[i][j] = distance between node (i + 1) and node (j + 1)
        total_distance_temp += adj_matrix[tour_temp[len(tour_temp) - 1] - 1][tour_temp[0] - 1] # distance from the last node to the first node
        # update if current tour has lower distance
        if total_distance_temp < total_distance:
            total_distance = total_distance_temp
            tour = tour_temp
            # store current best results into return_dict
            return_dict["tour"] = tour
            return_dict["total_distance"] = total_distance


def exact_algorithm(file_path, cutoff):
    """
    Brute-force solution to TSP with time cutoff. If the algorithm's running time exceeds cutoff, the algorithm will be terminated and return the solution found so far.

    Parameters:
        file_path (str): The path to the TSP file.
        cutoff (int): Maximum seconds the algorithm is allowed to execute.

    Updates:
        (tour, total_distance): tour is a list of node id (the first and last node are different). total_distance is the distance travelled in the tour.
    """

    manager = multiprocessing.Manager()
    return_dict = manager.dict()

    process = multiprocessing.Process(target=exact_algorithm_worker, args=(file_path, return_dict))
    process.start()
    process.join(timeout=cutoff)

    if process.is_alive():
        print(f"Execution exceeded cutoff time of {cutoff} seconds.")
        process.terminate()
    return return_dict.get("tour"), return_dict.get("total_distance")


'''
# local test
if __name__ == '__main__':
    # no time out
    exact_Cincinnati = exact_algorithm('DATA/Cincinnati.tsp', 40)
    print(f"Best tour: {exact_Cincinnati[0]}")
    print(f"Total distance: {exact_Cincinnati[1]}")

    # time out
    exact_NYC = exact_algorithm('DATA/NYC.tsp', 40)
    print(f"Best tour: {exact_NYC[0]}")
    print(f"Total distance: {exact_NYC[1]}")
'''