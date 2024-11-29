#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
File: approximate_algorithm.py
Description:
    This script implements the evolutionary algorithm solution to TSP.

Author: Kaichen Shen
Email: kshen77@gatech.edu
"""


from utility import read_tsp, construct_adj_matrix, math, np, random, copy, multiprocessing


def create_initial_population(num_nodes, population_size):
    """
    Create initial population for evolutionary algorithm.

    Parameters:
        num_nodes (int): Number of locations in the TSP file.
        population_size (int): Number of tours in tours_population.

    Returns:
        tours_population: A list of tours.
    """

    tours_population = []
    for i in range(population_size):
        tour = random.sample(list(range(num_nodes)), num_nodes)
        tours_population.append(tour)
    return tours_population


def fitness(adj_matrix, tours):
    """
    Compute the fitness of given tour. Fitness is inverse of distance.

    Parameters:
        adj_matrix: A list of lists representing a graph.
        tours: A list of tours.

    Returns:
        tours_fitness: A list of fitness of the tour in tours_population.
    """

    tours_fitness = []
    for tour in tours:
        total_distance = 0
        for i in range(len(tour) - 1):
            total_distance += adj_matrix[tour[i]][tour[i + 1]]
        total_distance += adj_matrix[tour[len(tour) - 1]][tour[0]]
        tour_fitness = 1 / total_distance # fitness is inverse of distance
        tours_fitness.append(tour_fitness)
    return tours_fitness


def normalize(tours_fitness):
    """
    Normalize the fitness of given tours' fitness, with highest fitness equals to 1 and lowest fitness equals to 0.

    Parameters:
        tours_fitness: A list of tours' fitness.

    Returns:
        normalized_tours_fitness: A list of tours' normalized fitness.
    """

    max_fitness = np.max(tours_fitness)
    min_fitness = np.min(tours_fitness)
    max_min_diff = max_fitness - min_fitness
    normalized_tours_fitness = []
    for fitness in tours_fitness:
        normalized_tours_fitness.append((fitness - min_fitness) / max_min_diff)
    return normalized_tours_fitness


def select(tours, tours_fitness):
    """
    Randomly select tours based on their fitness scores. Each fitness score is the selection probability of each tour.

    Parameters:
        tours: A list of tours.
        tours_fitness: A list of fitness of the tour in tours.

    Returns:
        selected_tours: A list, which is a subset of tours.
    """

    selected_tours = []
    for idx, tour in enumerate(tours):
        prob = tours_fitness[idx]
        random_number = random.random()
        if (random_number < prob):
            selected_tours.append(tour)
    return selected_tours


def crossover(tours):
    """
    Perform crossover for pairs of tours.

    Parameters:
        tours: A list of tours.

    Returns:
        new_tours: A list of crossovered tours.
    """

    new_tours = []

    # pair each tour with each neighbour
    tours_pairs = []
    for i in range(len(tours) - 1):
        tours_pair = (tours[i], tours[i + 1])
        tours_pairs.append(tours_pair)

    # each tour pair generated a new tour
    for tours_pair in tours_pairs:
        tour_1 = tours_pair[0]
        tour_2 = tours_pair[1]
        # randomly select a subset of the first parent tour and record the indices
        tour_1_subset_indices = random.sample(list(range(len(tour_1))), math.floor(len(tour_1) * 0.4))
        # generate child tour
        child_tour = [None] * len(tour_2)
        for idx in tour_1_subset_indices:
            child_tour[idx] = tour_1[idx]
        remaining_nodes = set(tour_2) - set(child_tour)
        remaining_nodes = list(remaining_nodes)
        idx = 0
        for i in range(len(child_tour)):
            if child_tour[i] is None:
                child_tour[i] = remaining_nodes[idx]
                idx += 1
        new_tours.append(child_tour)

    return new_tours


def mutate(tours):
    """
    Mutate tours by swapping two nodes with a low probability.

    Parameters:
        tours: A list of tours.

    Returns:
        new_tours: A list of mutated tours.
    """

    new_tours = []

    low_prob = 0.2
    for tour in tours:
        if random.random() < low_prob:
            tour_temp = copy.deepcopy(tour)
            idx1, idx2 = random.sample(range(len(tour_temp)), 2)
            tour_temp[idx2], tour_temp[idx1] = tour_temp[idx1], tour_temp[idx2]
            new_tours.append(tour_temp)
        else:
            new_tours.append(tour)

    return new_tours


def find_shortest_tour(adj_matrix, tours):
    """
    Find a tour of lowest total distance (highest fitness) given a list of tours.

    Parameters:
    adj_matrix: A list of lists representing a graph.
        tours: A list of tours.

    Returns:
        shortest_tour: A tours of lowest total distance (highest fitness).
        shortest_tour_total_distance: Total distance of the shortest tour.
    """

    tours_fitness = fitness(adj_matrix, tours)
    highest_fitness = -np.inf
    shortest_tour = []
    for idx, tour_fitness in enumerate(tours_fitness):
        if tour_fitness > highest_fitness:
            highest_fitness = tour_fitness
            shortest_tour = tours[idx]
    shortest_tour_total_distance = 0
    for i in range(len(shortest_tour) - 1):
        shortest_tour_total_distance += adj_matrix[shortest_tour[i]][shortest_tour[i + 1]]
    shortest_tour_total_distance += adj_matrix[shortest_tour[len(shortest_tour) - 1]][shortest_tour[0]]
    return shortest_tour, shortest_tour_total_distance


def evolutionary_algorithm_worker(file_path, return_dict):
    """
    Evolutionary algorithm to solve TSP.

    Parameters:
        file_path (str): The path to the TSP file.
        return_dict (dict): An empty dictionary.

    Updates:
        (tour, total_distance): tour is a list of node id (the first and last node are different). total_distance is the distance travelled in the tour.
        return_dict["tour"] = tour, return_dict["total_distance"] = total_distance
    """

    coordinates = read_tsp(file_path)
    adj_matrix = construct_adj_matrix(coordinates)

    # initialize a population of tours
    num_nodes = len(adj_matrix)
    tours_population = create_initial_population(num_nodes, 5 * num_nodes)

    # initialize exit condition
    exit_flag = 0
    shortest_tour = []
    #shortest_tour_total_distance = 0

    while (True):
        # compute the fitness of each tour in tours_population
        tours_fitness = fitness(adj_matrix, tours_population)
        # normalize fitness scores, each tour will have fitness between 0 and 1
        tours_fitness = normalize(tours_fitness)

        # randomly select tours based on their fitness scores
        selected_tours = select(tours_population, tours_fitness)

        # perform crossover for pairs of selected tours
        new_tours = crossover(selected_tours)

        # mutate new tours
        new_tours = mutate(new_tours)

        # update population
        tours_population = selected_tours + new_tours

        # record temporary results
        shortest_tour_temp, shortest_tour_total_distance_temp = find_shortest_tour(adj_matrix, tours_population)
        return_dict['tour'] = shortest_tour_temp
        return_dict['total_distance'] = shortest_tour_total_distance_temp

        # exit if tour with shortest distance doesn't change in 5 iterations
        if shortest_tour != shortest_tour_temp:
            exit_flag = 0
            shortest_tour = shortest_tour_temp
            #shortest_tour_total_distance = shortest_tour_total_distance_temp
        else:
            exit_flag += 1
        if exit_flag == 4:
            break
        

def evolutionary_algorithm(file_path, cutoff, seed):
    """
    Use evolutionary algorithm to solve TSP with time cutoff. If the algorithm's running time exceeds cutoff, the algorithm will be terminated and return the solution found so far.

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

    process = multiprocessing.Process(target=evolutionary_algorithm_worker, args=(file_path, return_dict))
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