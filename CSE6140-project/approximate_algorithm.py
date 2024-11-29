#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
File: approximate_algorithm.py
Description:
    This script implements the approximate solution to TSP.

Author: Kaichen Shen
Email: kshen77@gatech.edu
"""


from utility import read_tsp, construct_adj_matrix, math, np, copy


def has_cycle(adj_matrix, node_1_index, node_2_index):
    """
    Check whether adding an edge indicated by index will lead to a cycle in the graph represented by adj_matrix.

    Parameters:
        adj_matrix: A list of lists representing a graph.
        node_1_index: A integer representing node (node_1_index + 1).
        node_2_index: A integer representing node (node_2_index + 1).

    Returns:
        True if a cycle is created, False on the other hand.
    """

    num_nodes = len(adj_matrix)
    result = False

    # run bfs to find whether an edge exists between node 1 and node 2 before adding the edge
    frontier = [node_1_index]
    explored = []
    while len(frontier) != 0:
        i = frontier.pop()
        explored.append(i)
        if i == node_2_index:
            result = True
            break
        for j in range(num_nodes):
            if adj_matrix[i][j] != 0 and j not in explored:
                frontier.insert(0, j)

    return result


def MST_Krustal(adj_matrix):
    """
    Krustal algorithm to find one MST of a graph.

    Parameters:
        adj_matrix: A list of lists representing a graph.

    Returns:
        mst: adjacency matrix of the MST
    """

    # initialize mst
    mst = []
    num_nodes = len(adj_matrix)
    for i in range(num_nodes):
        mst.append([0] * num_nodes)

    # copy and modify the given adjacency matrix
    adj_matrix_copy = copy.deepcopy(adj_matrix)
    for i in range(num_nodes):
        adj_matrix_copy[i][i] = np.inf

    # Krustal algorithm
    num_iters = 0
    while (num_iters < num_nodes - 1):
        # indices of nodes corresponde to lowest weight edge
        index = np.argmin(adj_matrix_copy)
        node_1_index = math.floor(index / num_nodes)
        node_2_index = index % num_nodes
        # add the lowest weight to MST unless encounter cycle
        if not has_cycle(mst, node_1_index, node_2_index):
            mst[node_1_index][node_2_index] = adj_matrix_copy[node_1_index][node_2_index]
            mst[node_2_index][node_1_index] = mst[node_1_index][node_2_index]
            num_iters += 1
        # remove considered edge
        adj_matrix_copy[node_1_index][node_2_index] = np.inf
        adj_matrix_copy[node_2_index][node_1_index] = np.inf

    return mst


def approximate_algorithm(file_path):
    """
    Approximate solution to TSP.

    Parameters:
        file_path (str): The path to the TSP file.

    Returns:
        (tour, total_distance): tour is a list of node id (the first and last node are different). total_distance is the distance travelled in the tour.
    """

    coordinates = read_tsp(file_path)
    adj_matrix = construct_adj_matrix(coordinates)

    # initialize tour and total_distance
    tour = []
    total_distance = 0

    # Krustal algorithm to find one MST
    mst = MST_Krustal(adj_matrix)

    # perform a preorder walk of the MST (use dfs)
    frontier = [0]
    num_nodes = len(mst)
    while (len(frontier) != 0):
        node = frontier.pop()
        tour.append(node)
        for j in range(num_nodes):
            if mst[node][j] != 0 and j not in tour:
                frontier.append(j)

    # node indices in tsp files start from 1
    for i in range(len(tour)):
        tour[i] += 1
    # calculate total distance
    for i in range(len(tour) - 1):
        total_distance += adj_matrix[tour[i] - 1][tour[i + 1] - 1]
    total_distance += adj_matrix[tour[len(tour) - 1] - 1][tour[0] - 1]

    return (tour, total_distance)