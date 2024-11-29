#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
File: utility.py
Description:
    This script provides two utility functions.
        i) read_tsp: Read node coordinates and store the read information in a list.
        ii) construct_adj_matrix: Construct an adjacency matrix based on the result of read_tsp.
    This script also provides packages used by algorithm-related scripts.

Author: Kaichen Shen
Email: kshen77@gatech.edu
"""


import math
import itertools
import numpy as np
import multiprocessing
import copy
import random


def read_tsp(file_path):
    """
    Reads a TSP file and extracts node coordinates.

    Parameters:
        file_path (str): The path to the TSP file.

    Returns:
        coords: A list of tuples, each containing the node ID and its x, y coordinates. Example: [(1, x1, y1), (2, x2, y2), ...].
    
    Raises:
        ValueError: If the 'NODE_COORD_SECTION' symbol is not found in the TSP file.
    """
    
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # find first coordinate by locating "NODE_COORD_SECTION"
    node_coord_start = None
    for i, line in enumerate(lines):
        if line.strip() == 'NODE_COORD_SECTION':
            node_coord_start = i + 1
            break
    if node_coord_start is None:
        raise ValueError("NODE_COORD_SECTION symbol is not found in given tsp file")

    # record coordinates in a list
    coords = []
    for line in lines[node_coord_start:]:
        if line.strip() == 'EOF':
            break
        parts = line.split()
        node_id = int(parts[0])
        x = float(parts[1])
        y = float(parts[2])
        coords.append((node_id, x, y))

    return coords


def construct_adj_matrix(coords):
    """
    Use node coordinates to construct adjacency matrix.

    Parameters:
        coords: A list of tuples, each containing the node ID and its x, y coordinates. Same as function read_tsp's output.

    Returns:
        adj_matrix: A list of lists, adj_matrix[i][j] = adj_matrix[j][i] = Euclidean distance between node (i + 1) and node (j + 1).
    """

    adj_matrix = []

    for i in range(len(coords)):
        adj_matrix.append([])
        for j in range(len(coords)):
            if i < j:
                coord_1_x, coord_1_y = coords[i][1], coords[i][2]
                coord_2_x, coord_2_y = coords[j][1], coords[j][2]
                distance = math.sqrt((coord_1_x - coord_2_x) ** 2 + (coord_1_y - coord_2_y) ** 2)
                adj_matrix[i].append(round(distance))
            elif i == j:
                adj_matrix[i].append(0)
            elif i > j:
                adj_matrix[i].append(adj_matrix[j][i])

    return adj_matrix