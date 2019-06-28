'''
This file contains code to generate the Abrams-discretized model D_{n,k}Y, a model of the non-k-equal configuration space 
of the graph that looks like a Y. Most of the code is identical to code published in the appendix of Safia Chettih's PHD dissertation 
(Dancing In The Stars: Topology of Non-k-equal Configuration Spaces of Graphs), but some modifications were made to fix an issue 
where certain cells were not being counted because of additional movement "choices" that arise when k > 2.
Notably, the functions downstream_moves, simultaneous_moves, and downstream_cubes have been modified slightly, and new functions 
iterate_over_count and capacity have been added to perform deeper computations "in between" downstream_moves and simultaneous_moves.

TODO: Variables in the modified/new functions are likely to need more descriptive/meaningful names.
TODO: No attempts have been made to optimize this code for performance, and the reward for doing so would mean the ability to 
      generate complex for larger values of n (and k). It would be useful to look at abrams_y.py, which contains code to 
      generate the Abrams-discretized model D_{n,2}Y (non-2-equal only) that has some performance improvements implemented by 
      Langston Barrett (see https://github.com/siddharthist/computational-homology).
'''
from sage.all import *
import cubical_complex
import itertools
import logging
import collections

def generate_tree(n):
    '''
    We think of the Y-graph with numbered vertices, beginning from 0. The following are illustrations (from left to right)
    of n=1, n=2, and n=3:
     |         2              4
     0__       |              |
    /          1__0           3
              /               |
             3                2__1__0
                             /
                            5
                           /
                          6
    EXAMPLE (n=3):
    >>> generate_tree(n)
    ([(0, 2, 0), (0, 1, 0), (0, 0, 0), (1, 0, 0), (2, 0, 0), (0, 0, 1), (0, 0, 2)], ([1], [2], [3, 5], [4], [], [6], []))

    :param n: An integer representing the number of moving points (or "robots") living in our graph Y.
    :type n: int
    :return: A tuple, where the first element, lookup, is a list of points in 3-space coordinates, and the second element, tuple(tree), lists the possible "downstream" moves from a point to one greater than it.
    :rtype: tuple
    '''
    lookup = [(0, n-t-1, 0) for t in range(n-1)]
    lookup.append((0,0,0))
    lookup.extend([(t, 0, 0) for t in range (1, n)])
    lookup.extend([(0, 0, t) for t in range (1, n)])
    tree = []
    for point in range (n-1):
        tree.append([point+1])
    tree.append([n, 2*n-1])
    for point in range(n, 2*n-2):
        tree.append([point+1])
    tree.append([])
    for point in range(2*n-1, 3*n-3):
        tree.append([point+1])
    tree.append([])
    return (lookup, tuple(tree))

def no_k_equal(point_config, k):
    '''
    A helper function that checks if a particular point configuration is "allowed", i.e. is non-k-equal.
    
    EXAMPLE:
    In the point configuration [5,2,2], robot #0 is at vertex 5, while robots #1 and #2 are at vertex 2:

    >>> no_k_equal(config,3)
    True
    >>> no_k_equal(config,2)
    False

    :param point_config: A list representing how the robots are configured.
    :type point_config: list
    :param k: The integer k for the non-k-equal configuration space.
    :type k: int
    :return: True if no k points are equal in this point configuration, and False otherwise.
    :rtype: bool
    '''
    count_list = [point_config.count(i) < k for i in point_config]
    return count_list.count(False) == 0

def iterate_over_conf(I, n, k):
    '''
    Enumerate all the possible configurations of points at vertices, which will
    be the 0-cells. Gives lists of length n which consist of the positions
    (vertex labels) of different points.
    
    :param I: A tuple of possible "downstream" moves. This is obtained from the second element of the tuple returned by generate_tree.
    :type I: tuple
    :param n: An integer representing the number of moving points (or "robots") living in our graph Y.
    :type n: int
    :param k: The integer k for the non-k-equal configuration space.
    :type k: int
    :return: A generator which can be accessed by enumeration to give all the possible configurations for this space.
    :rtype: generator
    '''
    iterator_list = [xrange(len(I))] * n
    for point_config in  xmrange_iter(iterator_list):
        if no_k_equal(point_config, k): # no k points are equal
            yield point_config

def downstream_moves(point_config, I, k):
    '''
    Returns a list that encodes downstream movement information in this graph given a particular point configuration. 
    
    EXAMPLE:
    Imagine we are again dealing with D_{3,3}Y at the point configuration [5,2,2] where robot #0 is at vertex 5, 
    while robots #1 and #2 are at vertex 2:

    >>> (lookup, I) = generate_tree(3)
    >>> downstream_moves([5,2,2], I, 3)
    [None, None, {3: 2, 5: 1}, None, None, {6: 2}, None]

    (TODO: Is the dictionary at index 5 wrong? The point configuration [5,2,2] means that there is one robot at vertex 5, and so
    shouldn't the dictionary be {6:1}?)

    The element of the list at index i corresponds to points at vertex i, and so if no point is there then the value is None.
    Otherwise the list returns a dictionary where the key:value relation is 
    (point to move to):(number of points that are allowed to move there).
    It is useful here to think of a potential point downstream as a "bucket" of capacity (k-1). Thus, the dictionary above
    at position 2 tell us that from vertex 2, we could move both robots to vertex 3, but only one of them to vertex 5, since the 
    "bucket" at vertex 5 already has a robot living there in this point configuration.

    :param point_config: A list representing how the robots are configured.
    :type point_config: list
    :param I: A tuple of possible "downstream" moves. This is obtained from the second element of the tuple returned by generate_tree.
    :type I: tuple
    :param k: The integer k for the non-k-equal configuration space.
    :type k: int
    :return: A list the length of the number of vertices in the graph, with each element either None or a Dictionary representing moves.
    :rtype: list
    '''
    locationlist = config_to_locationlist(point_config, I)
    output = [None]*len(locationlist)
    multiplicity = [len(u) for u in locationlist]
    for p in set(point_config):
        available_points = multiplicity[p]
        downstream_points = [u for u in I[p] if multiplicity[u] < k-1]
        downstream = {}
        for point in downstream_points:
            vacant_spots = k - multiplicity[point] - 1
            downstream[point] = vacant_spots
        if downstream == {}:
            pass
        else:
            output[p] = downstream
    return output

def config_to_locationlist(config, I):
    '''
    :param config: A list representing how the robots are configured.
    :type config: list
    :param I: A tuple of possible "downstream" moves. This is obtained from the second element of the tuple returned by generate_tree.
    :type I: tuple
    :return: A "location list" with the ith element representing the robots located at vertex i.
    :rtype: list
    '''
    locationlist = []
    for vertex in range(len(I)):
        locationlist.append([])
    for point, location in enumerate(config):
        locationlist[location].append(point)
    return locationlist

def locationlist_to_config(mult, I, n):
    '''
    :param mult: A "location list" with the ith element representing the robots located at vertex i.
    :type mult: list
    :param I: A tuple of possible "downstream" moves. This is obtained from the second element of the tuple returned by generate_tree.
    :type I: tuple
    :param n: An integer representing the number of moving points (or "robots") living in our graph Y.
    :type n: int
    :return: A list representing how the robots are configured.
    :rtype: list
    '''
    config = [None] * n
    for location, pointlist in enumerate(mult):
        for point in pointlist:
            config[point] = location
    return config

def iterate_over_count(count, key_length):
    '''
    TODO: The parameter names are so given based on this function's usage in 'capacity', but could potentially be renamed for easier understanding.

    :param key_length: 
    :type key_length: int
    :param count: 
    :type count: int
    :return: A generator which can be accessed by enumeration to give the possible ways to arrange 'count' robots between 'key_length' vertices.
    :rtype: generator
    '''
    iterator_list = [xrange(count+1)] * key_length
    for partition in  xmrange_iter(iterator_list):
        if sum(partition) == count:
            yield partition

def capacity(point_config, I, k):
    '''
    TODO: This function probably needs renaming.

    A buffer between downstream_moves and simultaneous_moves, that can be thought of a "more detailed" version of the output of downstream_moves.
    EXAMPLE:
    Consider again D_{3,3}Y at the point configuration [5,2,2] where robot #0 is at vertex 5, 
    while robots #1 and #2 are at vertex 2:
    
    >>> (lookup, I) = generate_tree(3)
    >>> capacity([5,2,2], I, 3)
    [None, None, [{3: 1, 5: 1}, {3: 2}], None, None, [{6: 1}], None]

    The return value is a list that describes the capacities of vertexes downstream of the vertex i. For example, at vertex 2 there are
    2 robots, and we can either move we can either move one to vertex 3 and one to vertex 5 or move both robots to vertex 3. We can't move
    both robots to vertex 5 since there is already one living there in this configuration.

    :param point_config: A list representing how the robots are configured.
    :type point_config: list
    :param I: A tuple of possible "downstream" moves. This is obtained from the second element of the tuple returned by generate_tree.
    :type I: tuple
    :param k: The integer k for the non-k-equal configuration space.
    :type k: int
    :return: A list embedding the "bucket capacities" of downstream vertices.
    :rtype: list
    '''
    down_moves = downstream_moves(point_config, I, k)
    capacity_moves = [None for i in range(len(I))]
    for i in range(len(down_moves)):
        if down_moves[i] == None:
            pass
        else:
            capacity_dict = down_moves[i]
            key_list = capacity_dict.keys()
            moves_for_location = []
            value_sum = sum(capacity_dict.values())
            if value_sum < point_config.count(i):
                capacity_moves[i] = [capacity_dict]
            else:
                for partition in iterate_over_count(point_config.count(i), len(capacity_dict)):
                    assert(len(partition) == len(key_list))
                    move_allowed = True
                    for j in range(len(partition)):
                        if partition[j] > capacity_dict[key_list[j]]:
                            move_allowed = False
                            break
                    if move_allowed:
                        moves_for_partition = {}
                        for h in range(len(partition)):
                            if partition[h] != 0:
                                moves_for_partition[key_list[h]] = partition[h]
                        moves_for_location.append(moves_for_partition)
                if len(moves_for_location) > 0:
                    capacity_moves[i] = moves_for_location
    return capacity_moves    

def simultaneous_moves(point_config, I, k):
    '''
    :param point_config: A list representing how the robots are configured.
    :type point_config: list
    :param I: A tuple of possible "downstream" moves. This is obtained from the second element of the tuple returned by generate_tree.
    :type I: tuple
    :param k: The integer k for the non-k-equal configuration space.
    :type k: int
    :return: 
    :rtype: 
    '''
    down_moves = capacity(point_config, I, k)
    locationlist = config_to_locationlist(point_config, I)
    moves = []
    for start in range(len(I)):
        possible_moves = []
        if down_moves[start] == None:
            pass
        else:
            for perm_moves_dict in down_moves[start]:
                perm_moves_dict_keys = perm_moves_dict.keys()
                value_sum = sum(perm_moves_dict.values())
                for t in Permutations(locationlist[start], value_sum).list():
                    offset = 0
                    sub_perm_move = []
                    for i in range(len(perm_moves_dict_keys)):
                        value = perm_moves_dict[perm_moves_dict_keys[i]]
                        sub_perm_move.append([t[offset:value+offset], perm_moves_dict_keys[i]])
                        offset += value
                    possible_moves.append(sub_perm_move)
            moves.append(possible_moves)
    sim_moves = [i for i in xmrange_iter(moves)]

    # TODO: The triple for-loop below is a REALLY silly-looking step, but it is a quick fix to ensure that the output of this function has 
    # all its simultaneous moves together on the same "level". In sim_moves the moves are "grouped together" according to the vertices
    # that the robots currently occupy at a point, so the triple for-loop goes down to that level and removes the grouping.
    # EXAMPLE:
    # >>> sim_moves
    # [[[[[1], 3], [[2], 5]], [[[0], 6]]], [[[[2], 3], [[1], 5]], [[[0], 6]]], [[[[1, 2], 3]], [[[0], 6]]], [[[[2, 1], 3]], [[[0], 6]]]]
    # >>> new_sim_moves
    # [[[[1], 3], [[2], 5], [[0], 6]], [[[2], 3], [[1], 5], [[0], 6]], [[[1, 2], 3], [[0], 6]], [[[2, 1], 3], [[0], 6]]]

    # Notice the extra 'grouping' in sim_moves. Essentially both these lists encode the same information, but new_sim_moves has 
    # a slightly different format, so that it can be accepted by downstream_cubes to correctly build the cubical complex.
    # The 'TODO' here is to try and find a more elegant way to do this so that we don't need require the triple for-loop.
    new_sim_moves = [[] for i in range(len(sim_moves))]
    for i in range(len(sim_moves)):
        for j in range(len(sim_moves[i])):
                for h in range(len(sim_moves[i][j])):
                    new_sim_moves[i].append(sim_moves[i][j][h])
    return new_sim_moves

def downstream_cubes(point_config, I, lookup, k):
    '''
    Builds the highest-dimensional cubes (in the Abrams-discretized configuration space) that result from performing moves at the same time.

    EXAMPLE:
    >>> (lookup, I) = generate_tree(3)
    >>> downstream_cubes([5,2,2], I, lookup, 3)
    [[[0, 0], [0, 0], [1, 2], [0, 1], [0, 0], [0, 0], [0, 0], [0, 0], [0, 1]], [[0, 0], [0, 0], [1, 2], [0, 0], [0, 0], [0, 1], [0, 1], [0, 0], [0, 0]], [[0, 0], [0, 0], [1, 2], [0, 1], [0, 0], [0, 0], [0, 1], [0, 0], [0, 0]]]

    :param point_config: A list representing how the robots are configured.
    :type point_config: list
    :param I: A tuple of possible "downstream" moves. This is obtained from the second element of the tuple returned by generate_tree.
    :type I: tuple
    :param lookup: A list of points in 3-space coordinates.
    :type lookup: list
    :param k: The integer k for the non-k-equal configuration space.
    :type k: int
    :return: A list of maximal cubes resulting from simultaneous moves.
    :rtype: list
    '''
    cubes = []
    sim_moves = simultaneous_moves(point_config, I, k)
    for move in sim_moves:
        coord_list = [None] * len(point_config)
        for [points, place] in move:
            for point in points:
                coord_list[point] = place
        new_cube = []
        for (coordinate, point) in zip(coord_list, point_config):
            embedded_coords = lookup[point]
            if coordinate == None: # This point doesn't move in this cube.
                intervals = [[u, u] for u in embedded_coords]
            else:
                downstream_embedded = lookup[coordinate]
                intervals = [sorted([u, v]) for (u, v) in zip(embedded_coords, downstream_embedded)]
            new_cube.extend(intervals)

        #TODO: Is this if-statement necessary? Sage's CubicalComplex doesn't care if we have duplicate cubes, but perhaps 
        # this improves performance. Leave it in if so, but take it out if not. 
        if new_cube not in cubes:
            cubes.append(new_cube)
    return cubes

def the_complex(n, k):
    '''
    :param n: An integer representing the number of moving points (or "robots") living in our graph Y.
    :type n: int
    :param k: The integer k for the non-k-equal configuration space.
    :type k: int
    :return: The cubical complex D_{n,k}Y.
    :rtype: CubicalComplex
    '''
    (lookup, I) = generate_tree(n)
    cubes = []
    for point_config in iterate_over_conf(I, n, k):
        cubes.extend(downstream_cubes(point_config, I, lookup, k))
    return cubical_complex.CubicalComplex(cubes)
