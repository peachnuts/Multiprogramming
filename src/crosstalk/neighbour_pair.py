# ======================================================================
# Copyright LIRMM (02/2020)
# Contributor: Siyuan Niu (<siyuan.niu@lirmm.fr>)
# This software is governed by the CeCILL-B license under French law and
# abiding  by the  rules of  distribution of free software. You can use,
# modify  and/or  redistribute  the  software  under  the  terms  of the
# CeCILL-B license as circulated by CEA, CNRS and INRIA at the following
# URL "http://www.cecill.info".
#
# As a counterpart to the access to  the source code and rights to copy,
# modify and  redistribute granted  by the  license, users  are provided
# only with a limited warranty and  the software's author, the holder of
# the economic rights,  and the  successive licensors  have only limited
# liability.
#
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using, modifying and/or  developing or reproducing  the
# software by the user in light of its specific status of free software,
# that  may mean  that it  is complicated  to manipulate,  and that also
# therefore  means that  it is reserved for  developers and  experienced
# professionals having in-depth  computer knowledge. Users are therefore
# encouraged  to load and  test  the software's  suitability as  regards
# their  requirements  in  conditions  enabling  the  security  of their
# systems  and/or  data to be  ensured and,  more generally,  to use and
# operate it in the same conditions as regards security.
#
# The fact that you  are presently reading this  means that you have had
# knowledge of the CeCILL-B license and that you accept its terms.
# ======================================================================

import sys
sys.path.append('../')

from hardware.IBMQHardwareArchitecture import IBMQHardwareArchitecture
from hardware.distance_matrx import get_distance_matrix_swap_number
import numpy
import typing as ty
import math
from pathlib import Path
import json

def shortest_path_between_pair(pair1, pair2, distance_matrix):
    """
    Find the shortest path between two pairs
    Each pair has two CNOTs that have neighbour effect (crosstalk)
    """
    shortest_path = min(distance_matrix.item(pair1[0],pair2[0]),
                        distance_matrix.item(pair1[0],pair2[1]),
                        distance_matrix.item(pair1[1],pair2[0]),
                        distance_matrix.item(pair1[1],pair2[1]))
    return shortest_path

def find_crosstalk_pairs(edges: ty.List[ty.List],
                    distance_matrix: numpy.ndarray):
    """
    Find all the CNOT pairs that have neighbour effect (crosstalk)
    """
    crosstalk_pair_list = []
    for first_pair in edges:
        two_pairs = [first_pair]
        for second_pair in edges:
            if [second_pair, first_pair] not in crosstalk_pair_list:
                if shortest_path_between_pair(first_pair, second_pair, distance_matrix) == 1:
                    two_pairs.append(second_pair)
                    crosstalk_pair_list.append(two_pairs)
                    break
    return crosstalk_pair_list


def add_crosstalk_pairs(
        crosstalk_pair: ty.List,
        qubits_in_new_bin: ty.List,
        distance_matrix: numpy.ndarray,
):
    """
    Justify if the pair can be executed in parallel with the pairs in the list
    """
    qubits_used_crosstalk_pair = [qubit for pair in crosstalk_pair for qubit in pair]
    for qubit1 in qubits_used_crosstalk_pair:
        for qubit2 in qubits_in_new_bin:
            if distance_matrix.item(qubit1, qubit2) < 2:
                return False
    return True

def bin_packing(crosstalk_pairs: ty.List,
                distance_matrix: numpy.ndarray):
    """
    Use bin packing algorithm to find the combinations of pairs that can be executed in parallel
    """
    bins = []
    new_bin = []
    crosstalk_pairs_iterate = crosstalk_pairs.copy()
    numpy.random.shuffle(crosstalk_pairs_iterate)
    #print(crosstalk_pairs_iterate)
    qubits_in_new_bin = []
    while crosstalk_pairs_iterate:
        if not new_bin:
            new_bin.append(crosstalk_pairs_iterate[0])
            qubits_in_new_bin.extend([qubit for pair in crosstalk_pairs_iterate[0] for qubit in pair])
            #crosstalk_pairs_iterate.remove(crosstalk_pairs_iterate[0])
        for crosstalk_pair in crosstalk_pairs_iterate:
            if add_crosstalk_pairs(crosstalk_pair, qubits_in_new_bin, distance_matrix) == True:
                new_bin.append(crosstalk_pair)
                qubits_in_new_bin.extend([qubit for pair in crosstalk_pair for qubit in pair])
                #crosstalk_pairs_iterate.remove(crosstalk_pair)
        for pair in new_bin:
            crosstalk_pairs_iterate.remove(pair)
        bins.append(new_bin)
        new_bin = []
        qubits_in_new_bin = []

    return bins


def heuristic_bin_packing(crosstalk_pairs: ty.List,
                          distance_matrix: numpy.ndarray,
                          time: int =100):
    #Use heuristic algorithm to find the better combination
    best_bin_num = math.inf
    best_bins = []
    for i in range(time):
        bins = bin_packing(crosstalk_pairs, distance_matrix)
        if len(bins) < best_bin_num:
            best_bin_num = len(bins)
            best_bins = bins

    return best_bins


def edges_1d(hardware: IBMQHardwareArchitecture):
    """
    Only consider one direction of the CNOT gate
    eg. (0,1) = (1,0)
    """
    edges_all = []
    for i, j in hardware.edges:
        if [i, j] not in edges_all and [j, i] not in edges_all:
            edges_all.append([i, j])
    return edges_all

def test_neighbour_pair():
    hardware_name = "ibmq_toronto"
    hardware = IBMQHardwareArchitecture(hardware_name)
    distance_swap_number = get_distance_matrix_swap_number(hardware)
    edges_all = edges_1d(hardware)
    pairs = find_crosstalk_pairs(edges_all, distance_swap_number)
    bins = heuristic_bin_packing(pairs, distance_swap_number)
    for bin in bins:
        print(bin)
    hardware_neighbour_name = hardware_name + "_neighbour.json"
    output = str(Path(__file__).parent) + "/hardware_neighbour_pairs/" + hardware_neighbour_name
    # with open(output, "w") as file_object:
    #     for bin in bins:
    #         file_object.write(str(bin) + "\n")
    jsObj = json.dumps(bins)
    fileObject = open(output, 'w')
    fileObject.write(jsObj)
    fileObject.close()




test_neighbour_pair()
