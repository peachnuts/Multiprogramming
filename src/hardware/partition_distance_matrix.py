# ======================================================================
# Copyright LIRMM (12/2020)
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
from hardware.IBMQHardwareArchitecture import IBMQHardwareArchitecture
import typing as ty
import numpy as np

def adj_matrix_construct(
        hardware: IBMQHardwareArchitecture,
        original_distance_matrix: ty.Callable[
            [IBMQHardwareArchitecture], np.ndarray,
        ],
        partition: ty.List = None,
        ) -> np.ndarray:
    """
    Construct the adjacent matrix according to the given original distance matrix
    :param hardware: selected hardware
    :param partition: allocated partition
    :param original_distance_matrix: distance matrix calculated according to the architecture of the hardware,
           there are different types of distance matrices, such as swap number, error rate or the combination
           of them.
    :return: The constructed adjacent matrix
    """
    if partition == None:
        partition = [i for i in range(hardware.qubit_number)]
    partition_pairs = []
    for i in partition:
        for j in partition:
            if (i,j) in hardware.edges:
                partition_pairs.append((i,j))
    qubit_num = hardware.qubit_number
    adj_mat = np.zeros((qubit_num, qubit_num))

    distance_matrix = original_distance_matrix(hardware)
    for (i,j) in partition_pairs:
        adj_mat[i][j] = distance_matrix.item((i,j))
    return adj_mat


def partition_distance_matrix(qubit_num: int,
          adj_mat: np.ndarray
          ) -> np.ndarray:
    """
    Calculate the distance matrix after the partition by Floyd-Warshall algorithm.
    """
    distance_mat = np.zeros((qubit_num, qubit_num))

    for i in range(qubit_num):
        for j in range(qubit_num):
            if adj_mat.item((i,j)) != 0:
                distance_mat[i][j] = adj_mat.item((i,j))
            else:
                distance_mat[i][j] = 1000000000
        distance_mat[i][i] = 0

    for k in range(qubit_num):
        for i in range(qubit_num):
            for j in range(qubit_num):
                if distance_mat[i][j] > distance_mat[i][k] + distance_mat[k][j]:
                    distance_mat[i][j] = distance_mat[i][k] + distance_mat[k][j]

    return distance_mat

