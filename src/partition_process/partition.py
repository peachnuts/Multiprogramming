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
import numpy as np
import typing as ty
from qiskit import QuantumCircuit
from hardware.IBMQHardwareArchitecture import IBMQHardwareArchitecture
from networkx.algorithms.distance_measures import diameter

class Partition:
    def __init__(self, value):
        self.value = value

    def partition_connectivity_error_rate_greedy(self,
                                          hardware: IBMQHardwareArchitecture,
                                          circuit: QuantumCircuit,
                                          cnot_error_matrix: np.ndarray,
                                          readout_list: ty.List,
                                          crosstalk_pairs: ty.Dict = None,
                                          ):
        """
        Calculate the score of each partition candidate obtained by greedy algorithm.
        Score = L(longest_shortest_path) + Avg_CNOT * #CNOTs + Sum(R_Qi)
        :param hardware: hardware target
        :param circuit:  circuit for partition
        :param cnot_error_matrix: cnot error matrix of the hardware
        :param readout_list: list of readout error of physical qubit used in quantum hardware
        :param crosstalk_pairs: CNOT pairs that have strong crosstalk effect
        :return: score of the partition
        """

        longest_shortest_distance = diameter(self.value)
        cnot_error_rate = 0.0
        cnot_pair_num = 0
        for i in self.value:
            for j in self.value:
                crosstalk_effect = 0.0
                if (i, j) in hardware.edges:
                    if crosstalk_pairs and (i, j) in crosstalk_pairs.keys():
                        crosstalk_effect = crosstalk_pairs[(i, j)]
                    cnot_error_rate += cnot_error_matrix.item((i, j)) + crosstalk_effect
                    cnot_pair_num += 1
        cnot_error_average = cnot_error_rate / cnot_pair_num if cnot_error_rate != 0 else 0

        readout_error_rate = 0.0

        for i in self.value:
            readout_error_rate += readout_list[i]
        readout_error_average = readout_error_rate / len(self.value) if len(self.value) != 0 else 0

        circuit_qubit_num = circuit.cregs[0].size
        circuit_cnot_num = circuit.count_ops().get("cx", 0)

        self.longest_path = longest_shortest_distance
        self.cnot_error = cnot_error_average * circuit_cnot_num
        self.readout_error = readout_error_average * circuit_qubit_num

        self.fidelity = longest_shortest_distance + cnot_error_average * circuit_cnot_num + readout_error_average * circuit_qubit_num



    def partition_connectivity_error_rate_heuristic(self,
                                          hardware: IBMQHardwareArchitecture,
                                          circuit: QuantumCircuit,
                                          cnot_error_matrix: np.ndarray,
                                          readout_list: ty.List,
                                          crosstalk_pairs: ty.Dict = None,
                                          ):
        """
        Calculate the score of each partition candidate obtained by heuristic algorithm.
        Score = Avg_CNOT * #CNOTs + Sum(R_Qi)
        :param hardware: hardware target
        :param circuit:  circuit for partition
        :param cnot_error_matrix: cnot error matrix of the hardware
        :param readout_list: list of readout error of physical qubit used in quantum hardware
        :param crosstalk_pairs: CNOT pairs that have strong crosstalk effect
        :return: score of the partition
        """

        cnot_error_rate = 0.0
        cnot_pair_num = 0
        for i in self.value:
            for j in self.value:
                crosstalk_effect = 0.0
                if (i, j) in hardware.edges:
                    if crosstalk_pairs and (i, j) in crosstalk_pairs.keys():
                        crosstalk_effect = crosstalk_pairs[(i, j)]
                    cnot_error_rate += cnot_error_matrix.item((i, j)) + crosstalk_effect
                    cnot_pair_num += 1
        cnot_error_average = cnot_error_rate / cnot_pair_num if cnot_error_rate != 0 else 0

        readout_error_rate = 0.0

        for i in self.value:
            readout_error_rate += readout_list[i]
        readout_error_average = readout_error_rate / len(self.value) if len(self.value) != 0 else 0

        circuit_qubit_num = circuit.cregs[0].size
        circuit_cnot_num = circuit.count_ops().get("cx", 0)

        self.cnot_error = cnot_error_average * circuit_cnot_num
        self.readout_error = readout_error_average * circuit_qubit_num

        self.fidelity = cnot_error_average * circuit_cnot_num + readout_error_average * circuit_qubit_num




