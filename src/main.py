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
import networkx as nx
from hardware.distance_matrx import (
    get_distance_matrix_cnot_error_cost,
    get_qubit_readout_error,
)
from tools.read_circuit import read_benchmark_circuit
from mapping.mapping_transition import circuits_schedule
from partition_process.qubit_partition import partition_hardware_heuristic, partition_hardware
if __name__ == '__main__':
    hardware = IBMQHardwareArchitecture("ibmq_toronto")
    hardware_graph = nx.DiGraph(hardware._coupling_graph)
    cnot_error_matrix = get_distance_matrix_cnot_error_cost(hardware)
    readout_error = get_qubit_readout_error(hardware)
    #(pair):(influenced pair)
    # crosstalk_properties = {(2,3):{(5,8): 0.0317},
    #                         (5,8):{(2,3): 0.031},
    #                         (7,10):{(12,15): 0.0242},
    #                         (12,15):{(7,10): 0.0417},
    #                         (15,18):{(10,12): 0.0281},
    #                         }

    # A list of circuit that are supposed to be executed simultaneously
    circuits = []
    circuit1 = read_benchmark_circuit("3_17_13")
    circuit1.name = "3_17_13"
    circuit2 = read_benchmark_circuit("4mod5-v1_22")
    circuit2.name = "4mod5-v1_22"
    circuit3 = read_benchmark_circuit("mod5mils_65")
    circuit3.name = "mod5mils_65"
    circuit4 = read_benchmark_circuit("alu-v0_27")
    circuit4.name = "alu-v0_27"
    circuit5 = read_benchmark_circuit("decod24-v2_43")
    circuit5.name = "decod24-v2_43"
    circuit6 = read_benchmark_circuit("3_17_13")
    circuit6.name = "3_17_13_2"



    circuits.append(circuit1)
    circuits.append(circuit2)
    # circuits.append(circuit3)
    # circuits.append(circuit4)
    # circuits.append(circuit5)
    # circuits.append(circuit6)

    # sort circuit according to ascending order of CNOT density
    circuits = sorted(circuits, key=lambda x: x.count_ops().get("cx", 0) / x.cregs[0].size)
    print([circuit.name for circuit in circuits])

    # threshold of the partition fidelity difference when partitioning independently and simultaneously
    epslon = 0.1
    # weight parameter lambda set by user to weight between the CNOT error rate and readout error rate
    # for fidelity degree of the qubit
    weight_lambda = 2

    # QHSP: partition_hardware_heuristic
    # GSP: partition_hardware
    # Result includes PHA, QHSP/GSP, HA
    circuits_schedule(circuits,
                      hardware_graph,
                      hardware,
                      cnot_error_matrix,
                      readout_error,
                      partition_hardware_heuristic,
                      epslon,
                      weight_lambda,
                      )
