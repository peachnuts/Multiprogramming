# ======================================================================
# MIT License
#
# Copyright (c) [2020] [LIRMM]
# Contributor: Siyuan Niu (<siyuan.niu@lirmm.fr>)
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
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
    hardware = IBMQHardwareArchitecture("ibmq_manhattan")
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



    circuits.append(circuit1)
    circuits.append(circuit2)
    circuits.append(circuit3)
    circuits.append(circuit4)
    circuits.append(circuit5)


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
