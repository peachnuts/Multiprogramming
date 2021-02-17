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

from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister, QiskitError
import numpy as np
from hardware.IBMQHardwareArchitecture import IBMQHardwareArchitecture
import networkx as nx
from hardware.distance_matrx import (
    get_distance_matrix_cnot_error_cost,
    get_qubit_readout_error,
)
from partition_process.qubit_partition import partition_hardware_heuristic, partition_hardware


hardware = IBMQHardwareArchitecture("ibmq_manhattan")
hardware_graph = nx.DiGraph(hardware._coupling_graph)
cnot_error_matrix = get_distance_matrix_cnot_error_cost(hardware)
readout_error = get_qubit_readout_error(hardware)
from mapping.mapping_transition import circuits_schedule

def get_ucc_ansatz(theta):
    """
    Define the UCCSD ansatz circuit
    """
    circuit = QuantumCircuit(65,2)
    circuit.x(0)
    circuit.ry(theta, 1)
    circuit.cx(1, 0)
    return circuit


def measure_zi_and_iz(theta):
    circuit = get_ucc_ansatz(theta)
    return circuit


def measure_xx_and_yy(theta):
    circuit = get_ucc_ansatz(theta)
    circuit.cx(0, 1)
    circuit.h(0)
    return circuit

def measure_circuit_list(start, end, step):
    circuit_zi_iz = []
    circuit_xx_yy = []
    for theta in np.arange(start, end, step):
        circuit1 = measure_zi_and_iz(theta)
        circuit2 = measure_xx_and_yy(theta)
        circuit_zi_iz.append(circuit1)
        circuit_xx_yy.append(circuit2)
        print(theta)

    return circuit_zi_iz + circuit_xx_yy

epslon = 0.1
weight_lambda = 1


circuit_list = measure_circuit_list(0, np.pi, np.pi / 5)

ansatz_parameter = [5.9, 0.22, -6.1, -2.14, -2.14]

circuits_schedule(circuit_list,
                  hardware_graph,
                  hardware,
                  cnot_error_matrix,
                  readout_error,
                  partition_hardware_heuristic,
                  epslon,
                  weight_lambda,
                  ansatz_parameter,
                  )