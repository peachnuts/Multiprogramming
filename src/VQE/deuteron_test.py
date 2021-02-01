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