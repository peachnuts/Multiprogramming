# ======================================================================
# Copyright TOTAL / CERFACS / LIRMM (12/2020)
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
from qiskit import QuantumCircuit
from qiskit.converters import circuit_to_dag, dag_to_circuit
from qiskit.circuit.classicalregister import Clbit
from qiskit.dagcircuit.dagcircuit import DAGCircuit
from qiskit.circuit.quantumregister import Qubit
from mapping.iterative_mapping import iterative_mapping_algorithm
from mapping.initial_mapping_wrapper import initial_mapping
from mapping.initial_mapping_construct import cost
from HA.src.hamap import ha_mapping_paper_compliant
from partition_process.qubit_partition import (
    partition_hardware_heuristic,
    partition_circuits,
    largest_circuit_logical_degree,
    hardware_qubit_physical_degree,
)
import networkx as nx
from tools.submit import submit_circuits

import logging
import typing as ty
import numpy as np
import time

logger = logging.getLogger("mapping_transisiton")

def _modify_dag_circuit(circuit: QuantumCircuit, previous_qubits_used: int):
    """
    Modify the dag circuit in order to merge several dag circuits into one dag circuit
    """
    origin_dag = circuit_to_dag(circuit)
    new_dag = DAGCircuit()
    for qreg in origin_dag.qregs.values():
        new_dag.add_qreg(qreg)
    for creg in origin_dag.cregs.values():
        new_dag.add_creg(creg)
    for node in origin_dag.topological_op_nodes():
        new_dag.apply_operation_back(node.op,
                                       qargs=[Qubit(new_dag.qregs['q'], qarg.index + previous_qubits_used) for qarg in
                                              node.qargs],
                                       cargs=[Clbit(new_dag.cregs['c'], carg.index) for carg in
                                              node.cargs])
    new_circuit = dag_to_circuit(new_dag)
    return new_circuit

def multiprogram_initial_mapping(
        circuits: ty.List[QuantumCircuit],
        mappings: ty.List[ty.Dict[Qubit, int]]
) -> ty.Dict[Qubit, int]:
    """
    Contruct a complete initial mapping for multiprogramming process
    """

    multi_initial_mapping = dict()
    partition_qubit_number_diff = len(circuits[0].qubits)
    for index, circuit in enumerate(circuits):
        dag = circuit_to_dag(circuit)
        qubits_non_idle = [qubit for qubit in circuit.qubits if qubit not in dag.idle_wires()]
        partition_qubit_number_diff -= len(qubits_non_idle)
        for qubit in qubits_non_idle:
            multi_initial_mapping[qubit] = mappings[index][qubit]

    left_physical_qubit_list = []

    for i in range(len(circuits[0].qubits)):
        if i not in multi_initial_mapping.values():
            left_physical_qubit_list.append(i)

    j = 0
    for i, qubit in enumerate(circuits[0].qubits):
        if qubit in multi_initial_mapping.keys():
            continue
        else:
            multi_initial_mapping[qubit] = left_physical_qubit_list[j]
            j += 1

    return multi_initial_mapping

def cost_gate_num(quantum_circuit: QuantumCircuit):
    cx_num = quantum_circuit.count_ops().get("cx", 0)
    swap_num = quantum_circuit.count_ops().get("swap", 0) * 3
    ops_num = (
        cx_num + swap_num
    )
    return ops_num

def multiprogram_mapping(circuits: ty.List[QuantumCircuit],
                         hardware: IBMQHardwareArchitecture,
                         circuit_partitions: ty.List,
                         ):
    """
    Perform the qubit mapping algorithm for the multiprogramming mechanism.
    Include initial mapping generation and mapping transition.
    """
    circuit_partitions = [list(part.value) for part in circuit_partitions]
    print(circuit_partitions)

    # obtain the complete initial mapping of the merged circuit
    circuit_initial_mapping = dict()
    computed_initial_mappings = []
    update_circuits = []
    previous_qubit_used = 0

    num_cnots_circuits = sum([cost_gate_num(circuit) for circuit in circuits])

    for index, circuit in enumerate(circuits):
        circuit = _modify_dag_circuit(circuit, previous_qubit_used)
        update_circuits.append(circuit)
        computed_initial_mapping = initial_mapping(
            circuit, hardware, circuit_partitions[index], iterative_mapping_algorithm, cost, "sabre", 10,
            circuit_initial_mapping,
        )
        computed_initial_mappings.append(computed_initial_mapping)
        previous_qubit_used += len(circuit_partitions[index])
    merge_final_mapping = multiprogram_initial_mapping(update_circuits, computed_initial_mappings)

    # the result circuit of the merged circuit
    merge_circuit, merged_final_mapping = iterative_mapping_algorithm(
        update_circuits,
        hardware,
        merge_final_mapping,
        circuit_partitions,
        )
    num_cnots_merge_circuit = cost_gate_num(merge_circuit)
    num_additional_cnots = num_cnots_merge_circuit - num_cnots_circuits

    print(f"additional cnots is {num_additional_cnots}")

    initial_layout = merge_final_mapping.values()

    return initial_layout, merge_circuit


def circuits_schedule(circuits: ty.List[QuantumCircuit],
                      hardware_graph,
                      hardware: IBMQHardwareArchitecture,
                      cnot_error_matrix: np.ndarray,
                      readout_error: ty.List,
                      partition_method: ty.Callable[[
                           IBMQHardwareArchitecture,
                           nx.DiGraph,
                           QuantumCircuit,
                           np.ndarray,
                           ty.List,
                           ty.Set,
                           ty.List,
                           int,
                           ty.Dict,
                       ], ty.List],
                      epslon,
                      weight_lambda,
                      circuit_tag,
                      crosstalk_properties: ty.Dict=None):

    initial_layouts = []
    final_circuits = []
    partitions = []

    # Sort circuit according to ascending order of CNOT density
    circuits = sorted(circuits, key=lambda x: x.count_ops().get("cx", 0) / x.cregs[0].size)

    # Pick up K circuits that are able to be executed on hardware at the same time
    # sum(n_i) <= N (qubit number of hardware), 1 <= i <= K
    circuit_list = []
    qubit_circuit_sum = 0
    for circuit in circuits:
        qubit_circuit_sum += circuit.cregs[0].size
        if qubit_circuit_sum <= hardware.qubit_number:
            circuit_list.append(circuit)
        else:
            break

    qubit_physical_degree, largest_physical_degree = hardware_qubit_physical_degree(hardware)

    largest_logical_degrees = []
    for circuit in circuit_list:
        largest_logical_degrees.append(largest_circuit_logical_degree(circuit))

    # Partition independently (PHA algorithm)
    partition_fidelity_independent_list = []
    independent_partitions = []
    for circuit in circuit_list:
        independent_partition = partition_circuits([circuit],
                                                   hardware_graph,
                                                   hardware,
                                                   cnot_error_matrix,
                                                   readout_error,
                                                   qubit_physical_degree,
                                                   largest_physical_degree,
                                                   largest_logical_degrees,
                                                   weight_lambda,
                                                   partition_method,
                                                   )
        partitions.append(independent_partition[0].value)
        partition_fidelity_independent_list.append(independent_partition[0].fidelity)
        independent_partitions.append(independent_partition)

    # If K > 1, circuits are executed on the hardware simultaneously (Parallelism metric), K = length of circuit list
    while len(circuit_list) > 1:
        # Partition simultaneously (multiprogramming)
        start = time.time()
        multiple_partition = partition_circuits(circuit_list,
                                                hardware_graph,
                                                hardware,
                                                cnot_error_matrix,
                                                readout_error,
                                                qubit_physical_degree,
                                                largest_physical_degree,
                                                largest_logical_degrees,
                                                weight_lambda,
                                                partition_method,
                                                crosstalk_properties,
                                                )
        #print(f"time is {time.time() - start}")

        if not multiple_partition:
            circuit_list.pop()
            largest_logical_degrees.pop()
            partition_fidelity_independent_list.pop()
            partitions.pop()
            continue

        partition_fidelity_multiple = 0.0
        partition_fidelity_independent = sum(partition_fidelity_independent_list)
        for partition in multiple_partition:
            partition_fidelity_multiple += partition.fidelity


        # Post qubit partition process
        partition_fidelity_difference = abs(partition_fidelity_independent - partition_fidelity_multiple)
        #print("paritition fidelity difference is", partition_fidelity_difference)

        if partition_fidelity_difference < epslon:
            print(f"circuits that are executed simultaneously with threshold {epslon}")
            for idx, circuit in enumerate(circuit_list):
                print(f"circuit name : {circuit.name}")
                print("Independent partition (PHA)")
                initial_layout, final_circuit = multiprogram_mapping([circuit], hardware, independent_partitions[idx])
                initial_layouts.append(initial_layout)
                final_circuits.append(final_circuit)

            print("Simultaneous partition")
            initial_layout, final_circuit = multiprogram_mapping(circuit_list, hardware, multiple_partition)
            initial_layouts.append(initial_layout)
            final_circuits.append(final_circuit)
            partitions.append([partition.value for partition in multiple_partition])
            print(f"The number of circuits that are executed on hardware simultaneously is {len(circuit_list)}")
            break

        else:
            circuit_list.pop()
            largest_logical_degrees.pop()
            partition_fidelity_independent_list.pop()
            partitions.pop()
            continue

    # If only one circuit can be executed on the hardware, all the circuits should be executed independently (using HA)
    if len(circuit_list) == 1:
        # HA mapping
        for circuit in circuits:
            circuit_initial_mapping_ha = dict()
            computed_initial_mapping = initial_mapping(
                circuit, hardware, None, ha_mapping_paper_compliant, cost, "sabre", 10, circuit_initial_mapping_ha
            )

            mapped_circuit, final_mapping = ha_mapping_paper_compliant(
                circuit, hardware, computed_initial_mapping,
            )
            num_cnots_circuits = cost_gate_num(circuit)
            num_cnots_merge_circuit = cost_gate_num(mapped_circuit)
            num_additional_cnots = num_cnots_merge_circuit - num_cnots_circuits
            print(f"additional cnots is {num_additional_cnots}")
            initial_layouts.append(computed_initial_mapping.values())
            final_circuits.append(mapped_circuit)
            partitions.append(None)

    #quit()
    submit_circuits(hardware, initial_layouts, final_circuits, partitions, circuit_tag)
