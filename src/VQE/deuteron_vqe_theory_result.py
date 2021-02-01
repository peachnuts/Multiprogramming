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

from qiskit import QuantumCircuit, execute, Aer
import numpy as np
from tools.submit import zi_iz_measure_contribution, xx_yy_measure_contribution

def get_ucc_ansatz(theta):
    """
    Define the UCCSD ansatz circuit
    """
    circuit = QuantumCircuit(27,2)
    circuit.x(0)
    circuit.ry(theta, 1)
    circuit.cx(1, 0)
    return circuit


def measure_zi_and_iz(theta):
    circuit = get_ucc_ansatz(theta)
    circuit.measure(range(0, 2), range(0, 2))
    return circuit


def measure_xx_and_yy(theta):
    circuit = get_ucc_ansatz(theta)
    circuit.cx(0, 1)
    circuit.h(0)
    circuit.measure(range(0, 2), range(0, 2))
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

circuit_list = measure_circuit_list(0, np.pi, np.pi / 5)
backend = Aer.get_backend('qasm_simulator')

theory_results = []
ansatz_parameter = [5.9, 0.22, -6.1, -2.14, -2.14]

for num in range(len(circuit_list) // 2):
    count1 = execute(circuit_list[num], backend=backend, shots=8192).result().get_counts()
    zi_contribution, iz_contribution = zi_iz_measure_contribution(count1, 8192)
    count2 = execute(circuit_list[num + len(circuit_list) // 2], backend=backend, shots=8192).result().get_counts()
    xx_contribution, yy_contribution = xx_yy_measure_contribution(count2, 8192)
    result_simulator = ansatz_parameter[0] + ansatz_parameter[1] * zi_contribution + ansatz_parameter[
        2] * iz_contribution + \
                       ansatz_parameter[3] * xx_contribution + ansatz_parameter[4] * yy_contribution

    theory_results.append(result_simulator)

print("Theory result:", theory_results)
