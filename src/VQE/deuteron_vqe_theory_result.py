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
