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

import numpy as np
from qiskit import *
from qiskit import IBMQ
from pathlib import Path
import json

import qiskit.ignis.verification.randomized_benchmarking as rb
from qiskit.providers.aer.noise import NoiseModel
from qiskit.providers.aer.noise.errors.standard_errors import depolarizing_error, thermal_relaxation_error

IBMQ.load_account()

provider = IBMQ.get_provider(hub='ibm-q-france', group='univ-montpellier', project='reservations')

basis_gates = ['u1','u2','u3','cx']
shots = 1024

backend_name = "ibmq_toronto"
hardware_neighbour_file = str(Path(__file__).parent) + "/hardware_neighbour_pairs/" + backend_name + "_neighbour.json"

backend = provider.get_backend(backend_name)
# Number of qubits
nQ = 27
# Number of seeds
nseeds = 5
# Number of Cliffords in the sequence
nCliffs = np.arange(1, 100, 20)

print(hardware_neighbour_file)
with open(hardware_neighbour_file, 'r') as f:
     hardware_neighbour_date = json.load(f)


output_context = []

for bin in hardware_neighbour_date:
    qbits_pair_1 = []
    qbits_pair_2 = []
    for pairs in bin:
        qbits_pair_1 += [pairs[0]]
        qbits_pair_2 += [pairs[1]]
    print(qbits_pair_1)
    print(qbits_pair_2)

    print(qbits_pair_1, "---", qbits_pair_2)

    for kk in range(3):
        if kk == 0:
            nQ = 2 * len(qbits_pair_1)
            rb_pattern = qbits_pair_1
            length_multiplier = [1] * len(qbits_pair_1)
        elif kk == 1:
            nQ = 2 * (len(qbits_pair_1) + len(qbits_pair_2))
            rb_pattern = qbits_pair_1 + qbits_pair_2
            length_multiplier = [1] * (len(qbits_pair_1) + len(qbits_pair_2))
        elif kk == 2:
            nQ = 2 * len(qbits_pair_2)
            rb_pattern = qbits_pair_2
            length_multiplier = [1] * len(qbits_pair_2)


        rb_opts = {}
        rb_opts['length_vector'] = nCliffs
        rb_opts['nseeds'] = nseeds
        rb_opts['rb_pattern'] = rb_pattern
        rb_opts['length_multiplier'] = length_multiplier
        rb_circs, xdata = rb.randomized_benchmarking_seq(**rb_opts)

        result_list = []
        transpile_list = []

        for rb_seed, rb_circ_seed in enumerate(rb_circs):
            print('Compiling seed %d' % rb_seed)
            rb_circ_transpile = qiskit.transpile(rb_circ_seed, basis_gates=basis_gates)
            print('Simulating seed %d' % rb_seed)
            job = qiskit.execute(rb_circ_transpile, shots=shots, backend=backend)
            result_list.append(job.result())
            transpile_list.append(rb_circ_transpile)

        # gates_per_cliff = rb.rb_utils.gates_per_clifford(transpile_list, xdata[0], basis_gates, rb_opts['rb_pattern'][0])

        rbfit = rb.fitters.RBFitter(result_list, xdata, rb_opts['rb_pattern'])

        if kk == 0:
            output_context.append("-----independent(k=0)--------\n")
        elif kk == 1:
            output_context.append("----with neighbour effect(k=1)-----\n")
        else:
            output_context.append("-----independent(k=2)---------\n")

        for index, pattern in enumerate(rb_pattern):
            epc_0 = rbfit._fit[index]['epc']
            epg_0 = epc_0 / 1.5

            print("Error per Clifford = {}".format(epc_0))
            print("Error per Gate for CX {} = {}".format(pattern, epg_0))
            line = "Error per Gate for CX {} = {}\n".format(pattern, epg_0)
            output_context.append(line)


output_name = backend_name + "_crosstalk.txt"
output = str(Path(__file__).parent) + "/SRB_crosstalk_hardware/" + output_name

with open(output, "w") as file_object:
    for line_idx in output_context:
        file_object.write(line_idx)


