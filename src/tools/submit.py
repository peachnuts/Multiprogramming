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
from qiskit import Aer, IBMQ
from qiskit.providers.jobstatus import JobStatus
from tools.IBMQSubmitter import IBMQSubmitter

import logging
logger = logging.getLogger("submit")

def get_jobs(backend, tags):
    jobs = backend.jobs(
        job_tags=tags, job_tags_operator="AND", limit=1, descending=False
    )
    return jobs

def submit_circuits(hardware: IBMQHardwareArchitecture,
                   initial_layouts: ty.List,
                   final_circuits: ty.List,
                   partitions: ty.List,):

    print("Loading account...")
    IBMQ.load_account()
    provider = IBMQ.get_provider(
        hub="ibm-q-france", group="univ-montpellier", project="default"
    )

    backend = provider.get_backend(f"{hardware.name}")
    print(f"Running on {backend.name()}.")
    #backend = provider.get_backend("ibmq_qasm_simulator")

    circuit_tag = ["2circuitss"]
    submitter = IBMQSubmitter(backend, tags=circuit_tag)

    fidelities = []

    for index, circuit in enumerate(final_circuits):
        circuit.measure_active()
        submitter.add_circuit(
            circuit, initial_layouts[index], backend
        )
    print(f"Submitting {len(submitter)} circuits...")
    submitter.submit()

    jobs = get_jobs(backend, circuit_tag)
    for i, job in enumerate(jobs):
        if job.status() == JobStatus.ERROR:
            logger.error("job error")
            exit(1)
        qasm_simulator = Aer.get_backend('qasm_simulator')
        qobj = job.qobj()
        ideal_results = qasm_simulator.run(qobj).result()
        ideal_counts = ideal_results.get_counts()
        hardware_results = job.result()
        hardware_counts = hardware_results.get_counts()
        circuit_number = len(qobj.experiments)
        for num in range(circuit_number):
            ideal_result = list(ideal_counts[num].keys())[0]
            if partitions[num] is None or not isinstance(partitions[num], list):
                fidelities.append(
                    hardware_counts[num].get(ideal_result, 0) / ideal_counts[num][ideal_result]
                )
            else:
                fidelity_merge = []
                start = len(ideal_result)
                for partition in (partitions[num]):
                    start -=  len(partition)
                    end = start + len(partition)
                    fidelity = 0
                    for result, count in hardware_counts[num].items():
                        if result[start:end] == ideal_result[start:end]:
                            fidelity += count
                    fidelity = fidelity / ideal_counts[num][ideal_result]
                    fidelity_merge.append(fidelity)
                fidelities.append(fidelity_merge)
        print(fidelities)

