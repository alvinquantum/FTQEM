#!/usr/bin/env python3
from qiskit import QuantumCircuit
from numpy import pi
from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel, depolarizing_error
import math

def chunkwise(t, size=2):
    '''Rewrites an iterable into a list where each element is an iterable
    of a desired portion of the list. Ex use: [0, 1, 2, 3] -> [[0,1], [2,3]].'''
    return [t[i:i+size] for i in range(0, len(t), size)]

def to_percentages(result):
    total = sum(result.values())
    return {k: v/ total for k, v in result.items()}

def sort_dict(in_dict):
    return {k: in_dict[k] for k in sorted(in_dict.keys())}

def keys_to_binary(in_dict):
    return {bin(int(k, 2)):v for k,v in in_dict.items()}


def run_hardware(circ, backend, shots=10000):
    print(circ)
    bit_flip_code_job_set=backend.run(circ, name="bit_flip_concatenation_cnot", shots=shots, optimization_level=0)
    print(f"job_id: {bit_flip_code_job_set.job_id()}")
    result=bit_flip_code_job_set.result()
    counts=result.get_counts()
    return counts

def get_sso(dist1, dist2):
    '''Returns the square of the statistical overlap. 
    dist1 and dist2 are probability distributions.
    dist1: list
    dits2: list'''
    sum=0
    common_keys=dist1.keys() & dist2.keys()
    # print(f"dist1: {dist1}")
    # print(f"dist2: {dist2}")
    for key in common_keys:
        # print(f"key: {key}")
        # print(f"dist1[key]: {dist1[key]}")
        # print(f"dist2[key]: {dist2[key]}")
        sum+=math.sqrt(dist1[key]*dist2[key])
    return sum**2

def execute_circ_with_depol(circ, prob1, SHOTS=10000):
    # BACKEND=FakeGuadalupe()
    depol_noise_model=depol_noise(prob1)
    # print(depol_noise_model)
    # backend_sim=AerSimulator.from_backend(BACKEND)
    backend_sim=AerSimulator(noise_model=depol_noise_model)
    # grover_qc=transpile(grover_qc, backend=BACKEND, optimization_level=0)
    return backend_sim.run(circ, shots=SHOTS, optimization_level=0).result().get_counts()

def execute_circ_no_depol(circ, SHOTS=10000):
    # BACKEND=FakeGuadalupe()
    # depol_noise_model=depol_noise(prob1)
    # print(depol_noise_model)
    # backend_sim=AerSimulator.from_backend(BACKEND)
    backend_sim=AerSimulator()
    # grover_qc=transpile(grover_qc, backend=BACKEND, optimization_level=0)
    return backend_sim.run(circ, shots=SHOTS, optimization_level=0).result().get_counts()

def execute_circ_with_depol_1q_reptest(circ, prob1, SHOTS=10000):
    # BACKEND=FakeGuadalupe()
    depol_noise_model=depol_noise_1q_reptest(prob1)
    # print(depol_noise_model)
    # backend_sim=AerSimulator.from_backend(BACKEND)
    backend_sim=AerSimulator(noise_model=depol_noise_model)
    # grover_qc=transpile(grover_qc, backend=BACKEND, optimization_level=0)
    return backend_sim.run(circ, shots=SHOTS, optimization_level=0).result().get_counts()


def depol_noise(one_qerr_prob):
    noise_model=NoiseModel()
    err1=depolarizing_error(one_qerr_prob, 1)
    err2=depolarizing_error(10*one_qerr_prob, 2)
    noise_model.add_all_qubit_quantum_error(err1, ["s", "sdg", "h", "t", "tdg", "sx", "sxdg", "x", "y", "z", "rz"])
    noise_model.add_all_qubit_quantum_error(err2, ["cx"])
    return noise_model

def depol_noise_1q_reptest(one_qerr_prob):
    noise_model=NoiseModel()
    err1=depolarizing_error(one_qerr_prob, 1)
    err2=depolarizing_error(10*one_qerr_prob, 2)
    noise_model.add_all_qubit_quantum_error(err1, ["s", "sdg", "h", "t", "x", "y", "z", "id"])
    noise_model.add_all_qubit_quantum_error(err2, ["cx"])
    return noise_model

# Grover stuff ###########################################
def initialize_s(qc, qubits):
    """Apply a H-gate to 'qubits' in qc"""
    for q in qubits:
        qc.h(q)
    return qc

def grover_reps(reps):
    n = 2
    grover_circuit = QuantumCircuit(n)
    grover_circuit = initialize_s(grover_circuit, [0,1])
    for _ in range(reps):
        # grover_circuit.draw()
        # grover_circuit.draw()
        grover_circuit.cz(0,1) # Oracle
        # Diffusion operator (U_s)
        grover_circuit.h([0,1])
        grover_circuit.z([0,1])
        grover_circuit.cz(0,1)
        grover_circuit.h([0,1])

    grover_circuit.measure_all()


    return grover_circuit


def grover():
    n = 2
    grover_circuit = QuantumCircuit(n)
    grover_circuit = initialize_s(grover_circuit, [0,1])
    # grover_circuit.draw()
    grover_circuit.cz(0,1) # Oracle
    # grover_circuit.draw()
    # Diffusion operator (U_s)
    grover_circuit.h([0,1])
    grover_circuit.z([0,1])
    grover_circuit.cz(0,1)
    grover_circuit.h([0,1])

    grover_circuit.measure_all()

    return grover_circuit

# QFT stuff
def qft_rotations(circuit, n):
    """Performs qft on the first n qubits in circuit (without swaps)"""
    if n == 0:
        return circuit
    n -= 1
    circuit.h(n)
    for qubit in range(n):
        circuit.cp(pi/2**(n-qubit), qubit, n)
    # At the end of our function, we call the same function again on
    # the next qubits (we reduced n by one earlier in the function)
    qft_rotations(circuit, n)

def swap_registers(circuit, n):
    for qubit in range(n//2):
        circuit.swap(qubit, n-qubit-1)
    return circuit

def qft_init(nqubits):
    """QFT on the first n qubits in circuit with initialization."""
    # Create the circuit
    # nqubits = 3
    qc = QuantumCircuit(nqubits)
    for idx in range(nqubits):
        qc.x(idx)
    # # qc.draw()
    # qft_rotations(qc, n)
    # swap_registers(qc, n)
    return qc

def qft(qc, n):
    """QFT on the first n qubits in circuit"""
    # Encode the state 5
    # qc.x(0)
    # qc.x(2)
    # qc.draw()
    qft_rotations(qc, n)
    swap_registers(qc, n)
    return qc

def inverse_qft(circ, n):
    """Does the inverse QFT on the first n qubits in circuit"""
    # First we create a QFT circuit of the correct size:
    qft_circ = qft(QuantumCircuit(n), n)
    # Then we take the inverse of this circuit
    invqft_circ = qft_circ.inverse()
    # And add it to the first n qubits in our existing circuit
    circ.append(invqft_circ, circ.qubits[:n])
    circ.measure_all()
    return circ.decompose() # .decompose() allows us to see the individual gates