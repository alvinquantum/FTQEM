#!/usr/bin/env python3

# import circuit_cutter
# import mlrecon_methods as ml
import numpy as np
from copy import deepcopy
from qiskit import transpile, QuantumCircuit
from qiskit.circuit import QuantumRegister
from qiskit.circuit.classicalregister import ClassicalRegister
from qiskit import QuantumCircuit, QuantumRegister, transpile
from qiskit_aer import AerSimulator
from qiskit.converters import circuit_to_dag
from qiskit.dagcircuit.dagnode import DAGOpNode
import math
import utilities as commonUtl

np.set_printoptions(linewidth = 200)

def post_select(result, logical_bases={"00": "0", "11": "1"}):
    '''
    Post selects the results based on the logical outcomes. 
    '''
    pqubits_per_lqubit=7
    # print(pqubits_per_lqubit)
    # print(logical_bases)
    filtered_result={}
    sum=0
    for key, value in result.items():
        # We only need the data registers.
        key_split=key.split()
        if len(key_split)>1:
            key=key_split[1]
            # print("key: ", key)

        physical_keys=commonUtl.chunkwise(key, pqubits_per_lqubit)
        checkif=all(k1[:-1]=="0"*6 for k1 in physical_keys)
        # print("compute keys: ", physical_keys)
        if checkif:
            decoded_key=[k[-1] for k in physical_keys]
            concated_key="".join(decoded_key)
            # print("decoded key: ", decoded_key)
            # print("concat key: ", concated_key)
            filtered_result[concated_key]=filtered_result.get(concated_key, 0)+value
    for key, val in filtered_result.items():
        sum+=val
    # print("post select runs: ", sum)
    return filtered_result

def post_select_concat(result, concat_reps):
    '''
    Post selects the results based on the logical outcomes. 
    '''
    for _ in range(concat_reps+1):
        result=post_select(result)
    return result

class LOps:
    @staticmethod
    def hadamard_l(circ, qubit):
        for q in qubit:
            circ.h(q)

    @staticmethod
    def toff_l(circ, q1, q2, target):
        for idx, _ in enumerate(q1):
            circ.ccx(q1[idx], q2[idx], target[idx])

    @staticmethod
    def cnot_l(circ, q1, q2):
        for idx, _ in enumerate(q1):
            circ.cx(q1[idx], q2[idx])

    @staticmethod
    def sgate_l(circ, qubit):
        for q in qubit:
            circ.s(q)

    @staticmethod
    def sdggate_l(circ, qubit):
        for q in qubit:
            circ.sdg(q)

    # @staticmethod
    # def tgate_l(circ, qubit):
    #     pass

    @staticmethod
    def xgate_l(circ, qubit):
        for q in qubit:
            circ.x(q)

    @staticmethod
    def ygate_l(circ, qubit):
        for q in qubit:
            circ.y(q)

    @staticmethod
    def zgate_l(circ, qubit):
        for q in qubit:
            circ.z(q)

    @staticmethod
    def rzgate_l(circ, phi, qubit):
        # pass
        if float(2*phi/math.pi) % 1 == 0:
            for q in qubit:
                circ.rz(phi, q)

            # print("equal to zero: ", phi)
        else:
            # print(float(2*phi/math.pi) % 1)
            print("non transversal rz. ")
        

    @staticmethod
    def mgate_l(circ, qubit, reg):
        for idx, q in enumerate(qubit):
            circ.measure(q, reg[idx])

    @staticmethod
    def sxgate_l(circ, qubit):
        '''Logical Sqrt(x) gate.'''
        for q in qubit:
            circ.sx(q)

    @staticmethod
    def sxdggate_l(circ, qubit):
        '''Logical Sqrt(x) gate.'''
        for q in qubit:
            circ.sxdg(q)

def decode_steane(circ, lqubit):
    # for q in lqubit:
    circ.cx(lqubit[3],lqubit[6])
    circ.cx(lqubit[3],lqubit[5])
    circ.cx(lqubit[3],lqubit[4])
    circ.cx(lqubit[2],lqubit[6])
    circ.cx(lqubit[2],lqubit[4])
    circ.cx(lqubit[2],lqubit[0])
    circ.cx(lqubit[1],lqubit[5])
    circ.cx(lqubit[1],lqubit[4])
    circ.cx(lqubit[1],lqubit[0])
    circ.cx(lqubit[0],lqubit[6])
    circ.cx(lqubit[0],lqubit[5])
    circ.h([lqubit[1],lqubit[2],lqubit[3]])


def to_encoded_circ_IBMBasis(circ):
    num_physq_per_logq=7
    num_lqubits= len(circ.qubits)
    num_pqubits=num_physq_per_logq*num_lqubits
    phys_qr=QuantumRegister(num_pqubits)
    phys_cr=ClassicalRegister(num_pqubits)
    phys_qc=QuantumCircuit(phys_qr, phys_cr)
    logical_qubits=commonUtl.chunkwise(phys_qr, num_physq_per_logq)
    logical_cregs=commonUtl.chunkwise(phys_cr, num_physq_per_logq)
    for qubit in logical_qubits:
        to_steane_logical0(phys_qc, qubit)
    phys_qc.barrier()

    # print(logical_qubits)
    qc_dag=circuit_to_dag(circ)
    layers= list(qc_dag.multigraph_layers())
    # Iterate through the layers and convert to the encoding.
    for layer in layers:
        for node in layer:
            if type(node)==DAGOpNode:
                if node.name=="barrier":
                    phys_qc.barrier()
                elif node.name=="cx":
                    # print("cx")
                    ctrl_idx=node.qargs[0].index
                    tar_idx=node.qargs[1].index
                    # print(logical_qubits[ctrl_idx])
                    # print(logical_qubits[tar_idx])
                    LOps.cnot_l(phys_qc, logical_qubits[ctrl_idx], logical_qubits[tar_idx])
                elif node.name =="h":
                    LOps.hadamard_l(phys_qc, logical_qubits[node.qargs[0].index])
                elif node.name=="rz":
                    # print("h")
                    # print(logical_qubits[node.qargs[0].index])
                    # print(node)
                    # print(node.op.params)
                    # print(type(node.op))
                    LOps.rzgate_l(phys_qc, node.op.params[0], logical_qubits[node.qargs[0].index])
                elif node.name=="sx":
                    # print("s")
                    LOps.sxgate_l(phys_qc, logical_qubits[node.qargs[0].index])
                elif node.name=="sxdg":
                    # print("sdg")
                    LOps.sxdggate_l(phys_qc, logical_qubits[node.qargs[0].index])
                elif node.name=="x":
                    # print("t")
                    LOps.xgate_l(phys_qc, logical_qubits[node.qargs[0].index])
                elif node.name=="measure":
                    # print("m")
                    # print(node.qargs)
                    # print(node.cargs)
                    continue
                    # LOps.mgate_l(phys_qc, logical_qubits[node.qargs[0].index], logical_cregs[node.cargs[0].index])
                else:
                    assert False, f"{node.name} gate not recognized."

    phys_qc.barrier()
    

    for qubit in logical_qubits:
        decode_steane(phys_qc, qubit)
    phys_qc.barrier()
    phys_qc.measure_all(add_bits=False)
    basis2 = ["cx", "id", "rz", "sx", "x"]

    return transpile(phys_qc, basis_gates=basis2, optimization_level=0)

def to_encoded_circ_CliffT(circ):
    num_physq_per_logq=7
    num_lqubits= len(circ.qubits)
    num_pqubits=num_physq_per_logq*num_lqubits
    phys_qr=QuantumRegister(num_pqubits)
    phys_cr=ClassicalRegister(num_pqubits)
    phys_qc=QuantumCircuit(phys_qr, phys_cr)
    logical_qubits=commonUtl.chunkwise(phys_qr, num_physq_per_logq)
    logical_cregs=commonUtl.chunkwise(phys_cr, num_physq_per_logq)
    for qubit in logical_qubits:
        to_steane_logical0(phys_qc, qubit)

    # print(logical_qubits)
    qc_dag=circuit_to_dag(circ)
    layers= list(qc_dag.multigraph_layers())
    # Iterate through the layers and convert to the encoding.
    for layer in layers:
        for node in layer:
            if type(node)==DAGOpNode:
                if node.name=="barrier":
                    phys_qc.barrier()
                elif node.name=="cx":
                    # print("cx")
                    ctrl_idx=node.qargs[0].index
                    tar_idx=node.qargs[1].index
                    # print(logical_qubits[ctrl_idx])
                    # print(logical_qubits[tar_idx])
                    LOps.cnot_l(phys_qc, logical_qubits[ctrl_idx], logical_qubits[tar_idx])
                elif node.name=="h":
                    # print("h")
                    # print(logical_qubits[node.qargs[0].index])
                    LOps.hadamard_l(phys_qc, logical_qubits[node.qargs[0].index])
                elif node.name=="s":
                    # print("s")
                    LOps.sgate_l(phys_qc, logical_qubits[node.qargs[0].index])
                elif node.name=="sdg":
                    # print("sdg")
                    LOps.sdggate_l(phys_qc, logical_qubits[node.qargs[0].index])
                elif node.name=="t":
                    # print("t")
                    LOps.tgate_l(phys_qc, logical_qubits[node.qargs[0].index])
                elif node.name=="measure":
                    continue
                    # print("m")
                    # print(node.qargs)
                    # print(node.cargs)
                    # LOps.mgate_l(phys_qc, logical_qubits[node.qargs[0].index], logical_cregs[node.cargs[0].index])
                else:
                    assert False, f"{node.name} gate not recognized."

    phys_qc.barrier()
    

    for qubit in logical_qubits:
        decode_steane(phys_qc, qubit)
    phys_qc.barrier()
    phys_qc.measure_all(add_bits=False)
    basis1 = ["h", "s", "sdg", "cx", "t", "cxx"]

    return transpile(phys_qc, basis_gates=basis1, optimization_level=0)

def concat_IBMBasis(circ, reps):
    for idx in range(reps+1): #+1 because reps 1 means level of concat
        circ=to_encoded_circ_IBMBasis(circ)
        # print(circ)
    return circ


def insert_error(circ, op, qubit_idx):
    # Trivial gates
    circ.barrier()
    assert op in ["x", "y", "z"], "Error op not x, y, or z."
    if op == "x":
        circ.x(qubit_idx)
    elif op == "y":
        circ.y(qubit_idx)
    else:
        circ.z(qubit_idx)

def create_cat(circ, q1, q2):
    # Cat state
    circ.h(q1)
    for q2 in range(q1+1, q2):
        circ.cx(q1,q2)

def create_cat_ancillas(circ):
    create_cat(circ, 7 , 11)
    create_cat(circ, 11 , 15)
    create_cat(circ, 15 , 19)
    create_cat(circ, 19 , 23)
    create_cat(circ, 23 , 27)
    create_cat(circ, 27 , 31)

def measure_stabil(circ):
    #stabilizer 1
    circ.barrier()
    circ.cx(7, 3)
    circ.cx(8, 4)
    circ.cx(9, 5)
    circ.cx(10, 6)
    #stabilizer 2
    circ.barrier()
    circ.cx(11, 1)
    circ.cx(12, 2)
    circ.cx(13, 5)
    circ.cx(14, 6)
    # stabilizer 3
    circ.barrier()
    circ.cx(15, 0)
    circ.cx(16, 2)
    circ.cx(17, 4)
    circ.cx(18, 6)
    # stabilizer 4
    circ.barrier()
    circ.cz(19, 3)
    circ.cz(20, 4)
    circ.cz(21, 5)
    circ.cz(22, 6)
    # stabilizer 5
    circ.barrier()
    circ.cz(23, 1)
    circ.cz(24, 2)
    circ.cz(25, 5)
    circ.cz(26, 6)
    # stabilizer 6
    circ.barrier()
    circ.cz(27, 0)
    circ.cz(28, 2)
    circ.cz(29, 4)
    circ.cz(30, 6)
    circ.barrier()
    # Undo controlled nots
    q1s=[7,11,15,19,23,27]
    for q1 in q1s:
        targets=reversed(list(range(q1+1, q1+4)))
        circ.cx(q1, targets)
        circ.h(q1)

def to_steane_logical0(circ, logical_qubit):
    q=logical_qubit
    circ.h([q[1],q[2],q[3]])
    circ.cx(q[1],q[0])
    circ.cx(q[3],q[5])
    circ.cx(q[2],q[6])
    circ.cx(q[1],q[4])
    circ.cx(q[2],q[0])
    circ.cx(q[3],q[6])
    circ.cx(q[1],q[5])
    circ.cx(q[6],q[4])

def generate_syndromes(circ, aregs):
    temp_circ=deepcopy(circ)
    temp_circ.barrier()
    error_ops=["x", "y", "z"]
    syndromes={}
    for idx in range(7):
        for op in error_ops:
            circ_copy=deepcopy(temp_circ)    
            insert_error(circ_copy, op, idx)
            measure_stabil(circ_copy)
            measure_ancillas(circ_copy, aregs)
            sim_ideal=AerSimulator()
            ideal_result = sim_ideal.run(circ_copy).result()
            ideal_dist=ideal_result.get_counts()
            syndromes.update({list(ideal_dist.keys())[0]: op+f"_{idx}"})
    return syndromes

def correct_steane_cond_gates(circ, syndromes_1qubit, cregs):
    circ.barrier()
    clbits=circ.clbits
    print(clbits)
    print(cregs)
    for syndrome, error_code in syndromes_1qubit.items():
        # if error_code:
        print(error_code)
        op, idx=error_code.split("_")
        idx=int(idx)
        # syndrome=int(syndrome, 2)
        print(syndrome)
        controls=control_idxs_from_syndrome(syndrome)
        if op=="x":
            circ.cx(controls, idx)#.c_if(cregs, syndrome)
        elif op=="y":
            circ.cy(controls, idx)#.c_if(cregs, syndrome)
        elif op=="z":
            circ.cz(controls, idx)#.c_if(cregs, syndrome)
        # else:
        #     print("not a single qubit error.")
    # Measure
    circ.barrier()

def measure_ancillas(circ, aregs):
    circ.barrier()
    creg=0
    for idx in range(0,24):
        circ.measure(aregs[idx],creg)
        creg+=1
    circ.barrier()

def control_idxs_from_syndrome(syndrome):
    return [idx+7 for idx, elem in enumerate(reversed(syndrome)) if elem=="1"]

def correct_steane(circ, sydromes_1qubit, syndrome):
    if syndrome =="0"*24:
        print("no error. ")
        return
    error_code=sydromes_1qubit.get(syndrome)
    if error_code:
        print(error_code)
        op, idx=error_code.split("_")
        idx=int(idx)
        if op=="x":
            circ.x(idx)
        elif op=="y":
            circ.y(idx)
        else:
            circ.z(idx)
    else:
        print("not a single qubit error.")