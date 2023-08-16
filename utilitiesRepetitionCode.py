#!/usr/bin/env python3
#initialization

from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister
from qiskit.converters import circuit_to_dag
from qiskit.dagcircuit.dagnode import DAGOpNode
from qiskit_aer.noise import depolarizing_error, QuantumError
from copy import deepcopy
import math
import utilities as commonUtl

def post_select_dsm(result, pqubits_per_lqubit):
    '''
    Post selects the results based on the logical outcomes. 
    '''
    # print(pqubits_per_lqubit)
    # print(logical_bases)
    filtered_result={}
    sum=0
    for key, value in result.items():
        # # We only need the data registers.
        # key_split=key.split()
        # if len(key_split)>1:
        #     key=key_split[0]
        #     # print("key: ", key)

        physical_keys=commonUtl.chunkwise(key, pqubits_per_lqubit)
        # Assuming that the data qubit after decoding is the topmost in a logical qubit, which is
        # the rightmost classical register. Thus, the ancillas in a logical qubit are everything before
        # the last physical qubit.
        checkif=all(k1[:-1] in ["0"*(pqubits_per_lqubit-1)] for k1 in physical_keys)
        # print("compute keys: ", physical_keys)
        if checkif:
            # the decoded value is then the last physical qubit.
            decoded_key=[k[-1] for k in physical_keys]
            concated_key="".join(decoded_key)
            # print("decoded key: ", decoded_key)
            # print("concat key: ", concated_key)
            filtered_result[concated_key]=filtered_result.get(concated_key, 0)+value
    for key, val in filtered_result.items():
        sum+=val
    # print("post select runs: ", sum)
    return filtered_result

def post_select(result, pqubits_per_lqubit, logical_bases={"00": "0", "11": "1"}):
    '''
    Post selects the results based on the logical outcomes. 
    '''
    # print(pqubits_per_lqubit)
    # print(logical_bases)
    filtered_result={}
    sum=0
    for key, value in result.items():
        # We only need the data registers.
        key_split=key.split()
        if len(key_split)>1:
            key=key_split[0]
            # print("key: ", key)

        physical_keys=commonUtl.chunkwise(key, pqubits_per_lqubit)
        checkif=all(k1 in list(logical_bases.keys()) for k1 in physical_keys)
        # print("compute keys: ", physical_keys)
        if checkif:
            decoded_key=[logical_bases[k] for k in physical_keys]
            concated_key="".join(decoded_key)
            # print("decoded key: ", decoded_key)
            # print("concat key: ", concated_key)
            filtered_result[concated_key]=filtered_result.get(concated_key, 0)+value
    for key, val in filtered_result.items():
        sum+=val
    # print("post select runs: ", sum)
    return filtered_result

def post_select_stab_meas(result, pqubits_per_lqubit, logical_bases={"00": "0", "11": "1"}):
    '''
    Post selects the results based on the logical outcomes. For standard stabilizer measurement in the 
    repetition code.
    '''
    # print(pqubits_per_lqubit)
    # print(logical_bases)
    filtered_result={}
    sum=0
    for key, value in result.items():
        # We only need the data registers.
        # print("key before: ", key)
        key_split=key.split()
        if len(key_split)>1:
            key=key_split[0]
            print("key: ", key_split)

        # The stabilizers are ZZ so the ancilla cat states are all 2.       
        physical_ancilla_keys=commonUtl.chunkwise(key, 2)
        print(physical_ancilla_keys)
        # check if the ancillas for each logical qubit are all zeros since we want to lie in the positive eigenspace.
        checkif=all(k1 in ["0"*2] for k1 in physical_ancilla_keys)
        # print("compute keys: ", physical_keys)
        if checkif:
            # decoded_key=[logical_bases[k] for k in physical_keys]
            # concated_key="".join(decoded_key)
            concated_key=key_split[1]
            # print("decoded key: ", decoded_key)
            print("concat key: ", concated_key)
            filtered_result[concated_key]=filtered_result.get(concated_key, 0)+value
    for key, val in filtered_result.items():
        sum+=val
    # print("post select runs: ", sum)
    return filtered_result

def correct(result, pqubits_per_lqubit):
    '''
    Post selects the results based on the logical outcomes. 
    '''
    # It's a normal distribution so calculate 1std deviation from the mean. Use it
    # as the bounds for the error correction.
    mean_of_x=pqubits_per_lqubit/2
    number_of_elems=2**pqubits_per_lqubit
    mean_of_xsquared=sum([x**2*math.comb(pqubits_per_lqubit, x) for x in range(0, pqubits_per_lqubit+1, 1)])/number_of_elems
    # print("num phys qubits: ", pqubits_per_lqubit)
    # print( "vals: ", [x**2*math.comb(pqubits_per_lqubit, x) for x in range(0, pqubits_per_lqubit+1, 1)])
    # print("range", list[range(0, pqubits_per_lqubit+1, 1)])
    # for val_temp in range(0, pqubits_per_lqubit+1, 1):
    #     print("xval: ", val_temp)
    #     print("combo: ", math.comb(pqubits_per_lqubit, val_temp))
    # print("mean of x^2: ", mean_of_xsquared)
    # print("mean of x: ", mean_of_x)
    # print("square of mean of x: ", mean_of_x**2)
    std_dev=math.sqrt(mean_of_xsquared-mean_of_x**2)

    print("std_dev: ", std_dev)
    print("lower condition: ", pqubits_per_lqubit/2-std_dev)
    print("upper condition: ", pqubits_per_lqubit/2+std_dev)
    print()
    # filtered_result={}
    # sum=0
    corrected_results=deepcopy(result)
    # print(type(items))
    for key, value in result.items():
        # We only need the data registers.
        key_split=key.split()
        if len(key_split)>1:
            key=key_split[1]
            # print("key: ", key)

        # print("split: ", key_split)
        physical_keys=commonUtl.chunkwise(key, pqubits_per_lqubit)
        # print(physical_keys)
        for idx, val in enumerate(physical_keys):
            weight=val.count("1")

            # Using 1 std deviation
            if weight<=pqubits_per_lqubit/2-std_dev: #and weight< math.floor((pqubits_per_lqubit-1)/2):
                physical_keys[idx]="0"*pqubits_per_lqubit
            elif weight>=pqubits_per_lqubit/2+std_dev: #and weight> pqubits_per_lqubit-math.floor((pqubits_per_lqubit-1)/2):
                physical_keys[idx]="1"*pqubits_per_lqubit
            else:
                physical_keys[idx]=val

        # print("physical key", physical_keys)
        # print(type(physical_keys[0]))
        concated_key="".join(physical_keys)
        popped_val=corrected_results.pop(key)
        # if popped_val != value:
        #     print("key ", key)
        #     print("popped ", popped_val)
        #     print("value ", value)
        corrected_results[concated_key] = corrected_results.get(concated_key, 0)+popped_val
    return corrected_results

def post_select_ft_pcs(result, pqubits_per_lqubit, logical_bases={"00": "0", "11": "1"}):
    '''
    Post selects the results based on the logical outcomes. 
    '''
    # print(result)
    filtered_result={}
    sum=0
    for key, value in result.items():
        # We only need the data registers.
        key_split=key.split()
        # print(int(key_split[0]))
        # print(type(key_split[0]))
        
        if len(key_split)>1:
            if int(key_split[0])==1: #PCS check
                continue
            key=key_split[1]
            # print("key: ", key)

        physical_keys=commonUtl.chunkwise(key, pqubits_per_lqubit)
        # print("physical key: ", physical_keys)
        checkif=all(k1 in list(logical_bases.keys()) for k1 in physical_keys)
        # print("compute keys: ", physical_keys)
        if checkif:
            decoded_key=[logical_bases[k] for k in physical_keys]
            concated_key="".join(decoded_key)
            # print("decoded key: ", decoded_key)
            # print("concat key: ", concated_key)
            filtered_result[concated_key]=filtered_result.get(concated_key, 0)+value
    for key, val in filtered_result.items():
        sum+=val
    # print("post select runs: ", sum)
    return filtered_result


class LOps:
    @staticmethod
    def hadamard_l(circ, qubit):
        if len(qubit)==1:
            q1=qubit[0]
        else:
            q1=qubit[0]
            remain_qubits=qubit[1:]
            for q2 in remain_qubits:
                circ.cx(q1, q2)
            circ.h(q1)
            for q2 in remain_qubits:
                circ.cx(q1, q2)

    @staticmethod
    def hadamard_fault_tol_pcs_l(circ,  qubit_t, qubit_pcs, hadamard_pcs_cr):
        # circ.reset(qubit_pcs)
        # circ.h(qubit_pcs)
        circ.cx(qubit_pcs, qubit_t[-1])
        circ.barrier()
        q1=qubit_t[0]
        remain_qubits=qubit_t[1:]
        for q2 in remain_qubits:
            circ.cx(q1, q2)
        circ.h(q1)
        for q2 in remain_qubits:
            circ.cx(q1, q2)
        circ.barrier()
        circ.cx(qubit_pcs, qubit_t[-1])
        # circ.h(qubit_pcs)
        # circ.measure(qubit_pcs, hadamard_pcs_cr)


    @staticmethod
    def hadamard_fault_tol_l(circ, qubit_p, qubit_m, qubit_t, hadamard_reg):
        circ.barrier()
        circ.reset(qubit_m)
        circ.reset(qubit_p)
        Encodings.create_lminus_state(circ, qubit_m)
        # circ.barrier()
        Encodings.create_lplus_state(circ, qubit_p)
        circ.barrier()
        LOps.toff_l(circ, qubit_p, qubit_t, qubit_m)
        # circ.barrier()
        LOps.hadamard_l(circ, qubit_t)
        LOps.cnot_l(circ, qubit_t, qubit_p)
        LOps.hadamard_l(circ, qubit_t)
        LOps.mgate_l(circ, qubit_t, hadamard_reg)
        # with circ.if_test((hadamard_reg, 0b11)):
        #     xgate_l(circ, qubit_t)

    @staticmethod
    def toff_l(circ, q1, q2, target):
        for idx, _ in enumerate(q1):
            circ.ccx(q1[idx], q2[idx], target[idx])

    @staticmethod
    def cnot_l(circ, q1, q2):
        for idx, _ in enumerate(q1):
            circ.cx(q1[idx], q2[idx])

    @staticmethod
    def cnot_l_not_fault_tol(circ, q1, q2):
        for idx, _ in enumerate(q1):
            circ.cx(q1[0], q2[idx])

    @staticmethod
    def sgate_l(circ, qubit):
        circ.s(qubit[0])

    @staticmethod
    def sgate_l_reptest(circ, qubit):
        circ.s(qubit[0])
        for idx in qubit[1:]:
            circ.id(idx)

    @staticmethod
    def sdggate_l(circ, qubit):
        circ.sdg(qubit[0])
    
    @staticmethod
    def sdggate_l_reptest(circ, qubit):
        circ.sdg(qubit[0])
        for idx in qubit[1:]:
            circ.id(idx)

    @staticmethod
    def tgate_l(circ, qubit):
        circ.t(qubit[0])

    @staticmethod
    def tdggate_l(circ, qubit):
        circ.tdg(qubit[0])

    @staticmethod
    def tgate_l_reptest(circ, qubit):
        circ.t(qubit[0])
        for idx in qubit[1:]:
            circ.id(idx)

    @staticmethod
    def xgate_l(circ, qubit):
        for idx, _ in enumerate(qubit):
            circ.x(qubit[idx])

    @staticmethod
    def ygate_l(circ, qubit):
        if len(qubit)==1:
            circ.y(qubit[0])
        else:
            circ.y(qubit[0])
            for idx, _ in enumerate(qubit[1:]):
                circ.x(qubit[idx])

    @staticmethod
    def zgate_l(circ, qubit):
        circ.z(qubit[0])

    @staticmethod
    def rzgate_l(circ, phi, qubit):
        circ.rz(phi, qubit[0])

    @staticmethod
    def mgate_l(circ, qubit, reg):
        for idx, q in enumerate(qubit):
            circ.measure(q, reg[idx])

    @staticmethod
    def sxgate_l(circ, qubit):
        '''Logical Sqrt(x) gate.'''
        q1=qubit[0]
        if len(qubit)==1:
            circ.sx(q1)
        else:
            remain_qubits=qubit[1:]
            for q2 in remain_qubits:
                circ.cx(q1, q2)
            circ.sx(q1)
            for q2 in remain_qubits:
                circ.cx(q1, q2)

    @staticmethod
    def sxdggate_l(circ, qubit):
        '''Logical Sqrt(x) gate.'''
        q1=qubit[0]
        if len(qubit)==1:
            circ.sxdg(q1)
        else:
            remain_qubits=qubit[1:]
            for q2 in remain_qubits:
                circ.cx(q1, q2)
            circ.sxdg(q1)
            for q2 in remain_qubits:
                circ.cx(q1, q2)

class LOpsWithDepol:
    @staticmethod
    def hadamard_l(circ, qubit, quant_errs):
        q1=qubit[0]
        if len(qubit)==1:
            circ.h(q1)
            circ.append(quant_errs[0], [q1])
        # print(q1)
        # circ.append(quant_errs[0], qubit[0])
        else:
            remain_qubits=qubit[1:]
            for q2 in remain_qubits:
                circ.cx(q1, q2)
                circ.append(quant_errs[1], [q1, q2])
            circ.h(q1)
            circ.append(quant_errs[0], [q1])
            for q2 in remain_qubits:
                circ.cx(q1, q2)
                circ.append(quant_errs[1], [q1, q2])


    @staticmethod
    def hadamard_fault_tol_l(circ, qubit_p, qubit_m, qubit_t, hadamard_reg, quant_errs):
        circ.barrier()
        circ.reset(qubit_m)
        circ.reset(qubit_p)
        Encodings.create_lminus_state(circ, qubit_m)
        # circ.barrier()
        Encodings.create_lplus_state(circ, qubit_p)
        circ.barrier()
        LOpsWithDepol.toff_l(circ, qubit_p, qubit_t, qubit_m, quant_errs)
        circ.barrier()
        LOpsWithDepol.hadamard_l(circ, qubit_t, quant_errs)
        LOps.mgate_l(circ, qubit_t, hadamard_reg)
        circ.barrier()
        LOpsWithDepol.cnot_l(circ, qubit_t, qubit_p, quant_errs)
        # LOpsWithDepol.hadamard_l(circ, qubit_t, quant_errs)
        # with circ.if_test((hadamard_reg, 0b11)):
        #     xgate_l(circ, qubit_t)

    @staticmethod
    def toff_l(circ, q1, q2, target, quant_errs):
        for idx, _ in enumerate(q1):
            circ.ccx(q1[idx], q2[idx], target[idx])
            circ.append(quant_errs[2], [q1[idx], q2[idx], target[idx]])

    @staticmethod
    def cnot_l(circ, q1, q2, quant_errs):
        for idx, _ in enumerate(q1):
            circ.cx(q1[idx], q2[idx])
            circ.append(quant_errs[1], [q1[idx], q2[idx]])

    @staticmethod
    def sgate_l(circ, qubit, quant_errs):
        circ.s(qubit[0])
        circ.append(quant_errs[0], [qubit[0]])

    @staticmethod
    def sgate_l_reptest(circ, qubit):
        circ.s(qubit[0])
        for idx in qubit[1:]:
            circ.id(idx)

    @staticmethod
    def sdggate_l(circ, qubit, quant_errs):
        circ.sdg(qubit[0])
        circ.append(quant_errs[0], [qubit[0]])
    
    @staticmethod
    def sdggate_l_reptest(circ, qubit):
        circ.sdg(qubit[0])
        for idx in qubit[1:]:
            circ.id(idx)

    @staticmethod
    def tgate_l(circ, qubit, quant_errs):
        circ.t(qubit[0])
        circ.append(quant_errs[0], [qubit[0]])

    @staticmethod
    def tdggate_l(circ, qubit, quant_errs):
        circ.tdg(qubit[0])
        circ.append(quant_errs[0], [qubit[0]])

    @staticmethod
    def tgate_l_reptest(circ, qubit):
        circ.t(qubit[0])
        for idx in qubit[1:]:
            circ.id(idx)

    @staticmethod
    def xgate_l(circ, qubit, quant_errs):
        for idx, _ in enumerate(qubit):
            circ.x(qubit[idx])
            circ.append(quant_errs[0], [qubit[idx]])


    @staticmethod
    def ygate_l(circ, qubit, quant_errs):
        if len(qubit)==1:
            circ.y(qubit[0])
        else:
            circ.y(qubit[0])
            circ.append(quant_errs[0], [qubit[0]])
            for idx, _ in enumerate(qubit[1:]):
                circ.x(qubit[idx])
                circ.append(quant_errs[0], [qubit[idx]])


    @staticmethod
    def zgate_l(circ, qubit, quant_errs):
        circ.z(qubit[0])
        circ.append(quant_errs[0], [qubit[0]])

    @staticmethod
    def rzgate_l(circ, phi, qubit, quant_errs):
        circ.rz(phi, qubit[0])
        circ.append(quant_errs[0], [qubit[0]])

    @staticmethod
    def mgate_l(circ, qubit, reg):
        for idx, q in enumerate(qubit):
            circ.measure(q, reg[idx])

    @staticmethod
    def sxgate_l(circ, qubit, quant_errs):
        '''Logical Sqrt(x) gate.'''
        q1=qubit[0]
        if len(qubit)==1:
            circ.sxdg(q1)
            circ.append(quant_errs[0], [q1])
        else:
            remain_qubits=qubit[1:]
            for q2 in remain_qubits:
                circ.cx(q1, q2)
                circ.append(quant_errs[1], [q1, q2])
            circ.sx(q1)
            circ.append(quant_errs[0], [q1])
            for q2 in remain_qubits:
                circ.cx(q1, q2)
                circ.append(quant_errs[1], [q1, q2])

    @staticmethod
    def sxdggate_l(circ, qubit, quant_errs):
        '''Logical Sqrt(x) gate.'''
        q1=qubit[0]
        if len(qubit)==1:
            circ.sxdg(q1)
            circ.append(quant_errs[0], [q1])
        else:
            remain_qubits=qubit[1:]
            for q2 in remain_qubits:
                circ.cx(q1, q2)
                circ.append(quant_errs[1], [q1, q2])
            circ.sxdg(q1)
            circ.append(quant_errs[0], [q1])
            for q2 in remain_qubits:
                circ.cx(q1, q2)
                circ.append(quant_errs[1], [q1, q2])


class Encodings:
    @staticmethod
    def to_encoded_circ_CliffT_depol_stab_meas(circ, reps, p1_err, p2_err, p3_err):
        err1=depolarizing_error(p1_err, 1)
        err2=depolarizing_error(p2_err, 2)
        err3=depolarizing_error(p3_err, 3)
        quant_err1=QuantumError(err1).to_instruction()
        quant_err2=QuantumError(err2).to_instruction()
        quant_err3=QuantumError(err3).to_instruction()
        quant_errs=[quant_err1, quant_err2, quant_err3]
        num_lqubits= len(circ.qubits)
        num_pqubits=reps*num_lqubits
        phys_qr=QuantumRegister(num_pqubits)
        #ancillas for stabilizer measurements. the number of generators is reps-1. Each generator will be a ZZ. So we need 2*(num_lqubits-1)*num_lqubits
        # in total 
        num_gen=(reps-1)
        phys_ancil_qr1=QuantumRegister(2*num_gen*num_lqubits)
        #ancillas for observable measurement of logical sigma_z.
        phys_ancil_qr2=QuantumRegister(num_lqubits)
        phys_cr=ClassicalRegister(num_pqubits)
        #cr for stabilizer measurements. 
        phys_ancil_cr1=ClassicalRegister(2*num_gen*num_lqubits)
        #cr for the observable measurements.
        phys_ancil_cr2=ClassicalRegister(num_lqubits)

        phys_qc=QuantumCircuit(phys_qr, phys_ancil_qr1, phys_ancil_qr2, phys_cr, phys_ancil_cr2, phys_ancil_cr1)
        logical_qubits=commonUtl.chunkwise(phys_qr, reps)
        #Each generator will be a ZZ so each lqubit requires 2*num_gen physical ancillas.
        logical_ancil_qubits1=commonUtl.chunkwise(phys_ancil_qr1, 2*num_gen)
        logical_cregs=commonUtl.chunkwise(phys_cr, reps)
        logical_ancil_cregs1=commonUtl.chunkwise(phys_ancil_cr1, 2*num_gen)

        Encodings.init_ancillas(phys_qc, logical_ancil_qubits1, quant_errs)

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
                        LOpsWithDepol.cnot_l(phys_qc, logical_qubits[ctrl_idx], logical_qubits[tar_idx], quant_errs)
                    elif node.name=="h":
                        # print("h")
                        # print(logical_qubits[node.qargs[0].index])
                        # print(node)
                        # print(node.op.params)
                        # print(type(node.op))
                        LOpsWithDepol.hadamard_l(phys_qc,logical_qubits[node.qargs[0].index], quant_errs)
                    elif node.name=="s":
                        # print("s")
                        LOpsWithDepol.sgate_l(phys_qc, logical_qubits[node.qargs[0].index], quant_errs)
                    elif node.name=="sdg":
                        # print("sdg")
                        LOpsWithDepol.sdggate_l(phys_qc, logical_qubits[node.qargs[0].index], quant_errs)
                    elif node.name=="t":
                        # print("t")
                        LOpsWithDepol.tgate_l(phys_qc, logical_qubits[node.qargs[0].index], quant_errs)
                    elif node.name=="tdg":
                        # print("t")
                        LOpsWithDepol.tdggate_l(phys_qc, logical_qubits[node.qargs[0].index], quant_errs)
                    elif node.name=="x":
                        # print("t")
                        LOpsWithDepol.xgate_l(phys_qc, logical_qubits[node.qargs[0].index], quant_errs)
                    elif node.name=="y":
                        # print("t")
                        LOpsWithDepol.ygate_l(phys_qc, logical_qubits[node.qargs[0].index], quant_errs)
                    elif node.name=="z":
                        # print("t")
                        LOpsWithDepol.zgate_l(phys_qc, logical_qubits[node.qargs[0].index], quant_errs)
                    elif node.name=="measure":
                        # we measure at the end after the stabilizer and the observable measurements.
                        continue
                        # LOpsWithDepol.mgate_l(phys_qc, logical_qubits[node.qargs[0].index], logical_cregs[node.cargs[0].index])
                    else:
                        assert False, f"{node.name} gate not recognized."

        phys_qc.barrier()
        Encodings.meas_stabs(phys_qc, logical_qubits, logical_ancil_qubits1, logical_ancil_cregs1, quant_errs)
        phys_qc.barrier()
        Encodings.meas_sigmaz(phys_qc, logical_qubits, phys_ancil_qr2, phys_ancil_cr2, quant_errs)
        phys_qc.barrier()
        for idx, log_qubit in enumerate(logical_qubits):
            phys_qc.measure(log_qubit, logical_cregs[idx])

        return phys_qc

    @staticmethod
    def to_encoded_circ_CliffT_dsm(circ, reps):
        num_lqubits= len(circ.qubits)
        num_pqubits=reps*num_lqubits
        phys_qr=QuantumRegister(num_pqubits)
        phys_cr=ClassicalRegister(num_pqubits)
        phys_qc=QuantumCircuit(phys_qr, phys_cr)
        logical_qubits=commonUtl.chunkwise(phys_qr, reps)
        logical_cregs=commonUtl.chunkwise(phys_cr, reps)

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
                        # print(node)
                        # print(node.op.params)
                        # print(type(node.op))
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
                    elif node.name=="tdg":
                        # print("t")
                        LOps.tdggate_l(phys_qc, logical_qubits[node.qargs[0].index])
                    elif node.name=="x":
                        # print("t")
                        LOps.xgate_l(phys_qc, logical_qubits[node.qargs[0].index])
                    elif node.name=="y":
                        # print("t")
                        LOps.ygate_l(phys_qc, logical_qubits[node.qargs[0].index])
                    elif node.name=="z":
                        # print("t")
                        LOps.zgate_l(phys_qc, logical_qubits[node.qargs[0].index])
                    elif node.name=="measure":
                        # print("m")
                        # print(node.qargs)
                        # print(node.cargs)
                        continue
                        # LOps.mgate_l(phys_qc, logical_qubits[node.qargs[0].index], logical_cregs[node.cargs[0].index])
                    else:
                        assert False, f"{node.name} gate not recognized."
        # Decode the logical states and measure. Note that to test the dsm scheme we post select when the ancillas are all zero.
        # We measure the entire logical qubit because at the very end we calculate the sso.
        for idx, qubit in enumerate(logical_qubits):
            phys_qc.cx(qubit[0], qubit[1:])
            phys_qc.barrier()
            phys_qc.measure(qubit, logical_cregs[idx])
        return phys_qc

    @staticmethod
    def to_encoded_circ_IBMBasis_dsm(circ, reps):
        num_lqubits= len(circ.qubits)
        num_pqubits=reps*num_lqubits
        phys_qr=QuantumRegister(num_pqubits)
        phys_cr=ClassicalRegister(num_pqubits)
        phys_qc=QuantumCircuit(phys_qr, phys_cr)
        logical_qubits=commonUtl.chunkwise(phys_qr, reps)
        logical_cregs=commonUtl.chunkwise(phys_cr, reps)

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
        # Decode the logical states and measure. Note that to test the dsm scheme we post select when the ancillas are all zero.
        # We measure the entire logical qubit because at the very end we calculate the sso.
        for idx, qubit in enumerate(logical_qubits):
            phys_qc.cx(qubit[0], qubit[1:])
            phys_qc.barrier()
            phys_qc.measure(qubit, logical_cregs[idx])
        return phys_qc

    @staticmethod
    def init_ancillas(phys_qc, logical_ancil_qubits, quant_errs):
        # Iterate through list of logical ancillas
        for lqubit in logical_ancil_qubits:
            # For each logical ancilla we divide it into pairs of two because they're all ZZ stabilizers.
            partitioned_lqubit=commonUtl.chunkwise(lqubit, 2)
            # Initialize.
            for pair_qubits in partitioned_lqubit:
                phys_qc.h(pair_qubits[0])
                phys_qc.append(quant_errs[0], [pair_qubits[0]])
                phys_qc.cx(pair_qubits[0], pair_qubits[1])
                phys_qc.append(quant_errs[1], [pair_qubits[0], pair_qubits[1]])


    @staticmethod
    def undo_ancillas(phys_qc, logical_ancil_qubits, logical_ancil_cregs, quant_errs):
        # Iterate through list of logical ancillas
        for lidx, lqubit in enumerate(logical_ancil_qubits):
            # For each logical ancilla we divide it into pairs of two because they're all ZZ stabilizers.
            partitioned_lqubit=commonUtl.chunkwise(lqubit, 2)
            # Decode cat states.
            for pair_qubits in partitioned_lqubit:
                phys_qc.cx(pair_qubits[0], pair_qubits[1])
                phys_qc.append(quant_errs[1], [pair_qubits[0], pair_qubits[1]])
                phys_qc.h(pair_qubits[0])
                phys_qc.append(quant_errs[0], [pair_qubits[0]])
            phys_qc.barrier()
            phys_qc.measure(lqubit, logical_ancil_cregs[lidx])

    @staticmethod
    def meas_stabs(phys_qc, logical_qubits, logical_ancil_qubits, logical_ancil_cregs, quant_errs):
        for lidx, lancil_qubit in enumerate(logical_ancil_qubits):
            ldata_qubit=logical_qubits[lidx]
            partitioned_ancil_lqubit=commonUtl.chunkwise(lancil_qubit, 2)
            # the number of pairs is the number of generators, which is equal to reps-1.
            for pair_idx, pair_ancil_qubits in enumerate(partitioned_ancil_lqubit):
                phys_qc.cz(pair_ancil_qubits[0], ldata_qubit[pair_idx])
                phys_qc.append(quant_errs[1], [pair_ancil_qubits[0], ldata_qubit[pair_idx]])
                phys_qc.cz(pair_ancil_qubits[1], ldata_qubit[pair_idx+1])
                phys_qc.append(quant_errs[1], [pair_ancil_qubits[1], ldata_qubit[pair_idx+1]])
                phys_qc.barrier()
        phys_qc.barrier()
        Encodings.undo_ancillas(phys_qc, logical_ancil_qubits, logical_ancil_cregs, quant_errs)

    @staticmethod
    def meas_sigmaz(phys_qc, logical_qubits, phys_ancil_qr, phys_ancil_cr, quant_errs):
        '''Helper method for to_encoded_circ_IBMBasis_depol_stab_meas.'''
        # Hadamards and noise
        phys_qc.h(phys_ancil_qr)
        for q in phys_ancil_qr:
            phys_qc.append(quant_errs[0], [q])
        # Cz and noise
        for idx, lqubit in enumerate(logical_qubits):
            phys_qc.cz(phys_ancil_qr[idx], lqubit[-1])
            phys_qc.append(quant_errs[1], [phys_ancil_qr[idx], lqubit[-1]])
        # Hadamards and noise
        phys_qc.h(phys_ancil_qr)
        for q in phys_ancil_qr:
            phys_qc.append(quant_errs[0], [q])
        # Measurement
        phys_qc.measure(phys_ancil_qr, phys_ancil_cr)

    @staticmethod
    def to_encoded_circ_IBMBasis_depol_stab_meas(circ, reps, p1_err, p2_err, p3_err):
        err1=depolarizing_error(p1_err, 1)
        err2=depolarizing_error(p2_err, 2)
        err3=depolarizing_error(p3_err, 3)
        quant_err1=QuantumError(err1).to_instruction()
        quant_err2=QuantumError(err2).to_instruction()
        quant_err3=QuantumError(err3).to_instruction()
        quant_errs=[quant_err1, quant_err2, quant_err3]
        num_lqubits= len(circ.qubits)
        num_pqubits=reps*num_lqubits
        phys_qr=QuantumRegister(num_pqubits)
        #ancillas for stabilizer measurements. the number of generators is reps-1. Each generator will be a ZZ. So we need 2*(num_lqubits-1)*num_lqubits
        # in total 
        num_gen=(reps-1)
        phys_ancil_qr1=QuantumRegister(2*num_gen*num_lqubits)
        #ancillas for observable measurement of logical sigma_z.
        phys_ancil_qr2=QuantumRegister(num_lqubits)
        phys_cr=ClassicalRegister(num_pqubits)
        #cr for stabilizer measurements. 
        phys_ancil_cr1=ClassicalRegister(2*num_gen*num_lqubits)
        #cr for the observable measurements.
        phys_ancil_cr2=ClassicalRegister(num_lqubits)

        phys_qc=QuantumCircuit(phys_qr, phys_ancil_qr1, phys_ancil_qr2, phys_cr, phys_ancil_cr2, phys_ancil_cr1)
        logical_qubits=commonUtl.chunkwise(phys_qr, reps)
        #Each generator will be a ZZ so each lqubit requires 2*num_gen physical ancillas.
        logical_ancil_qubits1=commonUtl.chunkwise(phys_ancil_qr1, 2*num_gen)
        logical_cregs=commonUtl.chunkwise(phys_cr, reps)
        logical_ancil_cregs1=commonUtl.chunkwise(phys_ancil_cr1, 2*num_gen)

        Encodings.init_ancillas(phys_qc, logical_ancil_qubits1, quant_errs)

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
                        LOpsWithDepol.cnot_l(phys_qc, logical_qubits[ctrl_idx], logical_qubits[tar_idx], quant_errs)
                    elif node.name=="rz":
                        # print("h")
                        # print(logical_qubits[node.qargs[0].index])
                        # print(node)
                        # print(node.op.params)
                        # print(type(node.op))
                        LOpsWithDepol.rzgate_l(phys_qc, node.op.params[0], logical_qubits[node.qargs[0].index], quant_errs)
                    elif node.name=="sx":
                        # print("s")
                        LOpsWithDepol.sxgate_l(phys_qc, logical_qubits[node.qargs[0].index], quant_errs)
                    elif node.name=="sxdg":
                        # print("sdg")
                        LOpsWithDepol.sxdggate_l(phys_qc, logical_qubits[node.qargs[0].index], quant_errs)
                    elif node.name=="x":
                        # print("t")
                        LOpsWithDepol.xgate_l(phys_qc, logical_qubits[node.qargs[0].index], quant_errs)
                    elif node.name=="measure":
                        # we measure at the end after the stabilizer and the observable measurements.
                        continue
                        # LOpsWithDepol.mgate_l(phys_qc, logical_qubits[node.qargs[0].index], logical_cregs[node.cargs[0].index])
                    else:
                        assert False, f"{node.name} gate not recognized."

        phys_qc.barrier()
        Encodings.meas_stabs(phys_qc, logical_qubits, logical_ancil_qubits1, logical_ancil_cregs1, quant_errs)
        phys_qc.barrier()
        Encodings.meas_sigmaz(phys_qc, logical_qubits, phys_ancil_qr2, phys_ancil_cr2, quant_errs)
        phys_qc.barrier()
        for idx, log_qubit in enumerate(logical_qubits):
            phys_qc.measure(log_qubit, logical_cregs[idx])

        return phys_qc

    @staticmethod
    def to_encoded_circ_CliffordT_ft_depol(circ, reps, one_qerr_prob):
        '''Encoded FT circ with depol added.'''
        # noise_model=NoiseModel()
        err1=depolarizing_error(one_qerr_prob, 1)
        err2=depolarizing_error(10*one_qerr_prob, 2)
        err3=depolarizing_error(5*one_qerr_prob, 3)
        quant_err1=QuantumError(err1).to_instruction()
        quant_err2=QuantumError(err2).to_instruction()
        quant_err3=QuantumError(err3).to_instruction()
        quant_errs=[quant_err1, quant_err2, quant_err3]

        num_lqubits= len(circ.qubits)
        num_pqubits=reps*num_lqubits
        phys_qr=QuantumRegister(num_pqubits)
        phys_cr=ClassicalRegister(num_pqubits)

        #Extra registers for hadamard
        hadamard_plus_qr=QuantumRegister(reps)
        hadamard_minus_qr=QuantumRegister(reps)
        hadamard_cr=ClassicalRegister(reps)

        phys_qc=QuantumCircuit(phys_qr, hadamard_plus_qr, hadamard_minus_qr, phys_cr, hadamard_cr)
        logical_qubits=commonUtl.chunkwise(phys_qr, reps)
        logical_cregs=commonUtl.chunkwise(phys_cr, reps)

        print(logical_qubits)
        qc_dag=circuit_to_dag(circ)
        layers= list(qc_dag.multigraph_layers())
        # Iterate through the layers and convert to the encoding.
        for layer in layers:
            for node in layer:
                if type(node)==DAGOpNode:
                    if node.name=="barrier":
                        phys_qc.barrier()
                    elif node.name=="rz":
                        # print("here")
                        # print(node.qargs[0].index)
                        # print(type(node.qargs[0].index))
                        # print(logical_qubits[node.qargs[0].index])
                        LOpsWithDepol.rzgate_l(phys_qc, node.cargs[0], logical_qubits[node.qargs[0].index], quant_errs)

                    elif node.name=="cx":
                        # print("cx")
                        ctrl_idx=node.qargs[0].index
                        tar_idx=node.qargs[1].index
                        # print(logical_qubits[ctrl_idx])
                        # print(logical_qubits[tar_idx])
                        LOpsWithDepol.cnot_l(phys_qc, logical_qubits[ctrl_idx], logical_qubits[tar_idx], quant_errs)
                    elif node.name=="h":
                        # print("h")
                        # print(logical_qubits[node.qargs[0].index])
                        # hadamard_l(phys_qc, logical_qubits[node.qargs[0].index])
                        # print("target before: ", logical_qubits[node.qargs[0].index])
                        # print("l minus before: ", hadamard_minus_qr)
                        # print("l plus before: ", hadamard_plus_qr)
                        #(circ, qubit_p, qubit_m, qubit_t, hadamard_reg)
                        LOpsWithDepol.hadamard_fault_tol_l(phys_qc, hadamard_plus_qr, hadamard_minus_qr, logical_qubits[node.qargs[0].index], hadamard_cr, quant_errs)
                        temp=deepcopy(logical_qubits[node.qargs[0].index])
                        logical_qubits[node.qargs[0].index]=hadamard_plus_qr
                        hadamard_plus_qr=temp
                        # print("target after: ", logical_qubits[node.qargs[0].index])
                        # print("l minus after: ", hadamard_minus_qr)
                        # print("l plus after: ", hadamard_plus_qr)

                    elif node.name=="s":
                        # print("s")
                        LOpsWithDepol.sgate_l(phys_qc, logical_qubits[node.qargs[0].index], quant_errs)
                    elif node.name=="sdg":
                        # print("sdg")
                        LOpsWithDepol.sdggate_l(phys_qc, logical_qubits[node.qargs[0].index], quant_errs)
                    elif node.name=="t":
                        # print("t")
                        LOpsWithDepol.tgate_l(phys_qc, logical_qubits[node.qargs[0].index], quant_errs)
                    elif node.name=="tdg":
                        # print("t")
                        LOpsWithDepol.tdggate_l(phys_qc, logical_qubits[node.qargs[0].index], quant_errs)
                    elif node.name=="measure":
                        # print("m")
                        # print(node.qargs)
                        # print(node.cargs)
                        LOpsWithDepol.mgate_l(phys_qc, logical_qubits[node.qargs[0].index], logical_cregs[node.cargs[0].index])
                    else:
                        assert False, f"{node.name} gate not recognized."
        return phys_qc

    @staticmethod
    def create_lplus_state(circ, qr):
        q1=qr[0]
        circ.h(q1)
        qr_temp=[elem for elem in qr]
        # print("all qr plus: ", qr_temp)
        for q2 in qr_temp[1:]:
            # print("qr after: ", qr_temp[1:])
            circ.cx(q1,q2)
    @staticmethod
    def create_lminus_state(circ, qr):
        q1=qr[0]
        circ.x(q1)
        circ.h(q1)
        qr_temp=[elem for elem in qr]
        # print("all qr minus: ", qr_temp)
        for q2 in qr_temp[1:]:
            # print("qr after: ", qr_temp[1:])
            circ.cx(q1,q2)

    @staticmethod
    def to_encoded_circ_test_1Qubit_Gate(circ, reps):
        num_lqubits= len(circ.qubits)
        num_pqubits=reps*num_lqubits
        phys_qr=QuantumRegister(num_pqubits)
        phys_cr=ClassicalRegister(num_pqubits)
        phys_qc=QuantumCircuit(phys_qr, phys_cr)
        logical_qubits=commonUtl.chunkwise(phys_qr, reps)
        logical_cregs=commonUtl.chunkwise(phys_cr, reps)

        # print(logical_qubits)
        qc_dag=circuit_to_dag(circ)
        layers= list(qc_dag.multigraph_layers())
        # Iterate through the layers and convert to the encoding.
        for layer in layers:
            for node in layer:
                if type(node)==DAGOpNode:
                    if node.name=="barrier":
                        phys_qc.barrier()

                    elif node.name=="s":
                        # print("s")
                        LOps.sgate_l_reptest(phys_qc, logical_qubits[node.qargs[0].index])
                    elif node.name=="sdg":
                        # print("sdg")
                        LOps.sdggate_l_reptest(phys_qc, logical_qubits[node.qargs[0].index])
                    elif node.name=="t":
                        # print("t")
                        LOps.tgate_l_reptest(phys_qc, logical_qubits[node.qargs[0].index])
                    elif node.name=="measure":
                        # print("m")
                        # print(node.qargs)
                        # print(node.cargs)
                        LOps.mgate_l(phys_qc, logical_qubits[node.qargs[0].index], logical_cregs[node.cargs[0].index])
                    else:
                        assert False, f"{node.name} gate not recognized."
        return phys_qc

    @staticmethod
    def to_encoded_circ_CliffT(circ, reps):
        num_lqubits= len(circ.qubits)
        num_pqubits=reps*num_lqubits
        phys_qr=QuantumRegister(num_pqubits)
        phys_cr=ClassicalRegister(num_pqubits)
        phys_qc=QuantumCircuit(phys_qr, phys_cr)
        logical_qubits=commonUtl.chunkwise(phys_qr, reps)
        logical_cregs=commonUtl.chunkwise(phys_cr, reps)

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
                    elif node.name=="tdg":
                        # print("t")
                        LOps.tdggate_l(phys_qc, logical_qubits[node.qargs[0].index])
                    elif node.name=="x":
                        # print("t")
                        LOps.xgate_l(phys_qc, logical_qubits[node.qargs[0].index])
                    elif node.name=="y":
                        # print("t")
                        LOps.ygate_l(phys_qc, logical_qubits[node.qargs[0].index])
                    elif node.name=="z":
                        # print("t")
                        LOps.zgate_l(phys_qc, logical_qubits[node.qargs[0].index])
                    elif node.name=="measure":
                        # print("m")
                        # print(node.qargs)
                        # print(node.cargs)
                        LOps.mgate_l(phys_qc, logical_qubits[node.qargs[0].index], logical_cregs[node.cargs[0].index])
                    else:
                        assert False, f"{node.name} gate not recognized."
        return phys_qc

    @staticmethod
    def to_encoded_circ_IBMBasis(circ, reps):
        num_lqubits= len(circ.qubits)
        num_pqubits=reps*num_lqubits
        phys_qr=QuantumRegister(num_pqubits)
        phys_cr=ClassicalRegister(num_pqubits)
        phys_qc=QuantumCircuit(phys_qr, phys_cr)
        logical_qubits=commonUtl.chunkwise(phys_qr, reps)
        logical_cregs=commonUtl.chunkwise(phys_cr, reps)

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
                        LOps.mgate_l(phys_qc, logical_qubits[node.qargs[0].index], logical_cregs[node.cargs[0].index])
                    else:
                        assert False, f"{node.name} gate not recognized."
        return phys_qc

    @staticmethod
    def to_encoded_circ_IBMBasis_with_depol(circ, reps, p1_err, p2_err, p3_err):
        err1=depolarizing_error(p1_err, 1)
        err2=depolarizing_error(p2_err, 2)
        err3=depolarizing_error(p3_err, 3)
        quant_err1=QuantumError(err1).to_instruction()
        quant_err2=QuantumError(err2).to_instruction()
        quant_err3=QuantumError(err3).to_instruction()
        quant_errs=[quant_err1, quant_err2, quant_err3]
        num_lqubits= len(circ.qubits)
        num_pqubits=reps*num_lqubits
        phys_qr=QuantumRegister(num_pqubits)
        phys_cr=ClassicalRegister(num_pqubits)
        phys_qc=QuantumCircuit(phys_qr, phys_cr)
        logical_qubits=commonUtl.chunkwise(phys_qr, reps)
        logical_cregs=commonUtl.chunkwise(phys_cr, reps)

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
                        LOpsWithDepol.cnot_l(phys_qc, logical_qubits[ctrl_idx], logical_qubits[tar_idx], quant_errs)
                    elif node.name=="rz":
                        # print("h")
                        # print(logical_qubits[node.qargs[0].index])
                        # print(node)
                        # print(node.op.params)
                        # print(type(node.op))
                        LOpsWithDepol.rzgate_l(phys_qc, node.op.params[0], logical_qubits[node.qargs[0].index], quant_errs)
                    elif node.name=="sx":
                        # print("s")
                        LOpsWithDepol.sxgate_l(phys_qc, logical_qubits[node.qargs[0].index], quant_errs)
                    elif node.name=="sxdg":
                        # print("sdg")
                        LOpsWithDepol.sxdggate_l(phys_qc, logical_qubits[node.qargs[0].index], quant_errs)
                    elif node.name=="x":
                        # print("t")
                        LOpsWithDepol.xgate_l(phys_qc, logical_qubits[node.qargs[0].index], quant_errs)
                    elif node.name=="measure":
                        # print("m")
                        # print(node.qargs)
                        # print(node.cargs)
                        LOpsWithDepol.mgate_l(phys_qc, logical_qubits[node.qargs[0].index], logical_cregs[node.cargs[0].index])
                    else:
                        assert False, f"{node.name} gate not recognized."
        return phys_qc

    
def testcirc_cx(reps):
    n = 2
    circ=QuantumCircuit(n)
    for _ in range(reps):
        circ.cx(0,1)
    circ.measure_all()
    return circ

def testcirc_cx_11init(reps):
    n = 2
    circ=QuantumCircuit(n)
    circ.x([0,1])
    for _ in range(reps):
        circ.cx(0,1)
    circ.measure_all()
    return circ

def testcirc_cx_10init(reps):
    n = 2
    circ=QuantumCircuit(n)
    circ.x([0,1])
    for _ in range(reps):
        circ.cx(0,1)
    circ.measure_all()
    return circ

def testcirc_cx_111init(reps):
    n = 3
    circ=QuantumCircuit(n)
    circ.x([0,1,2,])
    for _ in range(reps):
        circ.cx(0,1)
        circ.cx(1,2)
    circ.measure_all()
    return circ

def testcirc_cx_1111init(reps):
    n = 4
    circ=QuantumCircuit(n)
    circ.x([0,1,2,3])
    for _ in range(reps):
        circ.cx(0,1)
        circ.cx(1,2)
        circ.cx(2,3)
    circ.measure_all()
    return circ


def testcirc_cx_0000init(reps):
    n = 4
    circ=QuantumCircuit(n)
    for _ in range(reps):
        circ.cx(0,1)
        circ.cx(1,2)
        circ.cx(2,3)
    circ.measure_all()
    return circ

def testcirc_h(reps):
    circ=QuantumCircuit(1)
    for _ in range(reps):
        circ.h(0)
    circ.measure_all()
    return circ

def testcirc_s(reps):
    circ=QuantumCircuit(1)
    for _ in range(reps//2):
        circ.s(0)
        circ.sdg(0)
    circ.measure_all()
    return circ

def testcirc_t(reps):
    circ=QuantumCircuit(1)
    for _ in range(reps):
        circ.t(0)
    circ.measure_all()
    return circ