# FTQEM
Fault Tolerant Quantum Error Mitigation (FTQEM)
https://arxiv.org/abs/2308.05403

FTQEM protocol (works if you are below the error threshold for FTQEM):
(1) encode, (2) perform stabilizer quantum error detection
once at the end, and (3) concatenate to reduce the noise
to arbitrary levels.

This repo contains the simulations and hardware executions for the FTQEM paper.
Tested with:
Python 3.11.2
{'qiskit-terra': '0.23.3', 
'qiskit-aer': '0.12.0', 
'qiskit-ignis': None, 
'qiskit-ibmq-provider': '0.20.2', 
'qiskit': '0.42.1', 
'qiskit-nature': None, 
'qiskit-finance': None, 
'qiskit-optimization': None, 
'qiskit-machine-learning': None}