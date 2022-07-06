from qiskit import QuantumCircuit
from qiskit.dagcircuit import DAGCircuit
from qiskit.converters import circuit_to_dag
from qiskit.converters import dag_to_circuit

import itertools

import numpy as np
from qiskit.quantum_info.operators import Operator


# Takes the layer options and combines into sequences of 3
def combineLayers(layers):
    sequences = []
    for i in itertools.product(layers, repeat=3):
        sequences.append((i))
    return (sequences)


# Returns True if a match is found
# A match is a triplet where ABC = BCA, or ABC = CAB, but AB != BA and BC !=CB
# Significance of these is there exists commutation, but not pairwise commutation
def checkMatch(triplet):
    match = False

    # Creating operators A B and C from the passed quantum circuits
    C = Operator(triplet[0])
    B = Operator(triplet[1])
    A = Operator(triplet[2])

    # Check if there are pairs of the same gate
    if triplet[0] == triplet[1]:
        # print("First and second elements are the same")
        return match
    if triplet[1] == triplet[2]:
        # print("Second and third elements are the same")
        return match

    # Check if AB == BA
    BA = A.compose(B)
    AB = B.compose(A)
    if AB == BA:
        # print("AB equals BA - sequence has pairwise commutation")
        return match

    # Check if BC == CB
    BC = C.compose(B)
    CB = B.compose(C)
    if BC == CB:
        # print("BC equals CB - sequence has pairwise commutation")
        return match

    # Check if ABC == BCA
    ABC = BC.compose(A)
    BCA = A.compose(BC)
    if ABC == BCA:
        match = True
    #        print("Found a match!!!")

    # Check if ABC == CAB
    CAB = AB.compose(C)
    if ABC == CAB:
        match = True
    #       print("Found a match!!!!!!!")

    return match


def separateLayers(qc):
    #Input a quantum circuit object with three layers
    #Returns a list of the individual layers as circuits

    #2 Break circuit into layers
    ## 2a - Transform the circuit into a DAG
    dag = circuit_to_dag(qc)

    ## 2b - Get each layer of the DAG and make a circuit
    layer_circs = []
    i = 0
    for layer in dag.layers():
        i+=1
        subdag = layer['graph']
        q = dag_to_circuit(subdag)
        layer_circs.append(q)

    while (len(layer_circs)<3):
        new_dag = DAGCircuit()
        for qreg in dag.qregs.values():
            new_dag.add_qreg(qreg)
        q = dag_to_circuit(new_dag)
        layer_circs.append(q)
    return(layer_circs)