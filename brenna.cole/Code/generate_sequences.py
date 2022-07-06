from qiskit import QuantumCircuit

import numpy as np
from qiskit.quantum_info.operators import Operator

tof_gates = ["ccx0","ccx1","ccx2"]
cx_gates = ["cx01", "cx02", "cx10", "cx12", "cx20", "cx21"]
x_gates = ["x0", "x1", "x2"]

tof_gates_4 = ["t012", "t013", "t021","t023","t031", "t032","t123","t120","t130","t132","t230","t231"]
cx_gates_4 = ["cx01","cx02","cx03","cx10","cx12","cx13","cx20","cx21","cx23","cx30","cx31","cx32"
]
x_gates_4 = ["x0", "x1", "x2", "x3"]



def generate_circs():
    #returns a list of tuples, each tuple is the sequence of gates for the circuit
    circ_list  = []
    for i in tof_gates:
        for j in tof_gates:
            for k in tof_gates:
                circ = (i,j,k)
                circ_list.append(circ)
    return circ_list

def generate_CXcircs():
    #returns a list of tuples, each tuple is the sequence of gates for the circuit
    circ_list  = []
    for i in cx_gates:
        for j in cx_gates:
            for k in cx_gates:
                circ = (i,j,k)
                circ_list.append(circ)
    return circ_list

def generate_1Tof2CXcircs():
    # returns a list of tuples, each tuple is the sequence of gates for the circuit
    circ_list = []
    for i in tof_gates:
        for j in cx_gates:
            for k in cx_gates:
                circ1 = (i, j, k)
                circ_list.append(circ1)
                circ2 = (j, k, i)
                circ_list.append(circ2)
                circ3 = (k, j ,i )
                circ_list.append(circ3)
    return circ_list

def generate_2Tof1CXcircs():
    # returns a list of tuples, each tuple is the sequence of gates for the circuit
    circ_list = []
    for i in tof_gates:
        for j in tof_gates:
            for k in cx_gates:
                circ1 = (i, j, k)
                circ_list.append(circ1)
                circ2 = (j, k, i)
                circ_list.append(circ2)
                circ3 = (k, j, i)
                circ_list.append(circ3)
    return circ_list

def generate_3Xcircs():
    # returns a list of tuples, each tuple is the sequence of gates for the circuit
    circ_list = []
    for i in x_gates:
        for j in x_gates:
            for k in x_gates:
                circ = (i, j, k)
                circ_list.append(circ)
    return circ_list

def generate_2Tof1Xcircs():
    # returns a list of tuples, each tuple is the sequence of gates for the circuit
    circ_list = []
    for i in tof_gates:
        for j in tof_gates:
            for k in x_gates:
                circ1 = (i, j, k)
                circ_list.append(circ1)
                circ2 = (j, k, i)
                circ_list.append(circ2)
                circ3 = (k, j, i)
                circ_list.append(circ3)
    return circ_list

def generate_1Tof2Xcircs():
    # returns a list of tuples, each tuple is the sequence of gates for the circuit
    circ_list = []
    for i in tof_gates:
        for j in x_gates:
            for k in x_gates:
                circ1 = (i, j, k)
                circ_list.append(circ1)
                circ2 = (j, k, i)
                circ_list.append(circ2)
                circ3 = (k, j, i)
                circ_list.append(circ3)
    return circ_list

def generate_1Tof1CX1Xcircs():
    # returns a list of tuples, each tuple is the sequence of gates for the circuit
    circ_list = []
    for i in tof_gates:
        for j in cx_gates:
            for k in x_gates:
                circ1 = (i, j, k)
                circ_list.append(circ1)
                circ1a = (i, k, j)
                circ_list.append(circ1a)
                circ2 = (j, k, i)
                circ_list.append(circ2)
                circ2a = (j,i,k)
                circ_list.append(circ2a)
                circ3 = (k, j, i)
                circ_list.append(circ3)
                circ3a = (k,i,j)
                circ_list.append(circ3a)
    return circ_list



def generate_1Tof_1CX_1X_4Linescircs():
    # returns a list of tuples, each tuple is the sequence of gates for the circuit
    circ_list = []
    for i in tof_gates_4:
        for j in cx_gates_4:
            for k in x_gates_4:
                circ1 = (i, j, k)
                circ_list.append(circ1)
                circ1a = (i, k, j)
                circ_list.append(circ1a)
                circ2 = (j, k, i)
                circ_list.append(circ2)
                circ2a = (j,i,k)
                circ_list.append(circ2a)
                circ3 = (k, j, i)
                circ_list.append(circ3)
                circ3a = (k,i,j)
                circ_list.append(circ3a)
    return circ_list

def generate_3Tof_4Lines():
    # returns a list of tuples, each tuple is the sequence of gates for the circuit
    circ_list = []
    for i in tof_gates_4:
        for j in tof_gates_4:
            for k in tof_gates_4:
                circ = (i, j, k)
                circ_list.append(circ)
    return circ_list

def generate_3CX_4Lines():
    # returns a list of tuples, each tuple is the sequence of gates for the circuit
    circ_list = []
    for i in cx_gates_4:
        for j in cx_gates_4:
            for k in cx_gates_4:
                circ = (i, j, k)
                circ_list.append(circ)
    return circ_list

def generate_3X_4Lines():
    # returns a list of tuples, each tuple is the sequence of gates for the circuit
    circ_list = []
    for i in x_gates_4:
        for j in x_gates_4:
            for k in x_gates_4:
                circ = (i, j, k)
                circ_list.append(circ)
    return circ_list

def generate_1Tof_2X_4Lines():
    # returns a list of tuples, each tuple is the sequence of gates for the circuit
    circ_list = []
    for i in tof_gates_4:
        for j in x_gates_4:
            for k in x_gates_4:
                circ = (i, j, k)
                circ_list.append(circ)
    return circ_list

def generate_2Tof_1X_4Lines():
    # returns a list of tuples, each tuple is the sequence of gates for the circuit
    circ_list = []
    for i in tof_gates_4:
        for j in tof_gates_4:
            for k in x_gates_4:
                circ = (i, j, k)
                circ_list.append(circ)
    return circ_list


def generate_2Tof_1CX_4Lines():
    # returns a list of tuples, each tuple is the sequence of gates for the circuit
    circ_list = []
    for i in tof_gates_4:
        for j in tof_gates_4:
            for k in cx_gates_4:
                circ = (i, j, k)
                circ_list.append(circ)
    return circ_list

def generate_1Tof_2CX_4Lines():
    # returns a list of tuples, each tuple is the sequence of gates for the circuit
    circ_list = []
    for i in tof_gates_4:
        for j in cx_gates_4:
            for k in cx_gates_4:
                circ = (i, j, k)
                circ_list.append(circ)
    return circ_list

def generate_2X_1CX_4Lines():
    # returns a list of tuples, each tuple is the sequence of gates for the circuit
    circ_list = []
    for i in x_gates_4:
        for j in x_gates_4:
            for k in cx_gates_4:
                circ = (i, j, k)
                circ_list.append(circ)
    return circ_list

def generate_1X_2CX_4Lines():
    # returns a list of tuples, each tuple is the sequence of gates for the circuit
    circ_list = []
    for i in x_gates_4:
        for j in cx_gates_4:
            for k in cx_gates_4:
                circ = (i, j, k)
                circ_list.append(circ)
    return circ_list

    #Choose Gate 1
    #Choose Gate 2
    #Choose Gate 3



