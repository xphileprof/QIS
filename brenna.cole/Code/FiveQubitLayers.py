from qiskit import QuantumCircuit
from qiskit.quantum_info.operators import Operator
import numpy as np
import itertools

############################## CCX Layers #########################

def createParams_ccx_5():
    num_q = 5
    param_list = []

    for i in range(num_q):
        for j in range(num_q):
            if j != i:
                for k in range(j, num_q):
                    if k != i and k != j:
                        param_list.append([i, j, k])
    return (param_list)

# def makeQC_ccx_4(params):
#     """
#     Input:
#         array of form t,c1,c2
#     Returns:
#         QuantumCircuit object that has one ccx gate
#     """
#     qc = QuantumCircuit(4)
#     qc.ccx(params[1],params[2],params[0])
#     return qc
#
# param_list = createParams_ccx_4()
# ccx_layers = []
# for qc in param_list:
#     ccx_layers.append(makeQC_ccx_4(qc))
# for qc in ccx_layers:
#     q_circuits.append(qc)



############################## X Only Layers ######################
def createParams_x_5():
    """
    Returns:
        Array of arrays to be used as input for makeQC_x_4()
    """
    param_list = []
    num_q = 5
    for i in range(num_q):
        param_list.append([i])
        for j in range(i + 1, num_q):
            param_list.append([i, j])
            for k in range(j + 1, num_q):
                param_list.append([i, j, k])
                for l in range(k + 1, num_q):
                    param_list.append([i, j, k, l])
                    for m in range(l + 1, num_q):
                        param_list.append([i, j, k, l, m])
    return(param_list)

# def makeQC_x_4(x_array):
#     """
#     Returns:
#         QuantumCircuit object that has only x gates
#     """
#     qc = QuantumCircuit(4)
#     for x in x_array:
#         qc.x(x)
#     return qc
#
# x_array = []
# for p in params:
#     x_array.append(makeQC_x_4(p))
# for qc in x_array:
#    q_circuits.append(qc)

###################################### CCX and X Layers ###########################################
def createParams_ccx_x_5():
    """
    Returns:
        Array of arrays to be used as input for makeQC_ccx_x_4()
        The first three ints represent the ccx gate (target, control1, control2).
        The remaining ints represent the x gates
    """
    ccx_array = createParams_ccx_5()
    x_array = createParams_x_5()
    new = []
    for ccx in ccx_array:
        for x in x_array:
            match = False
            for qb in ccx:
                if qb in x:
                    print('breaking')
                    match = True
                    break

            print(match)
            if (match == False):
                n = ccx + x
                new.append(n)
    return (new)



############################# CX Only Layers ###########################
def createParams_cx_5():
    """
    Returns:
        Array of arrays to be used as input for makeQC_cx_5()
        There will either be 2 or 4 ints.
        They correspond to the cx gate as [target, control] or [target1, control1, target2, control2]
    """
    single_cx_array = []
    new = []
    num_q = 5
    for i in range(num_q):
        for j in range(num_q):
            if j != i:
                single_cx_array.append([i, j])

    for cx1 in single_cx_array:
        for cx2 in single_cx_array:
            for qb in cx1:
                match = False
                if qb in cx2:
                    match = True
                    break
            if (match == False):
                n = cx1 + cx2
                new.append(n)
    cx_array = single_cx_array + new

    return (cx_array)

#################################### CX and X Layers ##################################
def createParams_cx_x_5():
    """
    Returns:
        Tuple - First element is arrays with one cx and 1-3 x gates.
                - Second element is arrays with two cx and 1 x gates
        They correspond to the cx gate as [[target1, control1, x1, x2, x3], [target1, control1, target2, control2, x]]
    """
    cx_array_singles = createParams_cx_single_5()
    x_array = createParams_x_5()
    cx_array_mult = createParams_cx_mult_5()

    cx1_x = []
    cx2_x = []

    for cx in cx_array_singles:
        for x in x_array:
            match = False
            for qb in cx:
                if qb in x:
                    match = True
                    break

            if (match == False):
                n = cx + x
                cx1_x.append(n)

    for cx in cx_array_mult:
        for x in x_array:
            match = False
            for qb in cx:
                if qb in x:
                    match = True
                    break

            if (match == False):
                n = cx + x
                cx2_x.append(n)

    return ([cx1_x, cx2_x]) 

