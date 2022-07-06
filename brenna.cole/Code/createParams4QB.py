from qiskit import QuantumCircuit
from qiskit.quantum_info.operators import Operator


def createParams_ccx_4():
    num_q = 4
    param_list = []

    for i in range(num_q):
        for j in range(num_q):
            if j != i:
                for k in range(j, num_q):
                    if k != i and k != j:
                        param_list.append([i, j, k])
    return (param_list)


# def createParams_ccx_x_4():
#     """
#     Returns:
#         Array of arrays to be used as input for makeQC_ccx_x_4()
#     """
#     param_list = []
#     num_q = 4
#     for i in range(num_q):
#         pass
#     return(param_list)

def createParams_ccx_x_4():
    """
    Returns:
        Array of arrays to be used as input for makeQC_ccx_x_4()
        The first three ints represent the ccx gate (target, control1, control2).
        The remaining ints represent the x gate
    """
    new = []
    ccx_array = createParams_ccx_4()
    x_array = createParams_x_4()
    for ccx in ccx_array:
        for x in x_array:
            match = False
            for qb in ccx:
                if qb in x:
                    match = True
                    break
            if (match == False):
                n = ccx + x
                new.append(n)
    return (new)


def createParams_cx_4():
    """
    Returns:
        Array of arrays to be used as input for makeQC_cx_4()
        Each array contains two integers correlating to control qb and target qb
    """
    param_list = []
    num_q = 4
    for i in range(num_q):
        for j in range(num_q):
            if j != i:
                param_list.append([i, j])
    return(param_list)


# ##CX AND CX
# def createParams_cx_4():
#     """
#     Returns:
#         Array of arrays to be used as input for makeQC_cx_4()
#         There will either be 2 or 4 ints.
#         They correspond to the cx gate as [target, control] or [target1, control1, target2, control2]
#     """
#     single_cx_array = []
#     new = []
#     num_q = 5
#     for i in range(num_q):
#         for j in range(num_q):
#             if j != i:
#                 single_cx_array.append([i, j])
#
#     for cx1 in single_cx_array:
#         for cx2 in single_cx_array:
#             for qb in cx1:
#                 match = False
#                 if qb in cx2:
#                     match = True
#                     break
#             if (match == False):
#                 n = cx1 + cx2
#                 new.append(n)
#     cx_array = single_cx_array + new
#
#     return (cx_array)

def createParams_cx_mult_4():
    """
    Returns:
        Array of arrays to be used as input for makeQC_cx_mult_4()
        There will be 4 ints.
        They correspond to the cx gates as [target1, control1, target2, control2]
    """
    single_cx_array = createParams_cx_4()
    new = []
    num_q = 4

    for cx1 in single_cx_array:
        for cx2 in single_cx_array:
            for qb in cx1:
                match = False
                if qb in cx2:
                    match = True
                    break
            if (match == False):
                check = cx2 +cx1
                if check in new:
                    continue
                n = cx1 + cx2
                new.append(n)
    return (new)


def createParams_cx_x_4():
    """
    Returns:
        Array of integers [t,c,[x1,x2]]
    """
    num_qubits = 4
    param_list = []
    #Target
    for i in range(num_qubits):
        #Control
        for j in range(num_qubits):
            if j != i:
                #X1
                for k in range(num_qubits):
                    if k != j and k !=i:
                        #X2
                        for l in range(num_qubits):
                            if l != i and l != j and l != k:
                                #Append to param_list
                                param_list.append([i,j,[k]])
                                param_list.append([i,j,[l]])
                                param_list.append([i,j,[k,l]])
                        break
    return param_list

def createParams_x_4():
    """
    Returns:
        Array of arrays to be used as input for makeQC_x_4()
    """
    param_list = []
    num_q = 4
    for i in range(num_q):
        param_list.append([i])
        for j in range(i + 1, num_q):
            param_list.append([i, j])
            for k in range(j + 1, num_q):
                param_list.append([i, j, k])
                for l in range(k + 1, num_q):
                    param_list.append([i, j, k, l])
    return(param_list)


