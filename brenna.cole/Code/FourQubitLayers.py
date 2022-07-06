from qiskit import QuantumCircuit
from qiskit.quantum_info.operators import Operator

q_circuits = []
single_gate_layers = []

############################## CCX Layers #########################

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

def makeQC_ccx_4(params):
    """
    Input:
        array of form t,c1,c2
    Returns:
        QuantumCircuit object that has one ccx gate
    """
    qc = QuantumCircuit(4)
    qc.ccx(params[1],params[2],params[0])
    return qc

param_list = createParams_ccx_4()
ccx_layers = []
for qc in param_list:
    ccx_layers.append(makeQC_ccx_4(qc))
for qc in ccx_layers:
    q_circuits.append(qc)
############################### CCX and X Layers ##################
def createParams_ccx_x_4():
    """
    Returns:
        Array of arrays to be used as input for makeQC_ccx_x_4()
    """
    param_list = []
    num_q = 4
    for i in range(num_q):
        pass
    return(param_list)

def makeQC_ccx_x_4():
    """
    Returns:
        QuantumCircuit object that has ccx and x gates
    """
    qc = QuantumCircuit(4)

    return qc


############################# Single CX Layers ####################
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

def makeQC_cx_4(t, c):
    """
    Returns:
        QuantumCircuit object that has one cx gates
    """
    qc = QuantumCircuit(4)
    qc.cx(c,t)
    return qc

single_cx_array = []
params = createParams_cx_4()
for i in params:
    single_cx_array.append(makeQC_cx_4(i))
for qc in single_cx_array:
    q_circuits.append(qc)

############################### CX and CX Layers ##################
def makeQC_2cx_4(t1,c1,t2,c2):
    """
    Returns:
        QuantumCircuit object that has two cx gates
    """
    qc = QuantumCircuit(4)
    qc.cx(c1,t1)
    qc.cx(c2,t2)
    return qc

#TODO - automate this
two_cx_circs = []

two_cx_circs.append(makeQC_2CX(0,1,2,3))
two_cx_circs.append(makeQC_2CX(0,1,3,2))
two_cx_circs.append(makeQC_2CX(0,2,1,3))
two_cx_circs.append(makeQC_2CX(0,2,3,1))
two_cx_circs.append(makeQC_2CX(0,3,1,2))
two_cx_circs.append(makeQC_2CX(0,3,2,1))
two_cx_circs.append(makeQC_2CX(1,0,2,3))
two_cx_circs.append(makeQC_2CX(1,0,3,2))
two_cx_circs.append(makeQC_2CX(1,2,3,0))
two_cx_circs.append(makeQC_2CX(1,3,2,0))
two_cx_circs.append(makeQC_2CX(2,0,3,1))
two_cx_circs.append(makeQC_2CX(2,1,3,0))

for i in two_cx_circs:
    q_circuits.append(i)

############################### CX and X Layers ###################
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

def makeQC_cx_x_4(q_list):
    """
    Parameter: q_list = [t, c, x_array]
    Returns:
        QuantumCircuit object that has one cx gate and one or two x gates
    """
    qc = QuantumCircuit(4)
    qc.cx(q_list[1],q_list[0])
    for x in q_list[2]:
        qc.x(x)
    return qc

cx_x_array = []
params = createParams_cx_x_4()
for i in params:
    cx_x_array.append(makeQC_cx_x_4(i))
for qc in cx_x_array:
    q_circuits.append(qc)

############################## X Only Layers ######################
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

def makeQC_x_4(x_array):
    """
    Returns:
        QuantumCircuit object that has only x gates
    """
    qc = QuantumCircuit(4)
    for x in x_array:
        qc.x(x)
    return qc

x_array = []
for p in params:
    x_array.append(makeQC_x_4(p))
for qc in x_array:
    q_circuits.append(qc)

