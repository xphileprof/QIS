from qiskit import QuantumCircuit

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

def makeQC_x_4(x_array):
    """
    Returns:
        QuantumCircuit object that has only x gates
    """
    qc = QuantumCircuit(4)
    for x in x_array:
        qc.x(x)
    return qc


def makeQC_cx_4(params):
    """
    Returns:
        QuantumCircuit object that has one cx gates
    """
    qc = QuantumCircuit(4)
    qc.cx(params[1],params[0])
    return qc

def makeQC_2cx_4(params):
    """
    Input:
        List of format [t1,c1,t2,c2]
    Returns:
        QuantumCircuit object that has two cx gates
    """
    qc = QuantumCircuit(4)
    qc.cx(params[1],params[0])
    qc.cx(params[3], params[2])
    return qc

def makeQC_ccx_x_4(params):
    """
    Input:
        List of form [t, c1, c2, x]
    Returns:
        QuantumCircuit object that has ccx and x gates
    """
    qc = QuantumCircuit(4)
    qc.ccx(params[1],params[2],params[0])
    qc.x(params[3])
    return qc

def makeQC_cx_x_4(q_list):
    """
    Input:
        q_list = [t, c, x_array]. x_array may have one or two ints
    Returns:
        QuantumCircuit object that has one cx gate and one or two x gates
    """
    qc = QuantumCircuit(4)
    qc.cx(q_list[1],q_list[0])
    for x in q_list[2]:
        qc.x(x)
    return qc