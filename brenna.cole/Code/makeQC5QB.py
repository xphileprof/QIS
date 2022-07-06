from qiskit import QuantumCircuit

def makeQC_ccx_5(params):
    """
    Input:
        array of ints corresponding to t,c1,c2
    Returns:
        QuantumCircuit object that has one ccx gate
    """
    qc = QuantumCircuit(5)
    qc.ccx(params[1],params[2],params[0])
    return qc

def makeQC_x_5(x_array):
    """
    Input:
        array of ints corresponding to lines a x gate is applied to
    Returns:
        QuantumCircuit object that has only x gates
    """
    qc = QuantumCircuit(5)
    for x in x_array:
        qc.x(x)
    return qc

def makeQC_single_cx_5(params):
    """
    Input:
        array of ints corresponding to [t,c]
    Returns:
        QuantumCircuit object that has one cx gate
    """
    qc = QuantumCircuit(5)
    qc.cx(params[1],params[0])
    return qc

def makeQC_double_cx_5(params):
    """
    Input:
        array of ints corresponding to [t1,c1, t2, c2]
    Returns:
        QuantumCircuit object that has two cx gates
    """
    qc = QuantumCircuit(5)
    qc.cx(params[1],params[0])
    qc.cx(params[3], params[2])
    return qc

def makeQC_ccx_x_5(params):
    """
    Input:
        array of ints corresponding to t,c1,c2, x, x(opt)
    Returns:
        QuantumCircuit object that has one ccx gate and one or two x gates
    """
    qc = QuantumCircuit(5)
    qc.ccx(params[1],params[2],params[0])
    qc.x(params[3])
    if (len(params) == 5):
        qc.x(params[4])
    return qc

def makeQC_s_cx_x_5(params):
    """
    Input:
        array of ints corresponding to [t,c, x, x(opt), x(opt)]
    Returns:
        QuantumCircuit object that has one cx gate and 0-3 x gates
    """
    qc = QuantumCircuit(5)
    qc.cx(params[1],params[0])
    x_lines = params[2:]
    for x in x_lines:
        qc.x(x)
    return qc

def makeQC_d_cx_x_5(params):
    """
    Input:
        array of ints corresponding to [t1,c1, t2, c2, x]
    Returns:
        QuantumCircuit object that has two cx gates and one x gate
    """
    qc = QuantumCircuit(5)
    qc.cx(params[1],params[0])
    qc.cx(params[3], params[2])
    qc.x(params[4])
    return qc

def makeQC_ccx_cx_5(params):
    """
    Input:
        array of ints corresponding to t1, c1a, c1b, t2, c2
    Returns:
        QuantumCircuit object that has one ccx gate and one cx gate
    """
    qc = QuantumCircuit(5)
    qc.ccx(params[1],params[2],params[0])
    qc.cx(params[4], params[3])
    return qc