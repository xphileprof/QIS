import numpy as np
import itertools
from qiskit import QuantumCircuit


from qiskit.quantum_info.operators import Operator

from qiskit.compiler import transpile
from qiskit.transpiler import PassManager
from qiskit import QuantumCircuit
from qiskit.dagcircuit import DAGCircuit
from qiskit.converters import circuit_to_dag
from qiskit.tools.visualization import dag_drawer
from qiskit.converters import dag_to_circuit

def simplifyMatch(qc):
    layer_circs = []

    #1 Remove redundant gates - TODO - make a unique pass for this
    qc_opt = transpile(qc, optimization_level=3)
    return qc_opt

def separateLayers(qc_opt):
    #2 Break simplified circuit into layers
    ## 2a - Transform the optimized circuit into a DAG
    dag_opt = circuit_to_dag(qc_opt)

    ## 2b - Get each layer of the DAG and make a circuit
    layer_circs = []
    i = 0
    for layer in dag_opt.layers():
        print(i)
        i+=1
        subdag = layer['graph']
        q = dag_to_circuit(subdag)
        layer_circs.append(q)

    while (len(layer_circs)<3):
        new_dag = DAGCircuit()
        for qreg in dag_opt.qregs.values():
            new_dag.add_qreg(qreg)
        q = dag_to_circuit(new_dag)
        layer_circs.append(q)

    #3 With the new circuit layers, check if it still matches the criteria
    ## layer_circs is the triplet to pass into checkMatch

    return(layer_circs)

#Takes the layer options and combines into sequences of 3
def combineLayers(layers):
    sequences = []
    for i in itertools.product(layers, repeat=3):
        sequences.append((i))
    return(sequences)

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
        #print("First and second elements are the same")
        return match
    if triplet[1] == triplet[2]:
        #print("Second and third elements are the same")
        return match

    # Check if AB == BA
    BA = A.compose(B)
    AB = B.compose(A)
    if AB == BA:
        #print("AB equals BA - sequence has pairwise commutation")
        return match

    # Check if BC == CB
    BC = C.compose(B)
    CB = B.compose(C)
    if BC == CB:
        #print("BC equals CB - sequence has pairwise commutation")
        return match

    # Check if ABC == BCA
    ABC = BC.compose(A)
    BCA = A.compose(BC)
    if ABC == BCA:
        match = True
        print("Found a match!!!")

    # Check if ABC == CAB
    CAB = AB.compose(C)
    if ABC == CAB:
        match = True
        print("Found a match!!!!!!!")

    return match



def generate3QubitLayers():

    q_circuits = []

    ############################## X Only Layers #############
    xxx = QuantumCircuit(3)
    xxx.x(0)
    xxx.x(1)
    xxx.x(2)
    q_circuits.append((xxx))

    xxi = QuantumCircuit(3)
    xxi.x(0)
    xxi.x(1)
    q_circuits.append(xxi)

    ixx = QuantumCircuit(3)
    ixx.x(1)
    ixx.x(2)
    q_circuits.append(ixx)

    xix = QuantumCircuit(3)
    xix.x(0)
    xix.x(2)
    q_circuits.append(xix)

    xii = QuantumCircuit(3)
    xii.x(0)
    q_circuits.append(xii)

    ixi = QuantumCircuit(3)
    ixi.x(1)
    q_circuits.append(ixi)

    iix = QuantumCircuit(3)
    iix.x(2)
    q_circuits.append(iix)

    ############################### CX and X Layers ################
    tcx = QuantumCircuit(3)
    tcx.cx(1,0)
    tcx.x(2)
    q_circuits.append(tcx)

    txc = QuantumCircuit(3)
    txc.cx(2,0)
    txc.x(1)
    q_circuits.append(txc)

    ctx = QuantumCircuit(3)
    ctx.cx(0,1)
    ctx.x(2)
    q_circuits.append(ctx)

    xtc = QuantumCircuit(3)
    xtc.cx(2,1)
    xtc.x(0)
    q_circuits.append(xtc)

    cxt = QuantumCircuit(3)
    cxt.cx(0,2)
    cxt.x(1)
    q_circuits.append((cxt))

    xct = QuantumCircuit(3)
    xct.cx(1,2)
    xct.x(0)
    q_circuits.append(xct)

    ############################# CNOT Only Layers ############
    tci = QuantumCircuit(3)
    tci.cx(1,0)
    q_circuits.append(tci)

    tic = QuantumCircuit(3)
    tic.cx(2,0)
    q_circuits.append(tic)

    cti = QuantumCircuit(3)
    cti.cx(0,1)
    q_circuits.append(cti)

    itc = QuantumCircuit(3)
    itc.cx(2,1)
    q_circuits.append(itc)

    cit = QuantumCircuit(3)
    cit.cx(0,2)
    q_circuits.append(cit)

    ict = QuantumCircuit(3)
    ict.cx(1,2)
    q_circuits.append(ict)

    ############################## Toffoli Layers #############
    tcc = QuantumCircuit(3)
    tcc.ccx(1,2,0)
    q_circuits.append(tcc)

    ctc = QuantumCircuit(3)
    ctc.ccx(0,2,1)
    q_circuits.append(ctc)

    cct  = QuantumCircuit(3)
    cct.ccx(0,1,2)
    q_circuits.append(cct)

    return(q_circuits)




def generate2QubitLayers():
    q_circuits = []

    ####### X ##########
    xi = QuantumCircuit(2)
    xi.x(0)
    q_circuits.append(xi)

    xx = QuantumCircuit(2)
    xx.x(0)
    xx.x(1)
    q_circuits.append(xx)

    ix = QuantumCircuit(2)
    ix.x(1)
    q_circuits.append(ix)

    ######### CX ############
    cx = QuantumCircuit(2)
    cx.cx(0,1)
    q_circuits.append(cx)

    xc = QuantumCircuit(2)
    xc.cx(1,0)
    q_circuits.append(xc)

    return q_circuits

#Stores a quantum circuit to a file as a QASM string
def storeQASM(qc, filename):
    file = open(filename, "w+")
    file.write(qc.qasm())
    file.close()

#Reads a file containing a QASM string, returns as a quantum circuit
def readQASM(filename):
    file = open(filename, 'r')
    qc_qasm = file.read()
    qc = QuantumCircuit.from_qasm_str(qc_qasm)
    file.close()
    return qc

