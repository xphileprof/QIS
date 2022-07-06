from qiskit import QuantumCircuit

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