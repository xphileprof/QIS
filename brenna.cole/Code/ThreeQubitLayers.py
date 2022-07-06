import numpy as np
import itertools
import random
from generate_sequences import generate_circs, generate_CXcircs, generate_1Tof2CXcircs, generate_2Tof1CXcircs, generate_3Xcircs, generate_2Tof1Xcircs, generate_1Tof1CX1Xcircs, generate_3Tof_4Lines, generate_3X_4Lines, generate_1Tof_2X_4Lines, generate_2Tof_1X_4Lines
from generate_sequences import generate_2Tof_1CX_4Lines, generate_1Tof_2CX_4Lines, generate_3CX_4Lines, generate_1X_2CX_4Lines, generate_2X_1CX_4Lines, generate_1Tof_1CX_1X_4Linescircs
from gateDict import getDict, toffoliCommutingPairs
from checkSequence import removeRedundantGates, checkCommutation
#from checkMatch import checkMatch
from FourGateDict import getFourLineDict

from qiskit import QuantumCircuit
from qiskit.quantum_info.operators import Operator

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

#
# for qc in q_circuits:
#     print(qc.draw())
#     print("************************************")
#
# print(len(q_circuits))

def checkMatch(triplet):
    if triplet[0] == triplet[1]:
        print("First and second element are the same")


#Sequences is all possible combinations (including repititions) of three layers . 10648 in total
sequences = []
for i in itertools.product(q_circuits, repeat=3):
    sequences.append((i))

print(Operator(sequences[0][0]))

# for triplet in sequences:
#     checkMatch(triplet)

#Returns True if a match is found
#A match is a triplet where ABC = BCA, or ABC = CAB, but AB != BA and BC !=CB
#Significance of these is there exists commutation, but not pairwise commutation
def checkMatch(triplet):
    match = False

    #Creating operators A B and C from the passed quantum circuits
    A = Operator(triplet[0])
    B = Operator(triplet[1])
    C = Operator(triplet[2])

    #Check if there are pairs of the same gate
    if triplet[0] == triplet[1]:
        print("First and second elements are the same")
        return match
    if triplet[1] == triplet[2]:
        print("Second and third elements are the same")
        return match

    #Check if AB == BA
    BA = A.compose(B)
    AB = B.compose(A)
    if AB == BA:
        print("AB equals BA - sequence has pairwise commutation")
        return match

    #Check if BC == CB
    BC = C.compose(B)
    CB = B.compose(C)
    if BC == CB:
        print("BC equals CB - sequence has pairwise commutation")
        return match

    #Check if ABC == BCA
    ABC = BC.compose(A)
    BCA = A.compose(BC)
    if ABC == BCA:
        match = True
        print("Found a match!!!")

    #Check if ABC == CAB
    CAB = AB.compose(C)
    if ABC == CAB:
        match = True
        print("Found a match!!!!!!!")

    return match


# num = random.randint(0,10647)
# print(num)
# result = checkMatch(sequences[num])
#
# print(result)

matches = []

for triplet in sequences:
    result = checkMatch(triplet)
    print(result)
    if result:
        matches.append(triplet)


print(len(matches))

#This will store the three layer circuit
matched_sequences = []

for match in matches:
    qc_new = match[0] + match[1] + match[2]
    matched_sequences.append(qc_new)

i = 0
for circ in matched_sequences:
    print("Element number ", end=": ")
    print(i)
    print(circ)
    i += 1


# letters = ["a",'b','c']
#
# for i in itertools.product(letters, repeat=3):
#     print(i)

















# gate_dictionary = getFourLineDict()
# # print(gate_dictionary)
# # print(gate_dictionary["ccx0"])
#
# #print(checkCommutation(gate_dictionary["ccx0_ccx1"], gate_dictionary["ccx0"]))
# #
# circ_list = generate_1Tof_1CX_1X_4Linescircs()
# # print(circ_list)
#
# print(len(circ_list))
#
# circ_list = removeRedundantGates(circ_list)
# print(len(circ_list))
# #
# matches = []
#
# for i in circ_list:
#     print("Checking if matches criteria for ")
#     print(i, end=" : \n")
#     result = checkMatch(i, gate_dictionary)
#     if result:
#         matches.append(i)
#     print(checkMatch(i,gate_dictionary), end="\n\n")
# print(len(matches))




