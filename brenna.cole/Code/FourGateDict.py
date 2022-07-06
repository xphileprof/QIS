from qiskit import QuantumCircuit

import numpy as np
from qiskit.quantum_info.operators import Operator


# #Generating arrays for ccx gates with 4 lines
# def generateFourLineMatrices():
#     for i in range(3):
#         print(i)
#         for j in range(3):
#             if j != i:
#                 print(j)
#             for k in range(3):
#                 if k != i and k !=j:
#                     print(k)



#Toffoli gates on 4 lines
t012 = QuantumCircuit(4)
t012.ccx(0,1,2)
opT012 = Operator(t012).data.real

t013 = QuantumCircuit(4)
t013.ccx(0,1,3)
opT013 = Operator(t013).data.real

t023 = QuantumCircuit(4)
t023.ccx(0,2,3)
opT023 = Operator(t023).data.real

t021 = QuantumCircuit(4)
t021.ccx(0,2,1)
opT021 = Operator(t021).data.real

t031 = QuantumCircuit(4)
t031.ccx(0,3,1)
opT031 = Operator(t031).data.real

t032 = QuantumCircuit(4)
t032.ccx(0,3,2)
opT032 = Operator(t032).data.real

t123 = QuantumCircuit(4)
t123.ccx(1,2,3)
opT123 = Operator(t123).data.real

t120 = QuantumCircuit(4)
t120.ccx(1,2,0)
opT120 = Operator(t120).data.real

t130 = QuantumCircuit(4)
t130.ccx(1,3,0)
opT130 = Operator(t130).data.real

t132 = QuantumCircuit(4)
t132.ccx(1,3,2)
opT132 = Operator(t132).data.real

t230 = QuantumCircuit(4)
t230.ccx(2,3,0)
opT230 = Operator(t230).data.real

t231 = QuantumCircuit(4)
t231.ccx(2,3,1)
opT231 = Operator(t231).data.real

#X gates on 4 lines
x0 = QuantumCircuit(4)
x0.x(0)
opX0 = Operator(x0).data.real

x1 = QuantumCircuit(4)
x1.x(1)
opX1 = Operator(x1).data.real

x2 = QuantumCircuit(4)
x2.x(2)
opX2 = Operator(x2).data.real

x3 = QuantumCircuit(4)
x3.x(3)
opX3 = Operator(x3).data.real

#CX Gates on 4 Lines
cx01 = QuantumCircuit(4)
cx01.cx(0,1)
opCX01 = Operator(cx01).data.real

cx02 = QuantumCircuit(4)
cx02.cx(0,2)
opCX02 = Operator(cx02).data.real

cx03 = QuantumCircuit(4)
cx03.cx(0,3)
opCX03 = Operator(cx03).data.real

cx10 = QuantumCircuit(4)
cx10.cx(1,0)
opCX10 = Operator(cx10).data.real

cx12 = QuantumCircuit(4)
cx12.cx(1,2)
opCX12 = Operator(cx12).data.real

cx13 = QuantumCircuit(4)
cx13.cx(1,3)
opCX13 = Operator(cx13).data.real

cx20 = QuantumCircuit(4)
cx20.cx(2,0)
opCX20 = Operator(cx20).data.real

cx21 = QuantumCircuit(4)
cx21.cx(2,1)
opCX21 = Operator(cx21).data.real

cx23 = QuantumCircuit(4)
cx23.cx(2,0)
opCX23 = Operator(cx23).data.real

cx30 = QuantumCircuit(4)
cx30.cx(3,0)
opCX30 = Operator(cx30).data.real

cx31 = QuantumCircuit(4)
cx31.cx(3,1)
opCX31 = Operator(cx31).data.real

cx32 = QuantumCircuit(4)
cx32.cx(3,2)
opCX32 = Operator(cx32).data.real


gate_dictionary_4 = { "t012":opT012,
                      "t013":opT013,
                      "t021":opT021,
                      "t023":opT023,
                      "t031":opT031,
                      "t032":opT032,
                      "t123":opT123,
                      "t120":opT120,
                      "t130":opT130,
                      "t132":opT132,
                      "t230":opT230,
                      "t231":opT231,
                      "x0":opX0,
                      "x1":opX1,
                      "x2":opX2,
                      "x3":opX3,
                      "cx01":opCX01,
                      "cx02":opCX02,
                      "cx03":opCX03,
                      "cx10":opCX10,
                      "cx12":opCX12,
                      "cx13":opCX13,
                      "cx20":opCX20,
                      "cx21":opCX21,
                      "cx23":opCX23,
                      "cx30":opCX30,
                      "cx31":opCX31,
                      "cx32":opCX32 }


def getFourLineDict():
    return gate_dictionary_4