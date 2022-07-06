import numpy as np

a = np.matrix([[0],[1]])
print(a)

import numpy as np
from generate_sequences import generate_circs

from qiskit import QuantumCircuit
from qiskit.quantum_info.operators import Operator


a = np.matrix([[0],[1]])
print(a)


x = np.matrix([[0,1],[1,0]])
y = np.matrix([[0,-1j],[1j,0]])
z = np.matrix([[1,0],[0,-1]])

print(x,y,z)


xy = np.matmul(x,y)
print(xy)

xyz = np.matmul(xy,z)
print(xyz)

xz= np.matmul(x,z)
print(xz)

xzy = np.matmul(xz,y)
print(xzy)

print(xzy == xyz)