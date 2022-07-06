import warnings
warnings.filterwarnings('ignore', category=RuntimeWarning,
                        module='qiskit')
import datetime
import numpy as np
import itertools

from qiskit.quantum_info.operators import Operator

from qiskit.compiler import transpile
from qiskit.transpiler import PassManager
from qiskit import QuantumCircuit
from qiskit.dagcircuit import DAGCircuit
from qiskit.converters import circuit_to_dag
from qiskit.tools.visualization import dag_drawer
from qiskit.converters import dag_to_circuit

from qiskit.transpiler.passes import optimization
from qiskit.transpiler.passes.optimization import TemplateOptimization

from createParams5QB import *
from makeQC5QB import *
from seqFunctions import *

print("Hello")

#Set up files to write to
filename = "D:\\Users\\bcole\\Thesis\\ResultsRound2\\Results110-113.txt"
file = open(filename, "w+")
file.write("Hello \n\n")

#Create paramater arrays for each type of circuit
ccx_array = createParams_ccx_5()
x_array = createParams_x_5()
single_cx_array = createParams_cx_single_5()
double_cx_array = createParams_cx_mult_5()
ccx_x_array = createParams_ccx_x_5()
s_cx_x_array = createParams_cx_x_5()[0]
d_cx_x_array = createParams_cx_x_5()[1]
ccx_cx_array = createParams_ccx_cx_5()

# Create QuantumCircuit objects for each layer
layers = []
for l in ccx_array:
    layers.append(makeQC_ccx_5(l))
for l in x_array:
    layers.append(makeQC_x_5(l))
for l in single_cx_array:
    layers.append(makeQC_single_cx_5(l))
for l in double_cx_array:
    layers.append(makeQC_double_cx_5(l))
for l in ccx_x_array:
    layers.append(makeQC_ccx_x_5(l))
for l in s_cx_x_array:
    layers.append(makeQC_s_cx_x_5(l))
for l in d_cx_x_array:
    layers.append(makeQC_d_cx_x_5(l))
for l in ccx_cx_array:
    layers.append(makeQC_ccx_cx_5(l))
layer_len = len(layers)
print("Number of layers is: \n")
print(len(layers), end='\n\n')
file.write("Number of layers is:  " + str(layer_len) + '\n\n')


#Combine the layers into all possible sequences
file.write("Began combining layers at %s \n" %datetime.datetime.now())
print(("Combining Layers beginning at %s " %datetime.datetime.now()))

sequences = combineLayers(layers)

print(("Finished combining layers at %s " %datetime.datetime.now()))
print("Number of sequences to check is: ")
print(len(sequences))

file.write("Finished Combining Layers at %s \n" %datetime.datetime.now())
file.write("Number of sequences to check is: %s \n\n" %len(sequences))


# Check each sequence to see if it matches the triplet commutation property
# This will take a long time
file.write("Began checking original matches at %s \n" %datetime.datetime.now())
print(("Began checking original matches at %s " %datetime.datetime.now()))

matches = []
i = 110000000
for triplet in sequences[110000000:113000000]:
    print(i, end=' - ')
    i = i + 1
    if (checkMatch(triplet)):
        print("Found a match!")
        matches.append(triplet)
    print("Thus far have ", end="")
    print(len(matches), end=" matches\n")

print(("Finished checking original matches at %s " %datetime.datetime.now()))
print("Total number of matches is: ")
print(len(matches), end='\n\n')
file.write("Finished checking original matches at %s \n" %datetime.datetime.now())
file.write("Number of matches is: %s \n\n" %len(matches))


#Combine tuple into single circuit
print("Combining triplets into single circuit ")
matched_circs = []
for match in matches:
    qc_new = match[0] + match[1] + match[2]
    matched_circs.append(qc_new)
print("Number of Elements in matched_circs is: ")
print(len(matched_circs))


#Write matching circuits to a QASM string
#Store the  (non-processed) matching circuits in matched_circs to QASM strings
print("Writing circuit to QASM strings and storing as files")
i = 0
for qc in matched_circs:
    filename2 = "D:\\Users\\bcole\\Thesis\\FiveQBMatchesRound2_Orig\\check110_113\\match" + str(i) + ".txt"
    file2 = open(filename2, "w+")
    file2.write(qc.qasm())
    file2.close()
    i = i + 1

# Process the matching sequences by removing redundant gates
# Remove rendundant gates in each match with TemplatOptimization
print("Beginning to remove redundant gates at %s " %datetime.datetime.now())
pass_ = TemplateOptimization()
pm = PassManager(pass_)
simplified_matches = []
for qc in matched_circs:
    simplified_matches.append(pm.run(qc))

print("Number of simplified matches is: ")
print(len(simplified_matches))
print("Ending removing redundant gates at %s \n" %datetime.datetime.now())


#Break simplified circuits into layers
print("Breaking simplified circuits into layers at %s " %datetime.datetime.now())
simplified_layers = []
for qc in simplified_matches:
    simplified_layers.append((separateLayers(qc)))

print("Ended breaking simplified circuits into layers at %s " %datetime.datetime.now())
print("Total number of tuples in simplified_layers is: ", end="")
print(len(simplified_layers), end="\n\n")



#Run simplified layers through checkMatch
print("Running simplified circuits into checkMatch at %s " %datetime.datetime.now())
file.write("Running simplified circuits through checkMatch at %s \n" %datetime.datetime.now())

remaining_matches = []
for q in simplified_layers:
    if(checkMatch(q)):
        remaining_matches.append(q)
print("Total number of tuples in remaining matches is: ", end="")
print(len(remaining_matches))
file.write("Finished running simplified circuits through checkMatch at %s \n" %datetime.datetime.now())
file.write("Number of reduced matches is: %s \n\n" %len(remaining_matches))


#Combine reduced tuples into single circuit
print("Combining reduced triplets into single circuit at %s " %datetime.datetime.now())

reduced_circs = []
for match in remaining_matches:
    qc_new = match[0] + match[1] + match[2]
    reduced_circs.append(qc_new)
print("Number of Elements in reduced_circs is: ")
print(len(reduced_circs), end='\n\n')

#Write matching circuits to a QASM string
#Store the reduced and still valid matching circuits in matched_circs to QASM strings
print("Writing circuit to QASM strings and storing as files")
i = 0
for qc in reduced_circs:
    filename2 = "D:\\Users\\bcole\\Thesis\\FiveQBMatchesRound2_Reduced\\check110_113\\match" + str(i) + ".txt"
    file2 = open(filename2, "w+")
    file2.write(qc.qasm())
    file2.close()
    i = i + 1



file.close()

