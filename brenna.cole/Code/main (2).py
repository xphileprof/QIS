import warnings
warnings.filterwarnings('ignore', category=RuntimeWarning,
                        module='qiskit')

from qiskit import QuantumCircuit
from qiskit import qasm

#Read circuit in from QASM file
circs_from_qasm = []

i = 0
while i < 343171:
    filename = "D:\\Users\\bcole\\Thesis\\FiveQBMatchesRound2_Orig\\check70_90\\match" + str(i) + ".txt"
    file = open(filename, 'r')
    qc_qasm = file.read()
    qc = QuantumCircuit.from_qasm_str(qc_qasm)
    file.close()
    circs_from_qasm.append(qc)
    print("Read in match %s" %i)
    i = i + 1

#Check for single gate layers
three_gate_circs = []
match_nums = []
i = 0
for circ in circs_from_qasm:
    if circ.size() <= 3:
        three_gate_circs.append(circ)
        match_nums.append(i)
    i = i + 1
#print("Number of three gate circuits : %s " % (len(three_gate_circs)))

#Count ops and write result to the file
filename2 = "D:\\Users\\bcole\\Thesis\\FiveQBOrigOpCounts6.txt"
file2 = open(filename2, "w+")
file2.write("MatchNumber,Number ccx gates,Number cx gates,Number x gates, total num gates\n")

i = 0
for circ in circs_from_qasm:
    ops_dict = circ.count_ops()

    if 'ccx' in ops_dict.keys():
        ccx_num = ops_dict['ccx']
    else:
        ccx_num = 0

    if 'cx' in ops_dict.keys():
        cx_num = ops_dict['cx']
    else:
        cx_num = 0

    if 'x' in ops_dict.keys():
        x_num = ops_dict['x']
    else:
        x_num = 0

    num_gates = circ.size()

    #print("%s,%s,%s,%s \n" %(i, ccx_num, cx_num, x_num))
    file2.write("%s,%s,%s,%s,%s \n" %(i, ccx_num, cx_num, x_num, num_gates))
    print(i)
    i = i + 1

print("Number of circuits with a single gate per layer is: %s " % (len(three_gate_circs)))
if len(three_gate_circs) > 0:
    print(match_nums)
    for circ in three_gate_circs:
        print(circ)

