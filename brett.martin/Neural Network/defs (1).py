import ToricCircuit
import pandas as pd
import itertools
import concurrent.futures

def thread_func(num):
    name = str(num) + "Depth7Train.csv"
    make_csv_files(name)

def format_syndrome(syndrome, circuit):
    new_syn = []
    for syn in syndrome:
        new_syn.append(circuit.get_lattice_syn(syn))
    return new_syn

def get_all_syndrome(num_data_errs, depth):
    df = pd.DataFrame(
        columns=["X0", "Z1", "X2", "Z3", "Z4", "X5", "X6", "Z7"])
        #columns=["X0", "Z1", "X2", "X3", "Z4", "X5", "Z6", "Z7", "X8", "Z9", "X10", "Z11", "Z12", "X13", "Z14","X15", "Z16", "Z17", "X18", "X19", "Z20", "X21", "X22", "Z23"])
        #columns=["X0", "Z1", "X2", "X3", "Z4", "X5", "X6", "Z7", "X8", "Z9", "Z10", "X11", "Z12", "X13", "Z14", "X15",
                #"Z16", "Z17", "X18", "Z19", "X20", "Z21", "X22", "Z23", "Z24", "X25", "Z26", "X27", "Z28", "X29",
               #"Z30", "Z31", "X32", "Z33", "X34", "Z35", "X36", "Z37", "Z38", "X39", "X40", "Z41", "X42", "X43",
                #"Z44", "X45", "X46", "Z47"])
    df2 = pd.DataFrame(columns=["XSyn", "ZSyn"])
    xsyn = []
    zsyn = []
    df2_errs = pd.DataFrame()
    combos = []

    for i in range (num_data_errs+1):
        combos = combos + (list(itertools.combinations_with_replacement([1, 2, 3], i)))
    circuit = ToricCircuit.ToricCircuit(depth)

    qubit_list_x = []
    qubit_list_z = []
    qubit_list_y = []
    for x in range(depth):
        for y in range(depth):
            qubit_list_x.append('x' + str(x) + str(y))
            qubit_list_z.append('z' + str(x) + str(y))
            qubit_list_y.append('y' + str(x) + str(y))

    final_errs = []
    for errs in combos: #for each combination
        if len(errs) > 0:
                combos_x = []
                combos_y = []
                combos_z = []

        #for each type of error in the combination, put that error at all possible qubits
                combos_x = list(itertools.combinations(qubit_list_x, errs.count(1)))
                if len(combos_x) < 2:
                    combos_x = []

                combos_z = list(itertools.combinations(qubit_list_z, errs.count(2)))
                if len(combos_z) < 2:
                    combos_z = []

                combos_y = list(itertools.combinations(qubit_list_y, errs.count(3)))
                if len(combos_y) < 2:
                    combos_y = []

                #make product of all the lists
                temp_list = []
                for x_errs, z_errs in itertools.product(combos_x, combos_z):
                    temp_list.append(x_errs+z_errs)
                if len(temp_list) == 0 and len(combos_x) > 0:
                    temp_list = combos_x.copy()
                if len(temp_list) == 0 and len(combos_z) > 0:
                    temp_list = combos_z.copy()

                if len(combos_y) > 0:
                    total_list = []
                    for xz_errs, y_errs in itertools.product(temp_list, combos_y): #all three types of errors
                        total_list.append(xz_errs+y_errs)
                    if len(total_list) == 0: #only y errors
                        total_list = combos_y.copy()
                else: #only x and z errors
                    total_list = temp_list.copy()

                for c in total_list:
                    circuit.add_physical_errs(c)
                    syn = circuit.get_syndrome(0)
                    df = df.append(syn, ignore_index=True)

                    total_syn = circuit.get_syn()
                    xsyn.append(format_syndrome(total_syn[0], circuit))
                    zsyn.append(format_syndrome(total_syn[1], circuit))
                    df2_errs = df2_errs.append(circuit.get_physical_errs(),ignore_index=True)

                    circuit.clear_errors()


        else:
            df.append(circuit.get_syndrome(0), ignore_index=True)
    df.to_csv("depth" + str(depth) + "_all_combos.csv", index=False)
    df2 = pd.DataFrame({"XSyn":xsyn,"ZSyn":zsyn})
    df2 = df2_errs.join(df2)
    df2.to_csv("Graph_depth" + str(depth) + "_all_combos.csv", index=False)
    return df2


def make_csv_files(name):
    df = pd.DataFrame(columns=["X0", "Z1", "X2", "Z3", "Z4", "X5", "X6", "Z7"])
    #df = pd.DataFrame(columns=["X0","Z1","X2","X3","Z4","X5","Z6","Z7","X8","Z9","X10","Z11","Z12","X13","Z14","X15","Z16","Z17","X18","X19","Z20","X21","X22","Z23"])
    #df=pd.DataFrame(columns=["X0", "Z1", "X2", "X3", "Z4", "X5", "X6", "Z7", "X8", "Z9", "Z10", "X11", "Z12", "X13", "Z14", "X15", "Z16", "Z17", "X18", "Z19","X20", "Z21", "X22", "Z23", "Z24", "X25", "Z26", "X27", "Z28", "X29", "Z30", "Z31", "X32", "Z33", "X34", "Z35", "X36", "Z37", "Z38", "X39", "X40", "Z41", "X42", "X43", "Z44", "X45", "X46", "Z47"])
    df2 = pd.DataFrame()
    df2_errs = pd.DataFrame()
    for i in range (100000):
        circuit = ToricCircuit.ToricCircuit(7)
        circuit.add_random_error(3)
        syn = circuit.get_syndrome(0.0)

        total_syn = circuit.get_syn()
        x_syn = total_syn[0]
        z_syn = total_syn[1]
        #df2 = df2.append([format_syndrome(x_syn, circuit), format_syndrome(z_syn, circuit)], ignore_index=True)
        #df2_errs = df2_errs.append(circuit.get_phyiscal_errs(),ignore_index=True)


        df = df.append(syn, ignore_index=True)
        df = df.append(circuit.get_physical_errs(), ignore_index=True)

    df.to_csv(name, index=False)
    #df2 = df2_errs.join(df2)
    #df2.to_csv("Classic_decoding"+name, index=False)
