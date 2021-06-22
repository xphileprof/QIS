import ToricCircuit
import pandas as pd
import concurrent.futures
import defs
#use to generate whatever data I need


if __name__ == "__main__":
    print(defs.get_all_syndrome(2,3))
    """
    def format_syndrome(syndrome, circuit):
        new_syn = []
        for syn in syndrome:
            new_syn.append(circuit.get_lattice_syn(syn))
        return new_syn


    def make_csv_files(name):
        #df = pd.DataFrame(columns=["X0", "Z1", "X2", "Z3", "Z4", "X5", "X6", "Z7"])
        #df = pd.DataFrame(columns=["X0","Z1","X2","X3","Z4","X5","Z6","Z7","X8","Z9","X10","Z11","Z12","X13","Z14","X15","Z16","Z17","X18","X19","Z20","X21","X22","Z23"])
        df=pd.DataFrame(columns=["X0", "Z1", "X2", "X3", "Z4", "X5", "X6", "Z7", "X8", "Z9", "Z10", "X11", "Z12", "X13", "Z14", "X15", "Z16", "Z17", "X18", "Z19","X20", "Z21", "X22", "Z23", "Z24", "X25", "Z26", "X27", "Z28", "X29", "Z30", "Z31", "X32", "Z33", "X34", "Z35", "X36", "Z37", "Z38", "X39", "X40", "Z41", "X42", "X43", "Z44", "X45", "X46", "Z47"])
        df2 = pd.DataFrame()
        df2_errs = pd.DataFrame()
        for i in range (100):
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


    circuit = ToricCircuit.ToricCircuit(7)
    count = 0
    while():"""


