
# Inspired by rotated surface code of different depths, as described here:https://arxiv.org/pdf/1811.12456.pdf
import PhysicalQubits
import AncillaQubits
from random import *
import numpy as np
from collections import defaultdict, Counter


#representation of the rotated surface code for arbitary depth
#contains all the ancilla and physical qubits
#also includes methods to add error on both the physical qubits and the ancilla
class ToricCircuit:

    def __init__(self, depth):
        self.depth = depth
        self.qubits = {}
        self.ancilla = []

        #these three needed to format observations for MWPM decoding
        self.error_qubits_x = []
        self.error_qubits_z = []
        self.lattice_dict = {}

        #initialize qubits to zero
        for i in range(depth):
            for j in range(depth):
                qubit = PhysicalQubits.PhysicalQubit(i, j)
                self.qubits.update({str(i) + str(j): qubit})

        #set up the ancilla qubits
        self.create_ancilla()

    def add_random_error(self, num):
        # random flips with equal probability of x,z,or y on any qubits
        # could possibly flip already corrupted qubits to fix them
        # choose the number of physical qubits to put errors on
        error_types = ['x', 'y', 'z']
        for i in range(num):
            a = randint(0, self.depth-1)
            b = randint(0, self.depth-1)

            type = choice(error_types)
            self.qubits.get(str(a) + str(b)).flip_value(type)
            if type == "x" or type == "y":
                self.error_qubits_x.append(self.qubits.get(str(a) + str(b)))
            if type == "z" or type == "y":
                self.error_qubits_z.append(self.qubits.get(str(a) + str(b)))

    def get_lattice_syn(self, err):
        return self.lattice_dict[err]

    def get_syn(self):
        x_syn = []
        z_syn = []
        for a in self.ancilla:
            #If there is an ancilla measurement and no error on ancilla, then add it to the syndrome
            if a.measure() == -1 and a.error == False:
                if a.AType == 'x':
                    x_syn.append(a.name)
                else:
                    z_syn.append(a.name)
            #if there is no ancilla but an error on the ancilla, then add it to the syndrome
            if a.measure() == 1 and a.error == True:
                if a.AType == 'x':
                    x_syn.append(a.name)
                else:
                    z_syn.append(a.name)

        return [x_syn, z_syn]

    def get_syndrome(self, prob_err):

        return_dict = {}
        for a in self.ancilla:
            if prob_err > 0:
                rand = random()
                if rand <= prob_err:
                    a.err = not a.err
                    if a.measure() == 1:
                        return_dict.update({a.name:  -1})
                    else:
                        return_dict.update({a.name:  1})
                else:
                    return_dict.update({a.name: a.measure()})
            else:
                return_dict.update({a.name: a.measure()})

        else:
            errs = Counter()
            for x in self.error_qubits_x:
                errs.update({"X" + str(x.locationRow) + str(x.locationCol)})
            for z in self.error_qubits_z:
                errs.update({"Z" + str(z.locationRow) + str(z.locationCol)})

            errors_to_remove = []
            for e in errs: #for decoding, duplicate errors flip back the qubit to its intended state, so remove pairs
                if errs[e] %2 == 0:
                    errors_to_remove.append(e)
            for x in errors_to_remove:
                del errs[x]
            return_dict.update({"Labels": list(errs)})

        return return_dict

    def add_physical_errs(self, err_list):
        for err in err_list:
            if err[0] == 'x':
                self.qubits.get(str(err[1]) + str(err[2])).flip_value('x')
                self.error_qubits_x.append(self.qubits.get(str(err[1]) + str(err[2])))
            elif err[0] == 'z':
                self.qubits.get(str(err[1]) + str(err[2])).flip_value('z')
                self.error_qubits_z.append(self.qubits.get(str(err[1]) + str(err[2])))
            elif err[0] == 'y':
                self.qubits.get(str(err[1]) + str(err[2])).flip_value('x')
                self.error_qubits_x.append(self.qubits.get(str(err[1]) + str(err[2])))
                self.qubits.get(str(err[1]) + str(err[2])).flip_value('z')
                self.error_qubits_z.append(self.qubits.get(str(err[1]) + str(err[2])))
            else:
                print("invalid input in circuit")
        return


    def get_physical_errs(self):
        #get the labels (set of physical qubit errors) for a surface code
        return_dict = {}
        errs = []
        for x in self.error_qubits_x:
            errs.append("X" + str(x.locationRow) + str(x.locationCol))
        for z in self.error_qubits_z:
            errs.append("Z" + str(z.locationRow) + str(z.locationCol))
        return_dict.update({"Labels": errs})

        return return_dict

    def clear_errors(self):
        for err_x in self.error_qubits_x:
            self.qubits.get(str(err_x.locationRow) + str(err_x.locationCol)).clear()

        for err_z in self.error_qubits_z:
            self.qubits.get(str(err_z.locationRow) + str(err_z.locationCol)).clear()

        self.error_qubits_z.clear()
        self.error_qubits_x.clear()

    def correct_errs(self, err_list):
        for err in err_list:
            if err[0] == 'X':
                self.qubits.get(str(err[1]) + str(err[2])).flip_value('x')
                if self.qubits.get(str(err[1]) + str(err[2])) in self.error_qubits_x:
                    self.error_qubits_x.remove(self.qubits.get(str(err[1]) + str(err[2])))
                    print("Corrected qubit X" + str(err[1]) + str(err[2]))
                else:
                    self.error_qubits_x.append(self.qubits.get(str(err[1]) + str(err[2])))
                    print("Incorrectly flipped qubit X" + str(err[1]) + str(err[2]))
            elif err[0] == 'Z':
                self.qubits.get(str(err[1]) + str(err[2])).flip_value('z')
                if self.qubits.get(str(err[1]) + str(err[2])) in self.error_qubits_x:
                    self.error_qubits_x.remove(self.qubits.get(str(err[1]) + str(err[2])))
                    print("Corrected qubit Z" + str(err[1]) + str(err[2]))
                else:
                    self.error_qubits_x.append(self.qubits.get(str(err[1]) + str(err[2])))
                    print("Incorrectly flipped qubit Z" + str(err[1]) + str(err[2]))
            else:
                return

    def check_max_errs(self, max_errors):
        if len(self.error_qubits_x) + len(self.error_qubits_z) > max_errors:
            return False
        else:
            return True


    def create_ancilla(self):
        # This adds the ancilla with their correct physical qubits attached
        # only for depth of 3,5,7 code described here:https://arxiv.org/pdf/1811.12456.pdf
        count = 0

        for i in range(self.depth - 1):

            if i % 2 == 0:
                #start with z ancilla
                j = 0
                while j < self.depth:
                    if i == 0 and j != (self.depth - 1):
                        #ancilla located about the grid
                        newX = AncillaQubits.AncillaQubit("x", "X" + str(count))
                        newX.add_qubit(self.qubits.get(str(i) + str(j)))
                        newX.add_qubit(self.qubits.get(str(i) + str(j + 1)))
                        newX.setLocation(i,j)
                        self.ancilla.append(newX)
                        self.lattice_dict.update({"X" + str(count): (0, i - .5, j + .5)})
                        count += 1
                    if j == (self.depth - 1):
                        #at the right most edge ancilla
                        newZ = AncillaQubits.AncillaQubit("z", "Z" + str(count))
                        newZ.add_qubit(self.qubits.get(str(i) + str(j)))
                        newZ.add_qubit(self.qubits.get(str(i + 1) + str(j)))
                        newZ.setLocation(i+1, j)
                        self.ancilla.append(newZ)
                        self.lattice_dict.update({"Z" + str(count): (0, i + .5, j + .5)})
                        count += 1

                    #those two if statements check to see if the code is at the upper boundary and adds the ancilla at the top
                    else:
                        newZ = AncillaQubits.AncillaQubit("z", "Z" + str(count))
                        self.lattice_dict.update({"Z" + str(count): (0, i + .5, j + .5)})
                        newZ.add_qubit(self.qubits.get(str(i) + str(j))) #upper left corner
                        newZ.add_qubit(self.qubits.get(str(i) + str(j+1))) # upper right corner
                        newZ.add_qubit(self.qubits.get(str(i+1) + str(j)))  #lower left corner
                        newZ.add_qubit(self.qubits.get(str(i+1) + str(j + 1)))  #lower right corner
                        newZ.setLocation(i+1, j)
                        count += 1

                        newX = AncillaQubits.AncillaQubit("x", "X" + str(count))
                        self.lattice_dict.update({"X" + str(count): (0, i + .5, j + 1.5)})
                        newX.add_qubit(self.qubits.get(str(i) + str(j+1)))
                        newX.add_qubit(self.qubits.get(str(i) + str(j + 2)))
                        newX.add_qubit(self.qubits.get(str(i + 1) + str(j+1)))
                        newX.add_qubit(self.qubits.get(str(i + 1) + str(j + 2)))
                        newX.setLocation(i+1, j+1)
                        count += 1

                        self.ancilla.append(newZ)
                        self.ancilla.append(newX)

                    j = j+2
            else:
                #start witth x ancilla, new row
                j = 0
                while j < self.depth-1:
                    if j == 0:
                        #add z ancilla to the far left
                        newZ = AncillaQubits.AncillaQubit("z","Z" + str(count))
                        newZ.add_qubit(self.qubits.get(str(i) + str(j)))
                        newZ.add_qubit(self.qubits.get(str(i + 1) + str(j)))
                        newZ.setLocation(i+1, j)
                        self.ancilla.append(newZ)
                        self.lattice_dict.update({"Z" + str(count): (0, i + .5, j - .5)})
                        count += 1

                    newX = AncillaQubits.AncillaQubit("x", "X" + str(count))
                    newX.add_qubit(self.qubits.get(str(i) + str(j)))
                    newX.add_qubit(self.qubits.get(str(i) + str(j + 1)))
                    newX.add_qubit(self.qubits.get(str(i + 1) + str(j)))
                    newX.add_qubit(self.qubits.get(str(i + 1) + str(j + 1)))
                    newX.setLocation(i+1, j+1)
                    self.ancilla.append(newX)
                    self.lattice_dict.update({"X" + str(count): (0, i + .5, j + .5)})
                    count += 1

                    if i == (self.depth - 2) and j < (self.depth - 1):
                        #checks if this is the last row, if it is, then add the ancilla on the bottom
                        newX = AncillaQubits.AncillaQubit("x", "X" + str(count))
                        newX.add_qubit(self.qubits.get(str(i + 1) + str(j + 1)))
                        newX.add_qubit(self.qubits.get(str(i + 1) + str(j + 2)))
                        newX.setLocation(i+2, j+2)
                        self.ancilla.append(newX)
                        self.lattice_dict.update({"X" + str(count): (0, i + 1.5, j + 1.5)})
                        count += 1

                    newZ = AncillaQubits.AncillaQubit("z", "Z" + str(count))
                    newZ.add_qubit(self.qubits.get(str(i) + str(j + 1)))
                    newZ.add_qubit(self.qubits.get(str(i) + str(j + 2)))
                    newZ.add_qubit(self.qubits.get(str(i + 1) + str(j + 1)))
                    newZ.add_qubit(self.qubits.get(str(i + 1) + str(j + 2)))
                    newZ.setLocation(i+1, j+2)
                    self.lattice_dict.update({"Z" + str(count): (0, i + .5, j + 1.5)})
                    count += 1

                    self.ancilla.append(newZ)

                    j = j + 2
        return
