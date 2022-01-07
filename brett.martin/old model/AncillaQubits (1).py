import PhysicalQubits

#representation of an ancilla qubit
#contains all the physical qubits it is connected
#more important method is to determine which measurement it output, 1 or -1
#is either an x-type ancilla or z-type ancilla
class AncillaQubit:

    def __init__(self, type, name):
        self.AType = type
        self.qubits = []
        self.on = False
        self.error = False
        self.name = name
        self.row = 0
        self.col = 0

    def add_qubit(self, new_qubit):
        self.qubits.append(new_qubit)

    def setLocation(self, row, col):
        self.row = row
        self.col = col

    def connected(self, phys_qubit):
        if phys_qubit in self.qubits:
            return True
        else:
            return False

    #returns the measurement based on the error observed on its physical qubits
    def measure(self):
        total = 0
        if self.AType == "x":
            if len(self.qubits) == 2:
                if self.qubits[0].return_Xvalue() == self.qubits[1].return_Xvalue():
                    total = 1
                else:
                    total = -1
            if len(self.qubits) == 4:
                if (self.qubits[0].return_Xvalue() + self.qubits[1].return_Xvalue() + self.qubits[2].return_Xvalue() + self.qubits[3].return_Xvalue()) % 2 == 0:
                    total = 1
                else:
                    total = -1

        elif self.AType == "z":
            if len(self.qubits) == 2:
                if self.qubits[0].return_Zvalue() == self.qubits[1].return_Zvalue():
                    total += 1
                else:
                    total -= 1
            if len(self.qubits) == 4:
                if (self.qubits[0].return_Zvalue() + self.qubits[1].return_Zvalue() + self.qubits[2].return_Zvalue() + self.qubits[3].return_Zvalue()) % 2 == 0:
                    total = 1
                else:
                    total = -1
        return total