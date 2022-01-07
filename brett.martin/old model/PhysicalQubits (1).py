


class PhysicalQubit:
#relatively simple, only holds location and if it has undergone an x or z error

    def __init__(self, row, col):
        self.Xvalue = 0
        self.Zvalue = 0
        self.locationRow = row
        self.locationCol = col

    def flip_value(self, gate_type):
        if gate_type == 'x' or gate_type == 'y':
            if self.Xvalue == 1:
                self.Xvalue = 0
            else:
                self.Xvalue = 1

        if gate_type == 'z' or gate_type == 'y':
            if self.Zvalue == 1:
                self.Zvalue = 0
            else:
                self.Zvalue = 1

    def clear(self):
        self.Xvalue = 0
        self.Zvalue = 0

    def return_Xvalue(self):
        return self.Xvalue

    def return_Zvalue(self):
        return self.Zvalue