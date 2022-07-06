import numpy as np
#
#
def checkMatch(triplet):
    match = False

    # Creating operators A B and C from the passed quantum circuits
    C = Operator(triplet[0])
    B = Operator(triplet[1])
    A = Operator(triplet[2])

    # Check if there are pairs of the same gate
    if triplet[0] == triplet[1]:
        print("First and second elements are the same")
        return match
    if triplet[1] == triplet[2]:
        print("Second and third elements are the same")
        return match

    # Check if BA == AB
    BA = A.compose(B)
    AB = B.compose(A)
    if AB == BA:
        print("AB equals BA - sequence has pairwise commutation")
        return match

    # Check if CB == BC
    BC = C.compose(B)
    CB = B.compose(C)
    if BC == CB:
        print("BC equals CB - sequence has pairwise commutation")
        return match

    # Check if ABC == BCA
    ABC = BC.compose(A)
    BCA = A.compose(BC)
    if ABC == BCA:
        match = True
        print("Found a match!!!")

    # Check if ABC == CAB
    CAB = AB.compose(C)
    if ABC == CAB:
        match = True
        print("Found a match!!!!!!!")

    print("AB = ", end="")
    print(AB)
    print("BA = ", end="")
    print(BA)
    print("BC = ", end="")
    print(BC)
    print("CB = ", end="")
    print(CB)
    print("ABC = ", end="")
    print(ABC)
    print("CAB = ", end="")
    print(CAB)
    print("BCA = ", end="")
    print(BCA)

    return match


# def checkMatch(qc, dict):
#     match = False
#
#     #Set A,B,C to the matrix corresponding to the gate
#     A = dict[qc[0]]
#     B = dict[qc[1]]
#     C = dict[qc[2]]
#
#     AB = np.matmul(A,B)
#     BA = np.matmul(B,A)
#
#     BC = np.matmul(B,C)
#     CB = np.matmul(C,B)
#
#     ABC = np.matmul(AB,C)
#     CAB = np.matmul(C, AB)
#     BCA = np.matmul(BC,A)
#
#
#     #Check AB == BA
#     if np.array_equal(AB,BA):
#         print("AB eq BA")
#         return match
#
#     #Check BC == CB
#     elif np.array_equal(BC,CB):
#         print("BC eq CB")
#         return match
#
#     elif np.array_equal(ABC,CAB):
#         print("ABC eq CAB")
#         match = True
#
#     elif np.array_equal(ABC,BCA):
#         print("ABC eq BCA")
#         match = True
#
#     return match