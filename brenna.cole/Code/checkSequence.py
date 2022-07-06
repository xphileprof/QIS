import numpy as np

#Looking for sequences ABC s.t. ABC = CAB or ABC = BCA, with AB != BA and BC !=BC
def removeRedundantGates(circ_list):
    #Check if A=B or B=C. Remove if so since AB == BA or BC == CB (not calculating since they will compute to the identity)
    # print(circ_list)
    print("Number of elements in original circ list is: ")
    print(len(circ_list), end="\n\n")

    rendundant_gates = []

    for i in circ_list:
        if i[0] == i[1]:
            # print('First and second element are the same, removing...')
            # print(i)
            rendundant_gates.append(i)
            continue
        elif i[1] == i[2]:
            # print('Second and third element are the same, removing...')
            # print(i)
            rendundant_gates.append(i)

    # print(rendundant_gates)

    for i in rendundant_gates:
        circ_list.remove(i)


    # print(circ_list)
    print("Number of elements in updated circ list is: ")
    print(len(circ_list), end="\n\n")

    return circ_list



    # Check if AB == BA

    #Check if BC == CB


def checkCommutation(gate1, gate2):
    #Two gates A and B commute if AB == BA
    #Returns bool
    return(np.array_equal(gate1,gate2))


def checkABeqBA(circ):
    pass