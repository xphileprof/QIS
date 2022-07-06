def createParams_x_5():
    """
    Returns:
        Array of arrays to be used as input for makeQC_x_4()
    """
    param_list = []
    num_q = 5
    for i in range(num_q):
        param_list.append([i])
        for j in range(i + 1, num_q):
            param_list.append([i, j])
            for k in range(j + 1, num_q):
                param_list.append([i, j, k])
                for l in range(k + 1, num_q):
                    param_list.append([i, j, k, l])
                    for m in range(l + 1, num_q):
                        param_list.append([i, j, k, l, m])
    return (param_list)


def createParams_ccx_5():
    num_q = 5
    param_list = []

    for i in range(num_q):
        for j in range(num_q):
            if j != i:
                for k in range(j, num_q):
                    if k != i and k != j:
                        param_list.append([i, j, k])
    return (param_list)


def createParams_ccx_x_5():
    """
    Returns:
        Array of arrays to be used as input for makeQC_ccx_x_4()
        The first three ints represent the ccx gate (target, control1, control2).
        The remaining ints represent the x gates
    """
    new = []
    ccx_array = createParams_ccx_5()
    x_array = createParams_x_5()
    for ccx in ccx_array:
        for x in x_array:
            match = False
            for qb in ccx:
                if qb in x:
                    match = True
                    break
            if (match == False):
                n = ccx + x
                new.append(n)
    return (new)


def createParams_cx_mult_5():
    """
    Returns:
        Array of arrays to be used as input for makeQC_cx_5()
        There will either be 2 or 4 ints.
        They correspond to the cx gate as [target, control] or [target1, control1, target2, control2]
    """
    single_cx_array = createParams_cx_single_5()
    new = []
    num_q = 5

    for cx1 in single_cx_array:
        for cx2 in single_cx_array:
            for qb in cx1:
                match = False
                if qb in cx2:
                    match = True
                    break
            if (match == False):
                check = cx2 + cx1
                if check in new:
                    continue
                n = cx1 + cx2
                new.append(n)
    return (new)


def createParams_cx_single_5():
    """
    Returns:
        Array of arrays to be used as input for makeQC_cx_5()
        Note - makeQC fct may change depending on num of cx gates
        There will either be 2ints.
        They correspond to the cx gate as [target, control]
    """
    single_cx_array = []
    num_q = 5
    for i in range(num_q):
        for j in range(num_q):
            if j != i:
                single_cx_array.append([i, j])
    return (single_cx_array)


def createParams_cx_x_5():
    """
    Returns:
        Tuple - First element is arrays with one cx and 1-3 x gates.
                - Second element is arrays with two cx and 1 x gates
        They correspond to the cx gate as [[target1, control1, x1, x2, x3], [target1, control1, target2, control2, x]]
    """
    cx_array_singles = createParams_cx_single_5()
    x_array = createParams_x_5()
    cx_array_mult = createParams_cx_mult_5()

    cx1_x = []
    cx2_x = []

    for cx in cx_array_singles:
        for x in x_array:
            match = False
            for qb in cx:
                if qb in x:
                    match = True
                    break

            if (match == False):
                n = cx + x
                cx1_x.append(n)

    for cx in cx_array_mult:
        for x in x_array:
            match = False
            for qb in cx:
                if qb in x:
                    match = True
                    break

            if (match == False):
                n = cx + x
                cx2_x.append(n)

    return ([cx1_x, cx2_x])


def createParams_ccx_cx_5():
    """
    Returns:
        Array of arrays to be used as input for makeQC_ccx_cx_5()
        The first three ints represent the ccx gate (target, control1, control2).
        The 4th and 5th ints the cx gates (target, control)
    """
    new = []
    ccx_array = createParams_ccx_5()
    cx_array = createParams_cx_single_5()
    for ccx in ccx_array:
        for cx in cx_array:
            match = False
            for qb in ccx:
                if qb in cx:
                    match = True
                    break

            if (match == False):
                n = ccx + cx
                new.append(n)
    return (new)