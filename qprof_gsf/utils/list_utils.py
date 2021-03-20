
def list2_to_list(list2):
    """
    input: a list of list
    output: a list
    """

    out_list = []
    for list1 in list2:
        for el in list1:
            out_list.append(el)

    return out_list


def list3_to_list(list3):
    """
    input: a list of list of list
    output: a list
    """

    out_list = []
    for list2 in list3:
        for list1 in list2:
            out_list += list1

    return out_list
