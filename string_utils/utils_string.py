
def clean_string(current_string):

    current_string = current_string.replace("\t", "")
    current_string = current_string.replace("\r", "")
    current_string = current_string.replace("\n", "")

    return current_string

