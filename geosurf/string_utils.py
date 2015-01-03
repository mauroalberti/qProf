
def clean_string( str ):
    
    str = str.replace("\t","")
    str = str.replace("\r","")
    str = str.replace("\n","")
    
    return str