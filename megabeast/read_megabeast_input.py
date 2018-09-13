import numpy as np
from astropy.table import Table
import ast

import pdb

def read_megabeast_input(input_file):
    """
    Read in the file with MegaBEAST settings
    """
    
    # read everything in
    input_data = Table.read(input_file, format='ascii.no_header', delimiter='=', comment='#')
    
    # parse it
    input_vars = {}

    for i in range(len(input_data)):
        try:
            input_vars[input_data['col1'][i]] = ast.literal_eval(input_data['col2'][i])
        except:
            input_vars[input_data['col1'][i]] = str(input_data['col2'][i])


    # return it
    return input_vars


    
