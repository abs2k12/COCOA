import os
import string
import threading
import exceptions

def num (s):
    try:
        return int(s)
    except exceptions.ValueError:
        return float(s)

def read_lines(data):
    line_template = ["im", "r", "vr", "vt", "idd1", "idd2", "ikb", "ik1", "ik2", "sm1", "sm2", "slum1", "slum2", 
                    "rad1", "rad2", "spin1", "spin2", "iinter1", "iinter2", "a","ecc", "mv", "mbv", "mi", "mv1",
                    "mbv1", "mi1", "mv2", "mbv2", "mi2"]
    for line in data:    
        line = line.split()
        if line == []: continue
        line_variables = {}
    
        index = 0
        for var_name in line_template:
            line_variables[var_name] = num(string.strip(line[index]))
            index += 1

        yield line_variables
