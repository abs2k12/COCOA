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
    line_template = ["gridx", "gridy", "offx", "offy"]
    for line in data:    
        line = line.split()
        if line == []: continue
        line_variables = {}
    
        index = 0
        for var_name in line_template:
            line_variables[var_name] = num(string.strip(line[index]))
            index += 1

        yield line_variables
