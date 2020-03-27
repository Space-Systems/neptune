#!/usr/bin/env python3

''' 
    Takes a NEPTUNE variance output file and plots the individual components.
''' 

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import argparse
from datetime import datetime as dt

parser = argparse.ArgumentParser(description = __doc__)
parser.add_argument('-i', '--input-file', dest='var_file', required=True, help='The NEPTUNE variances output file to be used.')
args = parser.parse_args()

# read file to pandas data frame
data = pd.read_table(
            args.var_file, 
            comment='#', 
            header=None, 
            sep='\s+', 
            names=['date','time','mjd','rx','ry','rz','vx','vy','vz'], parse_dates=[[0,1]]
       )

data_labels = ['rx','ry','rz','vx','vy','vz']

data[data_labels] = data[data_labels].apply(np.sqrt)
data[data_labels] = data[data_labels].multiply(1000.0)

# strip MJD
data = data[['date_time', 'rx', 'ry', 'rz', 'vx', 'vy', 'vz']]

# now plot
data.plot(x='date_time', subplots=True, sharex=True, title='$1\sigma$ errors (r in m, v in m/s)')
plt.show()
