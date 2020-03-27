#!/usr/bin/env python3

''' 
    Takes a NEPTUNE kepler elements output file and plots the individual components.
''' 

import matplotlib.pyplot as plt 
import numpy as np
import pandas as pd
import argparse
from datetime import datetime as dt

parser = argparse.ArgumentParser(description = __doc__)
parser.add_argument('-i', '--input-file', dest='var_file', required=True, help='The NEPTUNE kepler elements output file to be used.')
args = parser.parse_args()

# read file to pandas data frame
data = pd.read_table(
            args.var_file, 
            comment='#', 
            header=None, 
            sep='\s+', 
            names=['date','time','mjd','sma','ecc','inc','raan','aop','tran','mean'], parse_dates=[[0,1]]
       )

# strip MJD
data = data[['date_time', 'sma', 'ecc', 'inc', 'raan', 'aop']]

# now plot
data.plot(x='date_time', subplots=True, sharex=True, title='Osculating Kepler elements (km, deg)')
plt.show()
