# reading .csv data from ansys and writing to .txt for python
print("convertCSVtoTxT()")
import time
import sys
sys.path.insert(0, '../')
import numpy as np
from tools.tools import *


directory_read="Al_Uniform_IC_BC_left500celsius_right_0_celsius/Al_Uniform_IC_BC_left500celsius_right_0_celsius.csv"
directory_write = "../data/ansys_linear0_Dirichlet500_0_Aluminium_T.txt"


convertCSVtoTxT(directory_read, directory_write)