# reading .csv data from ansys and writing to .txt for python
print("convertCSVtoTxT()")
import time
import sys
sys.path.insert(0, '../')
import numpy as np
from tools.tools import *


directory_read="Al_uniform_IC_BC_left0_right1000/Al_constant_IC_left0_right1000.csv"
directory_write = "../data/ansys_linear0_Dirichlet_Aluminium_T.txt"


convertCSVtoTxT(directory_read, directory_write)