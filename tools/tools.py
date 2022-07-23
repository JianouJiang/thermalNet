import os
# tool box:

# read parameters
default_parameters_directory = "../heat_transfer/parameters.txt"
def readParameters(parameters_directory):
  if parameters_directory==default_parameters_directory:
      print("reading default parameters from: "+parameters_directory)

  parameter_list = []  
  with open(parameters_directory) as fp:
    Lines = fp.readlines()

    for i in range(len(Lines)):
      line = Lines[i]
      if i==0: # skipping the first line: L, dx, dt...
        continue
      else:
        line_strip = line.replace('\n', "")
        line_strip = line.replace(' ', "")
      
        parameter_list.append(float(line_strip)) 
  L, dx, t_max, dt, _lambda1, _lambda2, number_of_ghost_points, num_of_timeSteps_for_plotting = parameter_list[0], parameter_list[1], parameter_list[2], parameter_list[3], parameter_list[4], parameter_list[5], parameter_list[6], parameter_list[7]
  print("L="+str(L)+" dx="+str(dx)+" t_max="+str(t_max)+" dt="+str(dt) + " _lambda1="+str(_lambda1)+" _lambda2="+str(_lambda2)+" number_of_ghost_points="+str(number_of_ghost_points))
  return L, dx, t_max, dt, _lambda1, _lambda2, number_of_ghost_points, num_of_timeSteps_for_plotting




# integrating from start to end
L=1
def integral(x, f, x_start=0, x_end=L): # by default
    dx = x[1] - x[0]
    sum_fx = 0
    for i in range(0,(len(f)-1)):
        xi = x[i]
        if x_start<=xi<=x_end:
            fi = (f[i]+f[i+1] )/2 # central differencing
            sum_fx = sum_fx + fi*dx
    return sum_fx
    
# open file in write mode and write data
def writeData(directory, ti, x, T, _lambda, mask=None): # note T and lambda here are lists, whereas time t is single value
    isFile = os.path.isfile(directory) 
    with open(directory, 'a') as fp:
        index = 0
        for Ti in T:
            xi = x[index]
            _lambda_i = _lambda[3][index] # only taking out lambda, but in the future, we can take out rho, cp, and k
            # write each ti, xi, Ti, lambda on a new line
            if mask==None:
                line = str(ti) + " " + str(xi) + " " + str(Ti) + " " + str(_lambda_i) + "\n"
                fp.write(line)
            else:
                mask_i = mask[index]
                line = str(ti) + " " + str(xi) + " " + str(Ti) + " " + str(_lambda_i) + " "+str(mask_i)+ "\n"
                fp.write(line)

            index = index + 1
    return
    
# open file in write mode and write data
def writeData2D(directory, ti, x, T, _lambda,mask=None): # note T and lambda here are lists, whereas time t is single value
    isFile = os.path.isfile(directory) 
    with open(directory, 'a') as fp:
        for i in range(len(T)):
            for j in range(len(T[0])):
                xij = x[i][j][0]
                yij = x[i][j][1]
                Tij = T[i][j]
                _lambda_ij = _lambda[3][i][j] # only taking out lambda, but in the future, we can take out rho, cp, and k
                if mask==None:
                    # write each ti, xi, Ti, lambda on a new line
                    line = str(ti) + " " + str(xij) + " " + str(yij) + " " + str(Tij) + " " + str(_lambda_ij) + "\n"
                    fp.write(line)
                else:
                    mask_ij = mask[i][j]
                  # write each ti, xi, Ti, lambda on a new line
                    line = str(ti) + " " + str(xij) + " " + str(yij) + " "  + str(Tij) + " " + str(_lambda_ij) + " " + str(mask_ij)+"\n"
                    fp.write(line)
    return
     
# reading .csv data from ansys and writing to .txt for python
import csv
def convertCSVtoTxT(directory_read, directory_write):
  i = 0
  with open(directory_read, 'r') as csvfile:
    with open(directory_write, 'a') as txtfile:
      data = csv.reader(csvfile, skipinitialspace=True)
      for row in data:
        if i==0: # skipping the first line: S (m) X Coordinate (m)  Y Coordinate (m)  Z Coordinate (m)  Temperature (°C)  Time (s)
          i = i + 1
          continue
        else:
          ti = row[5]
          xi = row[0]
          Ti = float(row[4])
          lambda_i = 1 # lambda_Aluminium(Ti)
          line = ti + " " + xi + " " + str(Ti) + " " + str(lambda_i)  + "\n"
          txtfile.write(line)
        i = i + 1  
  return

# deleteFile so that the data do not append to the same .txt over and over again
def deleteFile(fileName):
  if os.path.exists(fileName):
    os.remove(fileName)

  return