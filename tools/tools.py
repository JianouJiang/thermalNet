import os
# tool box:

# read parameters
def readParameters():
  parameters_directory = "parameters.txt"
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
  L, dx, t_max, dt, _lambda1, _lambda2, number_of_ghost_points = parameter_list[0], parameter_list[1], parameter_list[2], parameter_list[3], parameter_list[4], parameter_list[5], parameter_list[6]
  print("L="+str(L)+" dx="+str(dx)+" t_max="+str(t_max)+" dt="+str(dt) + " _lambda1="+str(_lambda1)+" _lambda2="+str(_lambda2)+" number_of_ghost_points="+str(number_of_ghost_points))
  return L, dx, t_max, dt, _lambda1, _lambda2, number_of_ghost_points


L, dx, t_max, dt, _lambda1, _lambda2, number_of_ghost_points = readParameters()


# integrating from start to end
def integral(f): 
    
    sum_fx = 0
    for i in range(0,(len(f)-1)):
        fi = (f[i]+f[i+1] )/2 # central differencing
        sum_fx = sum_fx + fi*dx
    return sum_fx
    
# open file in write mode and write data
def writeData(directory, ti, T, _lambda): # note T and lambda here are lists, whereas time t is single value
    isFile = os.path.isfile(directory) 
    with open(directory, 'a') as fp:
        index = 0
        for Ti in T:
            xi = -dx + dx * index
            _lambda_i = _lambda[index]
            # write each ti, xi, Ti, lambda on a new line
            line = str(ti) + " " + str(xi) + " " + str(Ti) + " " + str(_lambda_i) + "\n"
            fp.write(line)
            index = index + 1
    return
    
 
