# read parameters
def readParameters(parameters_directory):
  parameter_list = []  
  with open(parameters_directory) as fp:
    Lines = fp.readlines()
    for line in Lines:
      parameter_list.append(line[2])     
  L, dx, t_max, dt, _lambda1, _lambda2 = parameter_list[0], parameter_list[1], parameter_list[2], parameter_list[3], parameter_list[4], parameter_list[5]
  print("L="+str(L)+" dx="+str(dx)+" t_max="+str(t_max)+" dt="+str(dt) + " _lambda1="+str(_lambda1)+" _lambda2="+str(_lambda2))
  return L, dx, t_max, dt, _lambda1, _lambda2
