# the material(s) properties such as rho, Cp, k, h might change as functions of T, so we should not treat them as constants

# using simple linear interpolation to get the corresponding coefficients

# material 1: Aluminium
# density: ADD REFERENCE HERE! Could be a link or the name of the paper or commerical software...
def rho_Aluminium(T): # Kg/m3
  T_list = [300,400,500,600,700,800,900,1000,1100,1200] # 10 values here
  rho_list = [2700, 2675, 2650, 2625, 2600, 2575, 2550, 2400, 2350, 2300] # 10 values here as well
  rho = 0
  index = 0
  for i in range(len(T_list)):
    if T_list[i] <= T <= T_list[i+1]:
      index = i
      break
    elif T < T_list[0]:
      return rho_list[0]
    elif T > T_list[-1]:
      return rho_list[-1]
  rho = (T - T_list[index])/100 * (rho_list[index+1] - rho_list[index])
  
  return rho

# specific heat capacity: TODO: reference!
def Cp_Aluminium(T):  # J/Kg-K
  T_list = [300,400,500,600,700,800,900,1000,1100,1200] # 10 values here
  Cp_list = [900, 955.5, 994.8, 1034, 1088, 1136, 1190, 1238, 1281, 1321] # 10 values here as well
  Cp = 0
  index = 0
  for i in range(len(T_list)):
    if T_list[i] <= T <= T_list[i+1]:
      index = i
      break
    elif T < T_list[0]:
      return Cp_list[0]
    elif T > T_list[-1]:
      return Cp_list[-1]
  Cp = (T - T_list[index])/100 * (Cp_list[index+1] - Cp_list[index])

  return Cp


def k_Aluminium(T): # W/m-K
  T_list = [300,400,500,600,700,800,900,1000,1100,1200] # 10 values here
  k_list = [150, 160, 165, 170, 170, 170, 170, 170, 170, 170] # 10 values here as well
  k = 0
  index = 0
  for i in range(len(T_list)):
    if T_list[i] <= T <= T_list[i+1]:
      index = i
      break
    elif T < T_list[0]:
      return k_list[0]
    elif T > T_list[-1]:
      return k_list[-1]
  k = (T - T_list[index])/100 * (k_list[index+1] - k_list[index])

  return k

def _lambda_Aluminium(T):

  return k_Aluminium(T)/(Cp_Aluminium(T)*rho_Aluminium(T))


# heat transfer coefficient
def h_Aluminium(T):
  # TODO
  
  return

# TODO
# material 2: Inconel 800HT

def rho_Inconel800HT(T):
  
  return

def Cp_Inconel800HT(T):
  
  return


def k_Inconel800HT(T):
  
  
  return

def h_Inconel800HT(T):
  
  
  return



# TODO
# material 3: Ethane Coke

def rho_EthaneCoke(T):
  
  return

def Cp_EthaneCoke(T):
  
  return


def k_EthaneCoke(T):
  
  
  return

def h_EthaneCoke(T):
  
  
  return



# TODO
# material 4: Naphtha Coke

def rho_NaphthaCoke(T):
  
  return

def Cp_NaphthaCoke(T):
  
  return


def k_NaphthaCoke(T):
  
  
  return

def h_NaphthaCoke(T):
  
  
  return
