# the material(s) properties such as rho, Cp, k, h might change as functions of T, so we should not treat them as constants

# using simple linear interpolation to get the corresponding coefficients

# material 1: Aluminium
# density: ADD REFERENCE HERE! Could be a link or the name of the paper or commerical software...
def rho_Aluminium(T): # using simple linear interpolation to get the corresponding rho
  
  T_list = [300,400,500,600,700,800,900,1000,1100,1200] # 10 values here
  rho_list = [2700, 2675, 2650, 2625, 2600, 2575, 2550, 2400, 2350, 2300] # 10 values here as well
  rho = 0
  index = 0
  for i in range(len(T_list)):
    if T_list[i] <= T <= T_list[i+1]:
      index = i
      break
    elif T < T_list[0]:
      rho_min = 2700
      return rho_min
    elif T > T_list[-1]:
      rho_max = 
      return rho_max
  rho = (T - T_list[index])/100 * (rho_list[index+1] - rho_list[index])
  
  return rho

# specific heat capacity: TODO: reference!
def Cp_Aluminium(T): 
  T_list = [300,400,500,600,700,800,900,1000,1100,1200] # 10 values here
  rho_list = [...] # 10 values here as well
  
  index = 0
  for T in range(len(T_list)):
    if T_list[i] <= T <= T_list[i+1]:
      index = i
      break
    elif T
  rho = (T - T_list[index])/100 * (rho_list[index+1] - rho_list[index])
  
  return rho


def k_Aluminium(T):
  # TODO
  
  return

def _lambda_Aluminium(T):
  # TODO
  
  return


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
