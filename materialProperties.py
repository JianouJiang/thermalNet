# the material(s) properties such as rho, Cp, k, h might change as functions of T, so we should not treat them as constants

# using simple linear interpolation to get the corresponding coefficients

# material 1: Aluminium
# density
def rho_Aluminium(T): # using simple linear interpolation to get the corresponding rho
  
  T_list = [300,400,500,600,700,800,900,1000,1100,1200] # 10 values here
  rho_list = [...] # 10 values here as well
  
  index = 0
  for t in range(len(T_list)):
    if T_list[i] <= T <= T_list[i+1]:
      index = i
      break
  
  rho = (T - T_list[index])/100 * (rho_list[i+1] - rho_list[i])
  
  return rho

# specific heat capacity
def Cp_Aluminium(T): 
  # TODO
  
  return


def k_Aluminium(T):
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
