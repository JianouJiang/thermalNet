# importing libs
#from jianouTest import * # importing all functions from jianouTest.py
#from kaiyiTest import *
from tools.normalisation import *
import matplotlib.pyplot as plt

# TODO!
# from chennuoTest import * 
def kaiyiTest():
    print("Howdy friends, my name is Kaiyi Dai, a 3rd Earth & Planetary science undergraduate student in IC.")
    print("I'm interested in Mars rover missions and don't hesitate to speak to me if you enjoy planetary stuff too!")
    return

def xiaoyangTest():
    print("Hello! I'm Sean, year 2 student studying EIE!")
    print("I love playing guitar and are thrilled to work with you guys : )")
    return

def yanlinTest():
    print("Hey guys, I'm Yanlin or you may call me Philip, just graduated from Imperial College last year with Msci Physics.")
    print("Currently I'm a system testing engineer at Huawei.")
    return

def jianouTest():
  print("My name is Jianou, good to have everybody here! Testing plotTest() function below:")
  plotTest()
  return

def plotTest():
  print("--------plotTest()---------")
  pred_Benchmark = [1,2,3,4,5,6,7,8,9]
  predList_test_predictions = [9,8,7,6,5,4,3,2,1]
  plt.plot(pred_Benchmark, label='Benchmark')
  plt.plot(predList_test_predictions, label='thermalNet 3.0')
  plt.xlabel('X-axis') 
  plt.legend(loc='upper right')
  plt.title('Normalised Thickness H Comparison')
  #fig.savefig('Normalised_Thickness_H_Comparison.pdf', bbox_inches='tight', dpi=150)
  plt.show()
  return

print("Started testing all contributors...")
def test():
  
  print("--------test()---------")
  print("just testing from chennuo")
  
  print("First, let's welcome Jianou...")
  jianouTest()
  normalisation()
  
  print("Then, let's welcome Chen Nuo...")
  # ChenNuo's function below:
  # chennuoTest()
  
  print("Next, let's welcome Xiao Yang...")
  # Xiao Yang's function below:
  # xiaoyangTest()
  
  print("Now, let's welcome Kai Yi...")
  # Kai Yi's function below:
  # kaiyiTest()
  
  print("Again, let's welcome Zheng Yang...")
  # Zheng Yang's function below:
  # zhengyangTest()
  
  print("Finally, let's welcome Yan Lin...")
  # Yan Lin's function below:
  # yanlinTest()
     
  return

# testing if all contributors can push their codes.
test()
print("Finished testing all contributors :)")
