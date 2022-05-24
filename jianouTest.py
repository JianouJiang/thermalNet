import matplotlib.pyplot as plt


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
