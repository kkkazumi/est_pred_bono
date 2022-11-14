import matplotlib.pyplot as plt
import numpy as np

print("input filename")
filename=input()
data=np.loadtxt(filename,delimiter=",")
print(data.shape)

plt.plot(data[:,0],label="min")
plt.plot(data[:,1],label="mid")
plt.plot(data[:,2],label="max")
plt.legend()
plt.show()

