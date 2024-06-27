import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("1.txt")
### print (data.shape)

t,x,v = np.loadtxt("1.txt",unpack=True)   ### unpack=True は各列を個別に読み込む

plt.plot(t,x)