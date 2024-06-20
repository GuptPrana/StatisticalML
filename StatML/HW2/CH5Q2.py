import matplotlib.pyplot as plt
import numpy as np

n = np.arange(100000)+1
p = np.power(1-1/n, n)
p = 1 - p

plt.plot(n, p)
plt.show()