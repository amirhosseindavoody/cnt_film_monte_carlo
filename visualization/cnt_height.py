import numpy as np
import matplotlib.pyplot as plt

filename = "/Users/amirhossein/research/test/ycoordinates.dat"
data = np.loadtxt(filename)

plt.figure()
plt.plot(data, linewidth=3)
plt.title("mean y coordinate of carbon nanotubes")
plt.show()
