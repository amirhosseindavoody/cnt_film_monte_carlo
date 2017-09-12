# import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt


filename = "/Users/amirhossein/research/test/population_profile.dat"
data = np.loadtxt(filename)

plt.ion()
fig = plt.figure()
ax = fig.gca()
for i in range(0,data.shape[0]):
	ax.plot(data[i,1:],linewidth=2)
ax.plot(data[-1,1:],linewidth=4, color="black")

fig = plt.figure()
ax = fig.gca()
for i in range(1,data.shape[1]):
	ax.plot(data[:,0], data[:,i],linewidth=2)

input("Press Enter to exit...")
