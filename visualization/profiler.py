import numpy as np
import matplotlib.pyplot as plt

filename = "/Users/amirhossein/Desktop/test/population_profile.dat"
data = np.loadtxt(filename)
time = data[:,0]
data = data[:,1:]

# plt.style.use('classic')

# plot bare data
for i in range(0,data.shape[1]):
	plt.plot(time,data[:,i], linewidth=2)

# print(plt.style.available)


# calculate and plot box average of the data
box = np.ones(300)
for i in range(0,data.shape[1]):
	avg = np.convolve(data[:,i], box, mode='valid')/np.sum(box)
	plt.plot(time[:len(avg)],avg, linewidth=3)

plt.show()
