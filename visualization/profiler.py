import numpy as np
import matplotlib.pyplot as plt

# filename = "/Users/amirhossein/Desktop/test_new/population_profile.dat"
filename = "/Users/amirhossein/Desktop/test_new/region_current.dat"
data = np.loadtxt(filename)
time = data[:,0]
data = data[:,1:]

sum_data = np.sum(data,1)

# plt.style.use('classic')

# # plot bare data
# for i in range(0,data.shape[1]):
# 	plt.plot(time,data[:,i], linewidth=2)

# print(plt.style.available)


# calculate and plot box average of the data
box = np.ones(500)
for i in range(0,data.shape[1]):
	avg = np.convolve(data[:,i], box, mode='valid')/np.sum(box)
	plt.plot(time[:len(avg)],avg, linewidth=3)


# avg = np.convolve(sum_data, box, mode='valid')/np.sum(box)
# plt.plot(time[:len(avg)],avg, linewidth=3)
plt.plot(time,sum_data, linewidth=3, color='black')

plt.show()
