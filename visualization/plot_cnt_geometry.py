import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import sys
import os # this is for opening and closing directories.

plt.ion() # enables interactive mode for ploting. No plt.show() is needed!

fig = plt.figure()
ax = fig.gca(projection='3d')

directory = "/home/amirhossein/google_drive/research/exciton/data/cnt_mesh/output_results/"

for file in os.listdir(directory):
	if(file.startswith("cnt")):

		cnt_curve = np.loadtxt(directory+file, skiprows=0)

		x = cnt_curve[:,0]
		y = cnt_curve[:,1]
		z = cnt_curve[:,2]

		ax.plot(x, y, z, linewidth=2, color="blue")

ax.set_aspect('equal', 'datalim')

input("Press Enter to exit...")



# import matplotlib as mpl
# from mpl_toolkits.mplot3d import Axes3D
# import numpy as np
# import matplotlib.pyplot as plt
# import sys
# import itertools

# plt.ion() # enables interactive mode for ploting. No plt.show() is needed!

# fig = plt.figure()
# ax = fig.gca(projection='3d')

# directory = "/home/amirhossein/google_drive/research/exciton/data/cnt_mesh/output_results/"
# cnt_curve = np.loadtxt(directory+"bullet_physics_points.dat", skiprows=0)

# x = cnt_curve[:,0]
# y = cnt_curve[:,1]
# z = cnt_curve[:,2]

# ax.plot(x, y, z, linewidth=3)

# cnt_curve = np.loadtxt(directory+"pca_points.dat", skiprows=0)

# x = cnt_curve[:,0]
# y = cnt_curve[:,1]
# z = cnt_curve[:,2]

# for i in range(0,x.size,2):
# 	ax.plot(x[i:i+2], y[i:i+2], z[i:i+2], color="red", linewidth=3)


# input("Press Enter to exit...")