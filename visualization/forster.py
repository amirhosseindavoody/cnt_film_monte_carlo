import numpy as np
import matplotlib.pyplot as plt


a = 0.1
n = int(1e7)

t = np.random.random(n)
r = a/np.cbrt(1-t)


x = np.linspace(a,10*a, num=100)
y = 3*a**3/x**4

hist, edge = np.histogram(r,x)
hist = hist/hist[0]*y[0]

plt.figure()
plt.plot(x,y, linewidth=3, color='b')
plt.draw()

# plt.figure()
# plt.hist(r, bins=100, edgecolor="w")
plt.plot(edge[:-1],hist, linewidth=3, color='r')
plt.show()
