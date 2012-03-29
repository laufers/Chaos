import numpy as np
# alias for numpy
import matplotlib.pyplot as plt
# imports matlab like ploting and alaised

N = 10000
a = 10.
b = 69.
c = 8/3.
h = .001
x = np.zeros(N+1)
y = np.zeros(N+1)
z = np.zeros(N+1)

x1 = np.zeros(N+1)
y1 = np.zeros(N+1)
z1 = np.zeros(N+1)

x[0] = 10.
y[0] = 9.
z[0] = 9.

x1[0] = 10.00001 
y1[0] = 9.
z1[0] = 9.

for t in range(N):
    x[t + 1] = x[t] + h * a * (y[t] - x[t])
    y[t + 1] = y[t] + h * (x[t] * (b - z[t]) - y[t])
    z[t + 1] = z[t] + h *( (x[t] * y[t]) - (c * z[t]))
    x1[t + 1] = x1[t] + h * a * (y1[t] - x1[t])
    y1[t + 1] = y1[t] + h * (x1[t] * (b - z1[t]) - y1[t])
    z1[t + 1] = z1[t] + h *( (x1[t] * y1[t]) - (c * z1[t]))

# plt.plot(x,y)
# plt.figure()
# plt.plot(x)
# plt.plot(y)
# plt.plot(z)

# plt.plot(x,y)
plt.figure()
plt.subplot(3,1,1)
plt.plot(x, color = 'blue')
plt.plot(x1, color = 'red')

plt.subplot(3,1,2)
plt.plot(y, color = 'blue')
plt.plot(y1, color = 'red')

plt.subplot(3,1,3)
plt.plot(z, color = 'blue')
plt.plot(z1, color = 'red')

# 
#  
# fig = plt.figure()
# from mpl_toolkits.mplot3d import Axes3D
# ax = Axes3D(fig)
# ax.plot(x,y,z)

plt.show()