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
plotColor = [None] * (N + 1)

x[0] = 10.
y[0] = 9.
z[0] = 9.


for t in range(N):
    x[t + 1] = x[t] + h * a * (y[t] - x[t])
    y[t + 1] = y[t] + h * (x[t] * (b - z[t]) - y[t])
    z[t + 1] = z[t] + h *( (x[t] * y[t]) - (c * z[t]))
    if t > 5000 :
        plotColor[t] = 'red'
    else:
        plotColor[t] = 'blue'
    

plt.figure()
plt.plot(x,y)
plt.xlabel("Magnitude x")
plt.ylabel("Magnitude y")

plt.figure()
plt.plot(x,z)
plt.xlabel("Magnitude x")
plt.ylabel("Magnitude z")

plt.figure()
plt.plot(y,z)
plt.xlabel("Magnitude y")
plt.ylabel("Magnitude z")

plt.figure()
plt.plot(x, label = 'x')
plt.plot(y, label = 'y')
plt.plot(z, label = 'z')
plt.xlabel("Timestep")
plt.ylabel("Magnitude")
# plt.legend(loc = 'lower right')
plt.legend(bbox_to_anchor=(1.06,1), loc = 'upper center', ncol = 1)
# plt.show()

fig = plt.figure()
from mpl_toolkits.mplot3d import Axes3D
ax = Axes3D(fig, auto_add_to_figure=False)
fig.add_axes(ax)
ax.plot(x,y,z,alpha = float(t)/(N-1),color='green')

plt.show()