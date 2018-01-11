import numpy as np
# alias for numpy
import matplotlib.pyplot as plt
# imports matlab like ploting and alaised

N = 100000
dat = np.empty((6,4))

# Matrix
#a
dat[0, 0] = 0
dat[0, 1] = .85
dat[0, 2] = .2
dat[0, 3] = -.15
#b
dat[1, 0] = 0
dat[1, 1] = .04
dat[1, 2] = -.26
dat[1, 3] = .28
#c
dat[2, 0] = 0
dat[2, 1] = -.04
dat[2, 2] = .23
dat[2, 3] = .26
#d
dat[3, 0] = .16
dat[3, 1] = .85
dat[3, 2] = .22
dat[3, 3] = .24
#e
dat[4, 0] = 0
dat[4, 1] = 0
dat[4, 2] = 0
dat[4, 3] = 0
#f
dat[5, 0] = 0
dat[5, 1] = 1.6
dat[5, 2] = 1.6
dat[5, 3] = .44

x = np.zeros(N + 1)
y = np.zeros(N + 1)

for t in range(N):
    r = np.random.rand(1)
  
    if 0 < r <= 0.01 : use = 0
    elif 0.01 < r <= .85 : use = 1
    elif 0.85 < r <= 0.94 : use = 2
    elif 0.94 < r < 1 : use = 3
  
    x[t + 1] = dat[0, use] * x[t] + dat[1, use] * y[t] + dat[4, use]
    y[t + 1] = dat[2, use] * x[t] + dat[3, use] * y[t] + dat[5, use]

plt.scatter(x,y,c = "r", edgecolor = "b", marker = ",", s = 1)
plt.show()
