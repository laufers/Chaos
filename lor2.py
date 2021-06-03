# lorenz.py: Lorenz attractor
# 
# The equations governing "Lorenz attractor" motion are quite simple.
#   dx/dt = sigma*(y - x)
#   dy/dt = x*(r - z) - y
#   dz/dt = x*y - b*z
#
# The three parameters sigma, r, and b are called respectively the 
# Prandtl number, the Rayleigh ratio, and the geometric factor. But
# more to the point they have typical 'pretty' values of 10, 28, and 
# 8/3. It is suggested that one experiment with the r-value first
# and that in particular r = 99.96 is interesting.
#
# This program establishes a canvas (and could set up more of them
# in separate windows) and then simply enters a time loop triggered by
# clicking on the Go button at the bottom of the frame. This time loop 
# calculates new (x, y, z) locations from the above equations and plots 
# a dot at (x, y). Clicking Go again will pick up where the plotting 
# left off. There is also code for a separate trajectory with coordinate
# variables (u, v, w).
# 
#
# Remarks on expanding this program
# =================================
# The equations above generate a vector field of velocities with the
# dual nature of producing an orderly pattern of motion which also has
# an element of apparent randomness. This leads to interesting ideas for
# altering this program to investigate further. Here are a few that
# occur to me off the top of my head.
#
# First there is code below that is commented out that color-codes the dots 
# based on the signs of dx, dy, and dz. This could also be used to color-code 
# for local field gradient or curl or whatever. Taking the idea of gradient,
# it is easy to get an idea of how the Lorenz figure tends to map local space 
# by simply watching the trajectory with time. There is also clearly an area
# of divergence at the intersection of the two 'ears'. So field divergence
# might be an interesting thing to color-code. Color coding only the trajectory
# will illuminate divergence on the Lorenz figure (aka the Lorenz set, I think.)
#
# It might be interesting to look at the initial conditions and the way that
# starting points off the Lorenz set tend to fall in. This could be visualized
# by choosing a cube of space not too far away and dividing it into a 3D grid
# of starting points. We would then map the trajectories for only a short time,
# the time necessary to fall onto the Lorenz set. I'd guess that many areas 
# will fall uniformly and some areas might exhibit interesting falling tracks.
#
# One of the key ideas of chaos and strange attractors is the high degree 
# of sensitivity to initial conditions. Instead of a single trajectory, this 
# program could easily accommodate multiple trajectories. For example we could
# track 2 tracks that proceed from two points initially spaced close together.
# This code is written into the program but commented out.
# 
# One could measure the physical separation of the two locations after some time 
# interval and use this to characterize the degree of chaotic behavior as a function
# of starting location. Plus it could be fun to watch the two trajectories diverge. 
# In fact (along the lines of a Julia set fractal) another program could color
# code the initial location based on this degree of divergence after a set time limit. 
# A meta-loop canvassing a set of initial locations would then provide a different 
# sort of 3D visualization.
# 
# Another way of measuring divergence (instead of divergence after a fixed time) is 
# to measure the time required to reach a set divergence distance of two trajectories.
# This escapes a certain random aspect of the former approach.
#
# This idea of two trajectories starting from nearby points could be extended to a 
# little 3 x 3 x 3 cube of points with an appropriate extension of the idea of a 
# divergence metric, i.e. one would look at how far apart or distorted the
# initial cube gets up after a certain time, or the other way, measure how long
# it takes to reach a set degree of distortion. 
#
# Another way of investigating the divergence behavior is to use a single starting
# location and compare different choices of time increment. By the same rationale
# any difference in time increment will change the length of the position displacement
# and thereby produce (after one time step) the same effect of two start locations. 
#
# To Do
#   Put in a Scram button
#   Make the single canvas into three: xy, xz, yz
#   Make the r parameter an optional command line argument
#   Project the figure onto an arbitrary plane
#

# for access to command line arguments
import sys

# for graphics; canvas etc
from   Tkinter  import *

def deathHandler(event):
  """death-handling function."""
  print ("getting nuked. so long...")

# This method is called when the Go button is clicked.
def plot_lorenz():
  """Plot an ascii x-value file versus an ascii y-value file."""
  print ("plot_lorenz() running...")
  global x0, y0, z0
  global u0, v0, w0
  global xscale, xshift, yscale, yshift
  global sigma, r, b
  global dt, T
  FILL   = 'blue'
  FILL2  = 'red'
  for x in range(int(T/dt) + 1): 
    t  = float(x) * dt
    dx = dt * (sigma*(y0 - x0))
    dy = dt * (x0*(r-z0)-y0)
    dz = dt * (x0*y0 - b*z0)
    x1 = x0 + dx
    y1 = y0 + dy
    z1 = z0 + dz

    # second trajectory code
    du = dt * (sigma*(v0 - u0))
    dv = dt * (u0*(r-w0)-v0)
    dw = dt * (u0*v0 - b*w0)
    u1 = u0 + du
    v1 = v0 + dv
    w1 = w0 + dw

    # This section of code sets the dot color
    # based on some criteria.
#     if dx > 0:
#       if dy > 0:
#         if dz > 0:
#           FILL = 'blue'
#         else:
#           FILL = 'cyan'
#         FILL = 'red'
#       else:
#         if dz > 0:
#           FILL = 'green'
#         else:
#           FILL = 'yellow'
#         FILL = 'maroon'
#     else:
#       if dy > 0:
#         if dz > 0:
#           FILL = 'red'
#         else:
#           FILL = 'purple'
#         FILL = 'blue'
#       else:
#         if dz > 0:
#           FILL = 'orange'
#         else:
#           FILL = 'violet'
#         FILL = 'purple'


    xp = int(x1*xscale + xshift)
    yp = int(y1*yscale + yshift)
    c.create_line(xp,yp,xp+1,yp,fill=FILL)
    x0 = x1
    y0 = y1
    z0 = z1

    # second trajectory code
    up = int(u1*xscale + xshift)
    vp = int(v1*yscale + yshift)
    c.create_line(up,vp,up+1,vp,fill=FILL2)
    u0 = u1
    v0 = v1
    w0 = w1

    c.update()

# Initialize 
sigma        = 10.
r            = 28.
b            = 8./3.
x0           = 22.
y0           = -10.
z0           = 14.5

# Code for a second trajectory 
delta_x0     = .0001
delta_y0     = .0001
delta_z0     = .0001
u0           = x0 + delta_x0
v0           = y0 + delta_y0
w0           = z0 + delta_z0

# dt is the time increment, T the time interval per Go-button click
dt           = 0.0025
T            = 1000.

# These scale/shift values center and scale the figure on the canvas.
xscale       = 10.
xshift       = 400.
yscale       = 10.
yshift       = 300.

root         = Tk()
frame        = Frame(root)

# Go button starts/continues the motion
go_button    = Button(frame, text='Go')
go_button.pack(side = RIGHT)
go_button.config(command=plot_lorenz)

frame.pack(side=BOTTOM)
frame        = Frame(root)
c            = Canvas(root)
c.pack(side=TOP, fill=BOTH, expand=1)
root.mainloop()