# lorenz.py 
# Python package for the computational lab on the Lorenz Map used 
# in PHY 460/1460S (Non-linear Dynamics) offered by the University of Toronto Dept. of Physics. 
# Written by Peter Hitchcock, February 2010.
# Based on a Mathematica package (lorenz.m) written by David Harrison and C. G. Wong

import pylab as py
import numpy as np
from scipy import integrate

def lorenz(xs, t, r, p=10, b=8./3.):
# {{{
  ''' lorenz(xs, t, r, p, b)
    Returns RHS of the Lorenz equations for the given parameters. Used by trajectory().'''
  x, y, z = xs
  return p*(y - x), r*x - y - x*z, x*y - b*z
# }}}

def trajectory(r=22, p=10, b=8/3., x0=(0, 1., 0.2), tf=90, steps=1e4):
# {{{
  ''' trajectory([r, p, b, x0, tf, steps])
      Computes a solution to the Lorenz equations with parameters r, p, and b,
      and initial conditions x0, using scipy.integrate.odeint(). The trajectory
      is computed for 'steps' equally spaced time intervals until time tf. Note
      that this only determines the grid on which the solution is returned; the
      time step is determined automatically by the ODE solver. Returns (ts,
      sln), where sln is a tuple of the X, Y, and Z time series.'''
  ts = np.linspace(0, tf, steps)
  sln = integrate.odeint(lorenz, x0, ts, (r, p, b))

  return ts, (sln[:, 0], sln[:, 1], sln[:, 2])
# }}}

def make_title(r, **kwargs):
# {{{
  '''make_title(r, [p, b, x0])
  Returns a description of parameters used to compute the trajectory.'''
  x0 = kwargs.pop('x0', (0, 1., 0.2))
  p = kwargs.pop('p', None)
  b = kwargs.pop('b', None)

  title = r'$r = %.2f,\ ' % r
  if p is not None: title += 'p = %.1f,\ ' % p
  if b is not None: title += 'b = %.1f,\ ' % b
  title += 'X_0 = (%.1f,\ %.1f,\ %.1f)$' % (x0[0], x0[1], x0[2])

  return title
# }}}

def find_maxima(ts, sln, dim = 'Z'):
# {{{
  ''' find_maxima(ts, sln, [dim])
      Returns the time and location of the local maxima along the specified
      dimension ('X', 'Y' or 'Z'). Uses quadratic interpolation to better
      approximate the location of the maxima.'''
  # Extract specified dimension
  idim = {'X':0, 'Y':1, 'Z':2}
  xs = sln[idim[dim]]

  # Compute differences; find local maxima (clip last few points to ensure we can interpolate)
  xd = xs[1:] - xs[:-1]
  ixs = np.where((xd[:-1] > 0.) & (xd[1:] < 0.))[0]
  ixs = ixs[ixs < len(xs) - 2]

  # Use the Lagrangian polynomial interpolation formula to better approximate maxima
  l0 = np.poly1d([ 0.5, -1.5,  1. ])
  l1 = np.poly1d([-1,  2,  0])
  l2 = np.poly1d([0.5, -0.5, 0])
  L = lambda x1, x2, x3: float(x1) * l0 + float(x2) * l1 + float(x3) * l2

  # Find the roots of the derivative of the interpolating polynomial
  maxima = [np.real(np.roots(L(xs[i], xs[i+1], xs[i+2]).deriv())[0]) for i in ixs]

  # Interpolate t, x, y, and z
  dt = ts[1] - ts[0]
  mts = np.array([ts[i] + m*dt for i, m in zip(ixs, maxima)])

  x, y, z = sln
  mxs = np.array([L(x[i], x[i+1], x[i+2])(m) for i, m in zip(ixs, maxima)])
  mys = np.array([L(y[i], y[i+1], y[i+2])(m) for i, m in zip(ixs, maxima)])
  mzs = np.array([L(z[i], z[i+1], z[i+2])(m) for i, m in zip(ixs, maxima)])

  return mts, (mxs, mys, mzs)
# }}} 

def power_spectrum(ts, sln, dim = 'X', window=None):
 # {{{
  ''' power_spectrum(ts, sln, [dim, window])
      Returns the frequencies and power of the variable specified by dim from
      the given trajectory sln, windowed by a Hanning (sin) function. A specific
      segment of the trajectory can by specified by window.'''
  from numpy.fft import rfft

  # Pick the requested time series
  idim = {'X':0, 'Y':1, 'Z':2}
  xs = sln[idim[dim]]

  # Find indices of specified window (if any)
  if window is None:
    iwin = [0, len(xs) - 1]
  else:
    iwin = [np.argmin((ts - window[0])**2), np.argmin((ts - window[1])**2)]

  if iwin[1] - iwin[0] < 2:
    raise "Time window must include more than two points"

  # Window & remove the mean of the time series
  xwin = xs[iwin[0]:iwin[1] + 1]
  xwin = (xwin - np.mean(xwin)) * np.hanning(iwin[1] + 1 - iwin[0])

  # Compute fourier transform, power spectral density
  xks = rfft(xwin)
  ws = 2*np.pi*(np.arange(len(xks))) / (ts[iwin[1]] - ts[iwin[0]])
  psd = np.real(xks * np.conj(xks))

  return ws, psd
# }}}

def plot_timeseries(r, dim='Z', showmax=None, **kwargs):
# {{{
  ''' plot_timeseries(r, [showmax, dim, p, b, x0, tf, steps])
      Computes and plots time series of one of the three variables (dim = 'X',
      'Y', or 'Z').  r must be specified, initial conditions, p and b can
      optionally be specified, otherwise default values of x0 = (0, 1., 0.2),
      p=10 and b=8/3 are used. tf and steps specify how long to carry out the
      integration, and at what temporal resolution (see trajectory()).

      Optionally, the values at local maxima of any of three variables can by
      highlighted by specifying showmax (also 'X', 'Y', or 'Z'). By default no
      highlighing is done.  

      Examples:
      # Plots time series of X at r=10
      lorenz.plot_timeseries(10, 'X') 
      
      # Plots time series of Z at r=20, highlighting points where Z is at a
      # local maxima
      lorenz.plot_timeseries(20, 'Z', 'Z')  
                                            
      # Plots time series of Y at r=28, highlighting local points where Z is at
      # a local maxima
      lorenz.plot_timeseries(28, 'Y', 'Z')  
                                            
      # Plots time series of X at r=28, b=9/3., and initial conditions X, Y, Z =
      # (1., 2., 5.)
      lorenz.plot_timeseries(28, 'X', x0=(1., 2., 5.), b=9/3.).'''   
                                            
                                            
  # Compute solution
  ts, sln = trajectory(r, **kwargs)

  # Extract specified dimension
  idim = {'X':0, 'Y':1, 'Z':2}
  xs = sln[idim[dim]]

  py.figure(1)
  py.clf()
  py.plot(ts, xs, 'k', lw=1.0)
  py.xlabel('t')
  py.ylabel(dim)

  # If maxima are provided, plot them
  if showmax in idim.keys():
    rts, rsln = find_maxima(ts, sln, showmax)
    mxs = rsln[idim[dim]]
    py.plot(rts, mxs, 'bo', hold=True)

  py.title(make_title(r, **kwargs))

  py.draw()
# }}}

def plot_projection(r, absc='X', ord='Y', showmax=None, **kwargs):
# {{{
  ''' plot_projection(r, [absc, ord, showmax, p, b, x0, tf, steps])
      Computes and plots the projection of a trajectory in the plane specified
      by absc and ord.  r must be specified, initial conditions, p and b can
      optionally be specified, otherwise default values of x0 = (0, 1., 0.2),
      p=10 and b=8/3 are used. tf and steps specify how long to carry out the
      integration, and at what temporal resolution (see trajectory()).

      Optionally, the values at local maxima of any of three variables can by
      highlighted by specifying showmax (also 'X', 'Y', or 'Z'). By default no
      highlighing is done.  

      Examples:
      # Plots projection in XY plane for r = 10
      lorenz.plot_projection(10, 'X', 'Y') 
      
      # Plots projection in XZ plane, highlighting loxal maxima in 'Z', out to a
      # final time of 100
      lorenz.plot_projection(10, 'X', 'Z', tf=100.)'''

  # Compute solution
  ts, sln = trajectory(r, **kwargs)

  # Extract specified dimensions
  idim = {'X':0, 'Y':1, 'Z':2}
  xs = sln[idim[absc]]
  ys = sln[idim[ord]]

  py.figure(1)
  py.clf()
  py.plot(xs, ys, 'k', lw=1.0)
  py.xlabel(absc)
  py.ylabel(ord)

  # If maxima are requested, plot them
  if showmax in idim.keys():
    rts, rsln = find_maxima(ts, sln, showmax)
    mxs = rsln[idim[absc]]
    mys = rsln[idim[ord]]
    py.plot(mxs, mys, 'bo', hold=True)

  py.title(make_title(r, **kwargs))

  py.draw()
# }}}

def plot_3d(r, showmax=None, **kwargs):
# {{{
  ''' plot_3d(r, [showmax, p, b, x0, tf, steps])
      Computes and plots the trajectory in three dimensions. Requires the
      mplot3d toolkit to be available, which should be the case for versions of
      matplotlib after 0.99. r must be specified, initial conditions, p and b
      can optionally be specified, otherwise default values of x0 = (0, 1.,
      0.2), p=10 and b=8/3 are used. tf and steps specify how long to carry out
      the integration, and at what temporal resolution (see trajectory()).

      Optionally, the values at local maxima of any of three variables can by
      highlighted by specifying showmax (also 'X', 'Y', or 'Z'). By default no
      highlighing is done.  

      Examples:
      # Plots trajectory for r=28, computed on 1e4 time steps from t = 0 to 100
      lorenz.plot_3d(28, tf=100, steps=1e4) '''
  from mpl_toolkits.mplot3d import Axes3D
  idim = {'X':0, 'Y':1, 'Z':2}

  # Compute solution
  ts, sln = trajectory(r, **kwargs)

  # Extract specified dimensions
  xs, ys, zs = sln

  f = py.figure(1)
  f.clf()
  ax = Axes3D(f)
  ax.plot(xs, ys, zs, 'k', lw=1.0)
  ax.set_xlabel('X')
  ax.set_ylabel('Y')
  ax.set_zlabel('Z')
  f.text(0.05, 0.95, make_title(r, **kwargs), va='top', ha='left')

  # If maxima are requested, plot them
  if showmax in idim.keys():
    rts, rsln = find_maxima(ts, sln, showmax)
    mxs, mys, mzs = rsln
    ax.plot(mxs, mys, mzs, 'bo')

  py.draw()
# }}}

def plot_return_map(r, dim='Z', **kwargs):
# {{{
  ''' plot_return_map(r, [dim, p, b, x0, tf, steps])
      Computes trajectory and plots 1st return map defined by successive local
      maxima in the variable specified by dim (= 'X', 'Y', or 'Z'). r must be
      specified, initial conditions, p and b can optionally be specified,
      otherwise default values of x0 = (0, 1., 0.2), p=10 and b=8/3 are used.
      tf and steps specify how long to carry out the integration, and at what
      temporal resolution (see trajectory()). Also plots a dashed line with
      slope = 1.

      Examples:
      # Plots Z return map for r = 40
      lorenz.plot_return_map(40, 'Z') '''   

  # Compute solution
  ts, sln = trajectory(r, **kwargs)

  # Compute maxima
  rts, msln = find_maxima(ts, sln, dim)

  # Extract specified dimension
  idim = {'X':0, 'Y':1, 'Z':2}
  mzs = msln[idim[dim]]

  py.figure(2)
  py.clf()
  py.plot(mzs[:-1], mzs[1:], 'ko', ms=5.0)

  # Plot line with slope 1.0
  inear = np.where((mzs[1:]/mzs[:-1] > 0.95) & (mzs[1:]/mzs[:-1] < 1.05))[0]
  zm, zp = np.min(mzs[inear]), np.max(mzs[inear])
  py.plot([zm, zp], [zm, zp], 'k--', lw=2.0)
  py.title(dim + ' return map; ' + make_title(r, **kwargs))

  py.xlabel(r'$%s_k$' % dim, fontsize=16)
  py.ylabel(r'$%s_{k+1}$' % dim, fontsize=16)
  py.draw()
# }}}

def plot_lyapunov(r, x0=(0, 1, 0.2), dx=(1e-8, 0, 0), **kwargs):
# {{{
  ''' plot_lyapunov(r, [x0, dx, p, b, x0, tf, steps])
      Computes two trajectories with initial conditions x0 and x0 + dx, and
      plots the (natural) logarithm of the divergence between the two as a
      function of time. r must be specified; initial conditions, p and b can
      optionally be specified, otherwise default values of x0 = (0, 1., 0.2),
      p=10 and b=8/3 are used.  tf and steps specify how long to carry out the
      integration, and at what temporal resolution (see trajectory()). 

      Examples:
      # Plots divergence of two trajectories initially separated by (0, 1e-10,
      # 0) with r = 28
      lorenz.plot_lyapunov(28, dx=(0, 1e-10, 0) '''   

  # Compute first solution
  ts, sln = trajectory(r, x0 = x0, **kwargs)

  # Compute second solution
  tsp, slnp = trajectory(r, x0 = py.array(x0) + py.array(dx), **kwargs)

  # Compute lyapunov divergence
  dx0 = np.sqrt(np.sum(np.array(dx)**2))
  dxs = np.sqrt(np.sum((np.array(slnp) - np.array(sln))**2, axis=0))
  div = np.log(dxs/dx0)

  py.ioff()

  py.figure(4)
  py.clf()

  py.title(make_title(r, x0=x0, **kwargs) + '$,\ \delta X = (%.2g,\ %.2g,\ %.2g)$' % (dx[0], dx[1], dx[2]))
  py.plot(ts, div, 'k', lw=1.0)
  py.xlabel(r'$t$')
  py.ylabel(r'$\ln(|\frac{\delta X(t)}{\delta X(0)}|)$', fontsize='x-large')

  py.ion()
  py.show()
  py.draw()
# }}}

def plot_power_spectrum(r, dim='X', freqs=[0, 50], logy=True, twin=None, tf=150., **kwargs):
# {{{
  ''' plot_power_spectrum(r, [dim, freqs, logy, twin, x0, p, b, x0, tf, steps])
      Computes solution and power spectrum, and plots the spectral density for
      the given variable dim (= 'X', 'Y', or 'Z'). r must be specified; initial
      conditions, p and b can optionally be specified, otherwise default values
      of x0 = (0, 1., 0.2), p=10 and b=8/3 are used.  tf and steps specify how
      long to carry out the integration, and at what temporal resolution (see
      trajectory()). 

      freqs specifies the interval of frequencies to plot. By default the base-2
      logarithm of the spectral density is plotted; the raw spectral density can
      optionally be plotted by setting logy=False.

      You can specify a particular time period of the trajectory to compute the
      spectrum of using twin. By default Note that the spectral resolution is
      determined by the length of the trajectory computed (tf).  Unlike the
      other methods, tf is set by default to 150.

      By default this plot is *not* cleared; calling plot_power_spectrum
      multiple times will overplot successive spectra. A somewhat useful legend
      can be produced by calling legend() after the fact. Call clf() to clear
      the figure.

      Examples:
      # Plots power spectrum of X for an r=150 trajectory
      lorenz.plot_power_series(150, dim='X')

      # Plots power spectrum of an r=28 trajectory from time 10 to time 80
      lorenz.plot_power_series(150, twin=[10, 80])

      # Plots raw power spectrum of Z for an r=140 trajectory, computed out to
      # t=500
      lorenz.plot_power_series(140, dim='Z', logy=False, tf=500) '''   

  # Compute solution
  ts, sln = trajectory(r, **kwargs)

  # Pick the requested time series
  idim = {'X':0, 'Y':1, 'Z':2}
  xs = sln[idim[dim]]

  ws, psd = power_spectrum(ts, sln, dim, twin)
  py.figure(5)

  # Clip to requested frequency range
  iws = py.where((freqs[0] < ws) & (freqs[1] > ws))[0]

  if logy:
    psd = np.log2(psd[iws])
    ylbl = r'Spectral Density (dB / Hz)'
  else:
    psd = psd[iws]
    ylbl = r'Spectral Density (1/Hz)'

  py.plot(ws[iws], psd, lw=1.0, label='$r = %.1f$' % r, hold=True)
  py.xlabel(r'$\omega$')
  py.ylabel(ylbl)

  py.ion()
  py.show()
  py.draw()
# }}}

print ''' If you have imported this package in IDLE by opening this file and
choosing 'Run Module', the functions defined here will now be available to call
from the prompt. Try typing 'help(plot_timeseries)' to get information on how to
use them. If no plots appear after you run the code, try typing 'py.show()'.'''
