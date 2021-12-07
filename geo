import sympy as sym
import numpy as np
import scipy.integrate as sc

def Christoffel2nd_sphere():
    from sympy.abc import u,v
    from sympy import tan, cos ,sin
    # found using sympy.diffgeom
    return sym.Matrix([[(0, -2*tan(v)),  (0, 0)],[(sin(v)*cos(v), 0), (0, 0)]])
 
def display_sphere_with_geodesic(u0,s0,s1,ds): 
    C = Christoffel2nd_sphere()
    X = solve(C,u0,s0,s1,ds)
    import matplotlib.pylab as plt
    from mpl_toolkits.mplot3d import Axes3D
    N = X[:,0].shape[0]
    u,v = plt.meshgrid(np.linspace(0,2*np.pi,N),np.linspace(0,2*np.pi,N))
    x = np.cos(u)*np.cos(v)
    y = np.sin(u)*np.cos(v)
    z = np.sin(v)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.view_init(elev=50, azim=-120)
    # use transparent colormap
    import matplotlib.cm as cm
    theCM = cm.get_cmap()
    theCM._init()
    alphas = -.5*np.ones(theCM.N)
    theCM._lut[:-3,-1] = alphas
    ax.plot_surface(x,y,z,linewidth=0,cmap=theCM)
    plt.hold('on')

    # plot the parametrized data on to the sphere
    u,v = X[:,0], X[:,2]
    x = np.cos(u)*np.cos(v)
    y = np.sin(u)*np.cos(v)
    z = np.sin(v)

    ax.plot(x,y,z,'--r')
    from math import pi
    s1_ = s1/pi
    fig.suptitle('$s\in[%.1f\, , \,%2.1f\pi]$'%(s0,s1_))
    plt.show()
def f(y,s,C,u,v):
    y0 = y[0] # u
    y1 = y[1] # u'
    y2 = y[2] # v
    y3 = y[3] # v'
    dy = np.zeros_like(y)
    dy[0] = y1
    dy[2] = y3

    C = C.subs({u:y0,v:y2})

    dy[1] = -C[0,0][0]*dy[0]**2 -\
           2*C[0,0][1]*dy[0]*dy[2] -\
             C[0,1][1]*dy[2]**2
    dy[3] = -C[1,0][0]*dy[0]**2 -\
           2*C[1,0][1]*dy[0]*dy[2] -\
             C[1,1][1]*dy[2]**2
    return dy

def solve(C,u0,s0,s1,ds):
    s = np.arange(s0,s1+ds,ds)
    from sympy.abc import u,v
    return sc.odeint(f,u0,s,args=(C,u,v))

from math import pi
u0 = [0,.1,0,.1]
s0 = 0
s1 = 18*pi
ds = 0.15
display_sphere_with_geodesic(u0,s0,s1,ds)
