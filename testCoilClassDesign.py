# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 16:48:00 2022

@author: abate
"""
from CoilClassDesign import *
#parameters
radius = 70e-3
length = 140e-3
Nturns = 20
pts_per_turn=50

#defining grid of points for B
npoints = 100
xs = np.linspace(-radius+0.01, +radius-0.01, npoints)
ys = np.linspace(-radius+0.01, +radius-0.01, npoints)
zs = np.linspace(0+0.01, length-0.01, npoints)
PLANE_XZ = np.array([(x,0,z) for z in zs for x in xs])


##### circular coil
c1 = Coil()
c1.define_Circular_Solenoid(radius,length,Nturns,pts_per_turn)
c1.plot3DCoil()

c1.computeB(PLANE_XZ)

c1.plot3DBfield()

c1.plot2DBfield()

c1.plotContourBfield(50)

##### elliptic coil
c2 = Coil()
c2.define_Elliptic_Solenoid(a=radius, b=2*radius, length=length, turns=Nturns)
c2.plot3DCoil()
# c2.path
# c2.discretized_path
# c2.path is c2.dpath

c2.computeB(PLANE_XZ)

c2.plot3DBfield()

c2.plot2DBfield()

c2.plotContourBfield(50)

ys = np.linspace(-radius+0.01, +radius-0.01, npoints)
PLANE_XY= np.array([(x,y,length/2) for y in ys for x in xs])
c2.computeB(PLANE_XY)
c2.plot3DBfield()
c2.plot2DBfield()
c2.plotContourBfield(50)


n1target = 7
n2source = 5
L2 = computeInductance(c2.path, c2.path, n1target, n2source)
print('L(mH)={}'.format(L2*1000))

########## square coil superellipse
c3=Coil()
c3.define_SuperElliptic_Solenoid(a=radius,b=radius,exponent=4,length=length,turns=Nturns)      
c3.plot3DCoil()

c3.computeB(PLANE_XZ)
c3.plot3DBfield()
c3.plot2DBfield()
c2.plotContourBfield(50)

n1target = 7
n2source = 5
L3 = computeInductance(c3.path, c3.path, n1target, n2source)
print('L(mH)={0:.4}'.format(L3*1000))

########## square coil superellipse
c4=Coil()
c4.define_SuperElliptic_Solenoid(a=radius,b=radius,exponent=1,length=length,turns=Nturns)      
c4.plot3DCoil()

c4.computeB(PLANE_XZ)
c4.plot3DBfield()
c4.plot2DBfield()
c4.plotContourBfield(50)

n1target = 7
n2source = 5
L4 = computeInductance(c4.path, c4.path, n1target, n2source)
print('L(mH)={0:.4}'.format(L4*1000))

############# Axis ########################## ########################## 
sensorPosition = np.array([[0, 0, length/2]])
print('B axis={0:.4f}mT'.format(c1.interpBsensors(sensorPosition)))
print('B axis={0:.4f}mT'.format(c2.interpBsensors(sensorPosition)))
print('B axis={0:.4}mT'.format(c3.interpBsensors(sensorPosition)))
########################## ########################## ###################


# import scipy.io
# mat = scipy.io.loadmat('provaL.mat')
# circ = mat['coil_totoidal']

# computeInductance(circ,circ,7,5) #target, source, ntarget, nsource

from CoilClassDesign import *
## helmholtz coil

##xcoil
length = -9.526e-3
radius = 70.548e-3 #+ length/2
Nturns = 127
xcoil = Coil()
xcoil.define_Circular_Solenoid(R=radius, length=length, turns=Nturns, pts_per_turn=40)
xcoil.plot3DCoil()

xcoil.Set_Current(0.7)
#defining grid of points for B
npoints = 100
Sx = 36e-3
Sy = 5.3e-3
mid = radius/2
xs = np.linspace(-Sx/2, +Sx/2, npoints)
ys = np.linspace(-Sy/2, +Sy/2, npoints)
zs = np.linspace(0, mid, npoints)
PLANE_XZ = np.array([(x,0,z) for z in zs for x in xs])

ys = np.linspace(-radius+0.01, +radius-0.01, npoints)
PLANE_XY= np.array([(x,y,mid) for y in ys for x in xs])

xcoil.computeB(PLANE_XY)

xcoil.plot3DBfield()

xcoil.plot2DBfield()

xcoil.plotContourBfield(30)

sensorPosition = np.array([[0, 0, mid]])
print('----------------------Xcoil---------------')
print('B axis={0:.4f}mT'.format(xcoil.interpBsensors(sensorPosition)))

n1target = 7
n2source = 5
Lx = computeInductance(xcoil.path, xcoil.path, n1target, n2source)
print('L(mH)={}'.format(Lx*1000))


##ycoil
length = -12.99e-3
radius = 85.18e-3 #+ length/2
Nturns = 233
ycoil = Coil()
ycoil.define_Circular_Solenoid(R=radius, length=length, turns=Nturns, pts_per_turn=40)
ycoil.plot3DCoil()

ycoil.Set_Current(0.7)
#defining grid of points for B
npoints = 100
Sx = 40e-3
Sy = 5.8e-3
mid = radius/2
xs = np.linspace(-Sx/2, +Sx/2, npoints)
ys = np.linspace(-Sy/2, +Sy/2, npoints)
zs = np.linspace(0, mid, npoints)
PLANE_XZ = np.array([(x,0,z) for z in zs for x in xs])

ys = np.linspace(-radius+0.01, +radius-0.01, npoints)
PLANE_XY= np.array([(x,y,mid) for y in ys for x in xs])

ycoil.computeB(PLANE_XY)

ycoil.plot3DBfield()

ycoil.plot2DBfield()

ycoil.plotContourBfield(50)

sensorPosition = np.array([[0, 0, mid]])
print('----------------------Ycoil---------------')
print('B axis={0:.4f}mT'.format(ycoil.interpBsensors(sensorPosition)))

n1target = 7
n2source = 5
Ly = computeInductance(ycoil.path, ycoil.path, n1target, n2source)
print('L(mH)={}'.format(Ly*1000))

##zcoil
length = -6.0621e-3
radius = 41.417e-3 #+ length/2
Nturns = 53
zcoil = Coil()
zcoil.define_Circular_Solenoid(R=radius, length=length, turns=Nturns, pts_per_turn=40)
zcoil.plot3DCoil()

zcoil.Set_Current(0.5)
#defining grid of points for B
npoints = 100
Sx = 32e-3
Sy = 26e-3
mid = radius/2
xs = np.linspace(-Sx/2, +Sx/2, npoints)
ys = np.linspace(-Sy/2, +Sy/2, npoints)
zs = np.linspace(0, mid, npoints)
PLANE_XZ = np.array([(x,0,z) for z in zs for x in xs])

ys = np.linspace(-radius+0.01, +radius-0.01, npoints)
PLANE_XY= np.array([(x,y,mid) for y in ys for x in xs])

zcoil.computeB(PLANE_XY)

zcoil.plot3DBfield()

zcoil.plot2DBfield()

zcoil.plotContourBfield(30)

sensorPosition = np.array([[0, 0, mid]])
print('----------------------Zcoil---------------')
print('B axis={0:.4f}mT'.format(zcoil.interpBsensors(sensorPosition)))

n1target = 7
n2source = 5
Lz = computeInductance(zcoil.path, zcoil.path, n1target, n2source)
print('L(mH)={}'.format(Lz*1000))


Lxy = computeInductance(xcoil.path, ycoil.path, n1target, n2source)
print('Mxy(mH)={}'.format(Lxy*1000))


