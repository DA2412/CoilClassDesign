# -*- coding: utf-8 -*-
"""
Created on Sun Jun  5 22:46:14 2022

@author: abate
"""
from matplotlib import cm
from copy import deepcopy
import cmath
import matplotlib.pyplot as plt
import time
import numpy as np
from scipy.special import p_roots

class Coil:
    path = []
    def __init__(self, current=1, path=None,discretization_length=0.1):
        self.current = current
        self.path = path
        self.discretization_length=discretization_length
        self.B = None
        return
    
    def Set_Current(self, current):
        '''Sets current with absolute value modulus and phase angle (in radians)'''
        self.current = current
        return    
    
    def LinearPath(self,pt1=(0, 0, 0), pt2=(0, 0, 1)):
        self.path = np.array([pt1, pt2])
        return 

    def define_Circular_Solenoid(self,R=0.1, length=1, turns=30, pts_per_turn=20):
        '''
        R: radius
        N: turns
        '''       
        return Coil.define_Elliptic_Solenoid(self,a=R, b=R, length=length, turns=turns, pts_per_turn=pts_per_turn)

    def define_Elliptic_Solenoid(self,a, b, length, turns, pts_per_turn=20):
        '''
        a: minor semi-axis, along x
        b: major semi-axis, along y 
        '''            
        pitch = length/turns
        t = np.linspace(0, 2 * np.pi * turns, pts_per_turn * turns)
        X = a * np.sin(t)
        Y = b * np.cos(t)
        Z = t / (2 * np.pi) * pitch
        self.path = np.array([X, Y, Z]).T       
        return 

    def define_SuperElliptic_Solenoid(self,a=0.1, b=0.2,exponent=2, length=1, turns=30, pts_per_turn=20):
        '''
        a: minor semi-axis, along x
        b: major semi-axis, along y 
        '''            
        pitch = length/turns
        t = np.linspace(0, 2 * np.pi * turns, pts_per_turn * turns)
        X = ((np.abs(np.cos(t))) ** (2 / exponent)) * a * np.sign(np.cos(t))
        Y = ((np.abs(np.sin(t))) ** (2 / exponent)) * b * np.sign(np.sin(t))
        Z = t / (2 * np.pi) * pitch
        self.path = np.array([X, Y, Z]).T         
        return 
    
    def plot3DCoil(self, ax=None):
        X = self.path[:,0]
        Y = self.path[:,1]
        Z = self.path[:,2]

        if ax is None:
            fig = plt.figure()
            ax1 = plt.axes(projection='3d')
        else:
            ax1 = ax
            
        ax1.scatter3D(X, Y, Z,edgecolor='k',facecolor='k',alpha=0.1)
        ax1.plot3D(X, Y, Z, color='k')
        ax1.set_xlabel('x')
        ax1.set_ylabel('y')
        ax1.set_zlabel('z')
        ax1.set_aspect('auto')

     # Create cubic bounding box to simulate equal aspect ratio
        max_range = np.array([X.max()-X.min(), Y.max()-Y.min(), Z.max()-Z.min()]).max()
        Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(X.max()+X.min())
        Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(Y.max()+Y.min())
        Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(Z.max()+Z.min())
        # Comment or uncomment following both lines to test the fake bounding box:
        for xb, yb, zb in zip(Xb, Yb, Zb):
           ax1.plot([xb], [yb], [zb], 'w')

        return ax1

    @property
    def compute_discrete_path(self):
        '''
        calculate end points of segments of discretized path
        approximate discretization lenghth is given by self.discretization_length
        elements will never be combined
        elements longer that self.dicretization_length will be divided into pieces
        :return: discretized path as m x 3 numpy array
        '''
        try:
            return self.discrete_path
        except AttributeError:
            pass

        self.discrete_path = deepcopy(self.path)
        for c in range(len(self.discrete_path)-2, -1, -1):
            # go backwards through all elements
            # length of element
            element = self.discrete_path[c+1]-self.discrete_path[c]
            el_len = np.linalg.norm(element)
            npts = int(np.ceil(el_len / self.discretization_length))  # number of parts that this element should be split up into
            if npts > 1:
                # element too long -> create points between
                # length of new sub elements
                sel = el_len / float(npts)
                for d in range(npts-1, 0, -1):
                    self.discrete_path = np.insert(self.discrete_path, c+1, self.discrete_path[c] + element / el_len * sel * d, axis=0)

        return self.discrete_path

    @property
    def compute_IdL_midpoints(self):
        '''
        calculate discretized path elements dL and their center point midpoint
        :return: numpy array with I * dL vectors, numpy array of midpoint vectors (center point of element dL)
        '''
        npts = len(self.compute_discrete_path)
        if npts < 2:
            print("discretized path must have at least two points")
            return

        IdL = np.array([self.compute_discrete_path[c+1]-self.compute_discrete_path[c] for c in range(npts-1)]) * self.current
        midpoint = np.array([(self.compute_discrete_path[c+1]+self.compute_discrete_path[c])*0.5 for c in range(npts-1)])

        return IdL, midpoint

    def computeB(self, points):
        """
        calculate magnetic field B at given points
        :param points: numpy array of n points (x y z)
        :return: numpy array of n vectors representing the B field at given points
        """
        IdL, midpoints = self.compute_IdL_midpoints
        print("total number of segments: {}".format(len(IdL)))
        print("number of field points: {}".format(len(points)))
        print("total number of calculations: {}".format(len(points)*len(IdL)))

        #B = np.array([c2.dB(r, IdL, midpoints) for r in points])
        #B = np.array([c2.dB(r, IdL, midpoints) for r in points])
        B = np.array([self.computedB(r, IdL, midpoints) for r in points])
        self.points = points
        self.B = B*1000 #in mT
        self.modB = np.linalg.norm(B, axis=1) *1000 #in mT
        return 
    
    def computedB(self,r, IdL, r1):
        '''
        calculate magnetic field B for one point r in space
        :param r: 3 space B calculation
        :param IdL:  length vectors times current
        :param r1: midpoints
        :return: numpy array of 3 component vector of B multiplied by 1e7
        '''
        mu0 = 4*np.pi*1e-7;

        # calculate law of biot savart for all current elements at given point r
        r2 = r - r1
        r25 = np.linalg.norm(r2, axis=1)**3
        r3 = r2 / r25[:, np.newaxis]

        cr = mu0/(4*np.pi)*np.cross(IdL, r3)

        # sum contributions 
        s = np.sum(cr, axis=0)

        return s

    def interpBsensors(self,sensorPosition):
        from scipy.interpolate import griddata
        Bpoint = griddata((self.points[:,0],self.points[:,1],self.points[:,2]),
                       self.B, sensorPosition, method='nearest')
 
        return np.linalg.norm(Bpoint)
    
    def plot3DBfield(self, ax=None):
        X = self.path[:,0]
        Y = self.path[:,1]
        Z = self.path[:,2]

        if ax is None:
            fig = plt.figure()
            ax1 = plt.axes(projection='3d')
        else:
            ax1 = ax
        colorB_scaled = self.modB/np.amax(self.modB)
        p=ax1.scatter3D(self.points[:,0],self.points[:,1],self.points[:,2],
                        c = colorB_scaled,cmap=cm.coolwarm )
        plt.title('scale factor = {0:.4f} mT'.format(np.amax(self.modB)))
        ax1.set_xlabel('x')
        ax1.set_ylabel('y')
        ax1.set_zlabel('z') 
        fig.colorbar(p,ax=ax1)      
        self.plot3DCoil(ax1)

    def plot2DBfield(self, ax=None):
        X = self.path[:,0]
        Y = self.path[:,1]
        Z = self.path[:,2]

        if ax is None:
            fig = plt.figure()
            ax1 = plt.axes()
        else:
            ax1 = ax

        if np.all(self.points[:,1] == self.points[1,1]):
            a = np.unique(self.points[:, 0])
            b = np.unique(self.points[:, 2])
            x = self.path[:, 0]
            y = self.path[:, 2]
            npoints = len(a)        
            B_a = self.B[:,0].reshape(npoints,npoints)
            B_b = self.B[:,2].reshape(npoints,npoints)
            a_label = 'x'
            b_label = 'z'
        elif np.all(self.points[:,2] == self.points[2,2]):
            a = np.unique(self.points[:, 0])
            b = np.unique(self.points[:, 1]) 
            x = self.path[:, 0]
            y = self.path[:, 1]
            npoints = len(a)        
            B_a = self.B[:,0].reshape(npoints,npoints)
            B_b = self.B[:,1].reshape(npoints,npoints)
            a_label = 'x'
            b_label = 'y'
        else:
            a = np.unique(self.points[:, 1])
            b = np.unique(self.points[:, 2])
            npoints = len(a)        
            B_a = self.B[:,1].reshape(npoints,npoints)
            B_b = self.B[:,2].reshape(npoints,npoints)
            x = self.path[:, 1]
            y = self.path[:, 2]
            a_label = 'y'
            b_label = 'z'

        modB_res = self.modB.reshape([len(a), len(b)]).T  

        sp1 = ax1.streamplot(
            a, b, B_a, B_b,
            density=2,
            color=modB_res,
            linewidth=modB_res*max(1/self.modB),
            cmap='coolwarm',
        )
        plt.plot(x,y,'ko',alpha=.1)
        sonda = plt.Rectangle((np.mean(a)-37.7*1e-3/2, np.mean(b)-42*1e-3/2),37.7*1e-3,42*1e-3, ec="black",alpha=.2)
        ax1.add_patch(sonda)
        ax1.axis('scaled')
        ax1.set_xlabel(a_label)
        ax1.set_ylabel(b_label)
        fig.colorbar(sp1.lines, ax=ax1, label='[mT]')

    def plotContourBfield(self,ncont, ax=None):
        X = self.path[:,0]
        Y = self.path[:,1]
        Z = self.path[:,2]

        if ax is None:
            fig = plt.figure()
            ax1 = plt.axes()
        else:
            ax1 = ax

        if np.all(self.points[:,1] == self.points[1,1]):
            a = np.unique(self.points[:, 0])
            b = np.unique(self.points[:, 2])
            x = self.path[:, 0]
            y = self.path[:, 2]
            npoints = len(a)        
            B_a = self.B[:,0].reshape(npoints,npoints)
            B_b = self.B[:,2].reshape(npoints,npoints)
            a_label = 'x'
            b_label = 'z'
        elif np.all(self.points[:,2] == self.points[2,2]):
            a = np.unique(self.points[:, 0])
            b = np.unique(self.points[:, 1]) 
            x = self.path[:, 0]
            y = self.path[:, 1]
            npoints = len(a)        
            B_a = self.B[:,0].reshape(npoints,npoints)
            B_b = self.B[:,1].reshape(npoints,npoints)
            a_label = 'x'
            b_label = 'y'
        else:
            a = np.unique(self.points[:, 1])
            b = np.unique(self.points[:, 2])
            npoints = len(a)        
            B_a = self.B[:,1].reshape(npoints,npoints)
            B_b = self.B[:,2].reshape(npoints,npoints)
            x = self.path[:, 1]
            y = self.path[:, 2]
            a_label = 'y'
            b_label = 'z'
                    
        # remove big values close to the wire
        B = self.B
        cutoff = 10*1/self.modB*self.modB
        B[self.modB > cutoff] = [np.nan,np.nan,np.nan]
        modB = np.linalg.norm(B, axis=1)

        # X = np.unique(self.points[:, 0])
        # Z = np.unique(self.points[:, 2])  
        # npoints = len(X)
        modB_res = modB.reshape([len(a), len(b)]).T  

        cs = ax1.contour(a, b, modB_res, ncont,cmap=cm.cool)
        ax1.clabel(cs)
        # sonda = plt.Rectangle((np.mean(a)-37.7*1e-3/2, np.mean(b)-42*1e-3/2),37.7*1e-3,42*1e-3, ec="black",alpha=.2)
        # ax1.add_patch(sonda)
        ax1.axis('scaled')
        ax1.set_xlabel('x')
        ax1.set_ylabel('y')
        plt.plot(x,y,'ko',alpha=.1)
        ax1.axis('scaled')
        

def gaussWeights(n):
    [x,w] = p_roots(n+1)
    return x,w

def computeInductance(target,source,n1target,n2source):
    '''
    https://en.wikipedia.org/wiki/Inductance#cite_note-21
    Neumann, F. E. (1846). "Allgemeine Gesetze der inducirten elektrischen Ströme". Annalen der Physik und Chemie (in German). Wiley. 143 (1): 
    # '''
    # n1g=7#target
    # n2g=5#source
    [x1,w1] = gaussWeights(n1target)
    [x2,w2] = gaussWeights(n2source)

    path1 = target
    path2 = source
    
    # circ = mat['coil_totoidal']
    # path1=circ
    # path2 = circ
    # path1 = c2.path
    # path2 = c2.path
    
    n_1 = len(path1)#target
    n_2 = len(path2)#source
    
    # if path1.shape[0] %2 != 0: 
    #     path1 = np.vstack((path1, path1[-1,:]))    
    # if path2.shape[0] %2 != 0: 
    #     path2 = np.vstack((path2, path2[-1,:]))
    # mid1 = np.array([(path1[c+1]+path1[c])*0.5 for c in range(n_1-1)])
    # mid2 = np.array([(path2[c+1]+path2[c])*0.5 for c in range(n_2-1)])
    
    # dL1 = np.array([path1[c+1]-path1[c] for c in range(n_1-1)]) 
    # dL1_len = np.linalg.norm(dL1,axis=1)
    # w_dL1_len = np.array([w1*dL1_len[c]/2 for c in range(len(dL1_len))])

    points_1_target = np.zeros(((n1target+1)*(path1.shape[0]-1),path1.shape[1]))
    for ii in range(n_1-1):
        dL1 = path1[ii+1]-path1[ii]
        dL1_len = np.linalg.norm(dL1)
        w1_dL1 = w1*dL1_len/2
        mid1 = (path1[ii+1]+path1[ii])*0.5
        p1_g =  np.array([mid1 + .5*hh*dL1 for hh in x1])
        ind = np.arange((ii)*(n1target+1),(ii+1)*(n1target+1))
        points_1_target[ind,:] = p1_g
    
    current = 1
    source = path2
    target = points_1_target
    ngsource = n2source
    A = fun_A_3D_thin_coil(source,current,target,ngsource)

    M = 0
    for ii in range(len(path1)-1):
        dL1 = path1[ii+1]-path1[ii]
        vers_t = dL1/np.linalg.norm(dL1)
        ind = np.arange((ii)*(n1target+1),(ii+1)*(n1target+1))
        A_ii = A[ind,:]
        temp_int = np.dot(A_ii,vers_t)
        res_int =  np.dot(temp_int,w1_dL1)
        M = M + res_int
        
    return M
        
def fun_A_3D_thin_coil(source,current,target,ngsource):
    mu0 = 4*np.pi*1e-7
    [xx,ww] = gaussWeights(ngsource)
    #Computation of magnetic vector potential
    n_target = len(target)   
    A = np.zeros(target.shape)
    
    for ii in range(len(source)-1):
        dL_src = source[ii+1]-source[ii]
        dL_src_len = np.linalg.norm(dL_src)
        w_dL_src = ww*dL_src_len/2
        mid_src = (source[ii+1]+source[ii])*0.5
        
        norm_r_rprime = np.zeros((n_target,ngsource+1))
        for jj in range(ngsource+1):                
            p1_g =  np.array(mid_src + .5*xx[jj]*dL_src)
            norm_r_rprime_jj = np.array(np.sqrt((p1_g[0]-target[:,0])**2 +
                                             (p1_g[1]-target[:,1])**2 +
                                             (p1_g[2]-target[:,2])**2))
            norm_r_rprime[::,jj] = norm_r_rprime_jj
     
        temp_int = 1/norm_r_rprime      
        res_int =  np.dot(temp_int,w_dL_src)
        res_int = mu0*current*res_int/(4*np.pi)         
        vers_t = dL_src/np.linalg.norm(dL_src)
        A_ii = np.array([res_int*vers_t[0], res_int*vers_t[1], res_int*vers_t[2]]).T
        A = A + A_ii
            
    return A