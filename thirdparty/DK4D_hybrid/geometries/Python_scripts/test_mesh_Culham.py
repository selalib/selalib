import matplotlib.pyplot as mpp
import numpy as np
from Cequilibrum import Culham_computation
from pigasus.fit.curfit import curfit
from interpolator_Culham import curve_fit
from caid.cad_geometry import cad_geometry, cad_nurbs
from caid.utils.coons import coons, coonsInitialize

#--- Input variables ---
a = 1.
R0 = 5.
N = 100
B0 = 1.
q0 = 0.8
qa = 0.7
alpha = 1
p0 = 10**5
Nomega = 64
Nr = 10
pa = p0/10
h = a/N
r = h*np.arange((N+1))
q = q0 + (qa - q0)*(r/a)**2
p = (p0 - pa)*(1-(r/a)**2)**alpha + pa
Ea = 0.25
Ta = 0.1

#--------------------------------------
#--> Compute the Culham curves
#--------------------------------------
r     = h*np.arange(np.size(q))
omega = np.linspace(0, 2*np.pi, num=Nomega+1)
[Diff_Delta, Delta, E, Diff_E, T, Diff_T, f, g, C_E, C_T] = Culham_computation(q, p, R0, B0, Ea, Ta, h, r)
P = r**3 / (8 * R0**2) + r * Delta / (2*R0)

print "---> Compute (X,Y) positions"
Ncurves      = int(np.size(q)/Nr)
indx_rpoint  = range(Ncurves-1, np.size(q), Ncurves)
Nomega       = np.size(omega)
Xpoint       = np.zeros((Ncurves,Nomega),dtype=np.float64)
Ypoint       = np.zeros((Ncurves,Nomega),dtype=np.float64)
Zpoint       = np.zeros((Ncurves,Nomega),dtype=np.float64)
icurve = 0
for ipoint in indx_rpoint:
    for iomega in range(0,Nomega,1):
        omega_pt = omega[iomega]
        Xpoint[icurve,iomega] = (- E + r - P)[ipoint] * np.cos(omega_pt) - Delta[ipoint] + T[ipoint] * np.cos(2*omega_pt)
        Ypoint[icurve,iomega] = (r + E - P)[ipoint] * np.sin(omega_pt) - T[ipoint] * np.sin(2 * omega_pt)
    #end for
    icurve = icurve + 1
#end for

choice_algo = 1
if (choice_algo == 0): 
    #----------------------------------------------------
    #--> Construct the mesh by using a Coons algorithm
    #----------------------------------------------------
    #--> Construct the last closed surface
    print "---> Fit the last closed Culham curve"
    px      = 3
    nx      = 63  # nx+1 must be divisible by for to split the curve in 4
    METHODS = ["uniform", "chord", "centripetal"]
    method  = METHODS[1]
    alpha   = 1.
    geo     = curve_fit(Xpoint[-1,:],Ypoint[-1,:],nx,px,method,alpha,verbose=True)

    #--> Split the curve in 4 patchs
    geo.split(0,0.25,axis=0)
    geo.split(1,0.5,axis=0)
    geo.split(2,0.75,axis=0)

    #--> Apply the coons algorithm
    c0 = geo[0]
    c0.reverse(0)  # c0 and c1 must begin with the same point to give a direction (then not necessary for the other faces)
    c1 = geo[1]
    c2 = geo[2]
    c3 = geo[3]
    curves  = [[c0,c2],[c1,c3]]
    tol     = 1.e-7
    curves  = coonsInitialize(curves, tol=tol)
    nrb     = coons(curves)
    cad_nrb = cad_nurbs(nrb.knots, nrb.points, weights=nrb.weights)
    geo_out = cad_geometry()
    geo_out.append(cad_nrb)
    geo_out.plotMesh(MeshResolution=3)
    mpp.show()

elif ( choice_algo == 1 ):
    #------------------------------------------------------------
    #--> Construct the mesh aligned to the magnetic field lines
    #------------------------------------------------------------
    px       = 3
    nx       = 63  
    METHODS  = ["uniform", "chord", "centripetal"]
    method   = METHODS[1]
    alpha    = 1.
    list_crv = []
    for icurve in range(0,Ncurves,1):
        print "---> Fit curve " + str(icurve)
        geo = curve_fit(Xpoint[icurve,:],Ypoint[icurve,:],nx,px,method,alpha,verbose=False)
        list_crv.append(geo)
    #end for

    p1 = px
    p2 = px
    n1 = len(list_crv)
    n2 = list_crv[0][0].shape[0]

    t1 = np.array([0.]*(p1+1) + list(np.linspace(0.,1.,n1-p1-1+2)[1:-1]) + [1.]*(p1+1))
    t2 = list_crv[0][0].knots[0]

    knots   = [t1,t2]
    points  = np.zeros((n1,n2,3))
    weights = np.zeros((n1,n2))

    # ... loop over levels
    for i in range(0, n1):
        nrb = list_crv[i][0]
        # ... loop over curve's control points
        for j in range(0, n2):
            points[i,j,:] = nrb.points[j,:]
            weights[i,j]  = nrb.weights[j]
        # ...
    # ...

    cad_nrb = cad_nurbs(knots, points, weights=weights)
    cad_nrb.rational = True
    cad_nrb.orientation = [-1,1,1,-1]

    geo_out = cad_geometry()
    geo_out.append(cad_nrb)

    #---> How to find the internal face:
    #       geo_out[0].plotBoundariesInfo()
    #       mpp.show()

    face = 1   # Face corresponding to the internal boundary
    geo_new = geo_out.to5patchs(face)
    fig = geo_new.plotMesh(MeshResolution=3)
    mpp.axis('equal')
    mpp.show()

#endif    

