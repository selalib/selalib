#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 15:41:14 2019

@author: kako
"""

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 16:02:09 2019

@author: kako
"""
import numpy as np
import matplotlib.pyplot as plt

folder3000 = '/home/kako/run_dirs/gempic/strongA/weibel_191203/pyt/'
folder = '/home/kako/run_dirs/gempic/strongA/iaw_191213/'

# Reference solution
#ref = np.loadtxt(folder3000+'disgrad_f3_s2_thdiag.dat')
ref = np.loadtxt('/home/kako/run_dirs/gempic/strongA/weibel_191203/ref/strong_weibel_dg_f3_s3_thdiag.dat')
ref2 = np.loadtxt('/home/kako/run_dirs/gempic/strongA/iaw_191213/symplectic_f-1_s2_thdiag.dat')

sdeg_vec = [1, 2, 3, 4 ]
dfem_vec = [1,2,3]#[-1,1,2,3,4,5]
propagators = ['disgrad', 'symplectic']
propagators = ['symplectic']
def check_conservation(data ):
    print('Conservation Gauss:', np.max(np.abs(data[:,11]-data[0,11])))
    print('Conservation energy:', np.max(np.abs(data[:,8]-data[0,8])))
    print('Momentum error', np.max(np.abs(data[:,9])), np.max(np.abs(data[:,10])))
    #plt.semilogy(data[:,0],data[:,3])

#dfem_vec = [3,4]
#sdeg_vec = [3]
#dfem_vec = [1,2,3,4]
#propagators = ['symplectic']
#sdeg_vec = [3]
#dfem_vec = [3]
ran = range(0,40001)
for sdeg in sdeg_vec :
    for propagator in propagators :
        for dfem in dfem_vec :
            prefix = propagator +'_f'+str(dfem)+'_s' + str(sdeg)
            data = np.loadtxt(folder + prefix +'_thdiag.dat')
            plt.semilogy(data[ran,0],data[ran,6])
            print(prefix)
            #check_conservation(data)
       # plt.semilogy(ref[ran,0],ref[ran,3])#,ref2[ran,0],ref2[ran,3])
        plt.semilogy(ref2[ran,0],ref2[ran,6])
        plt.legend(dfem_vec)
        plt.title(str(sdeg)+' '+propagator)
        plt.show()

err = np.zeros([6,4])
ind = 0
for dfem in dfem_vec :
    for propagator in propagators :
        ins = 0
        for sdeg in sdeg_vec :
            prefix = propagator +'_f'+str(dfem)+'_s' + str(sdeg)
            data = np.loadtxt(folder + prefix +'_thdiag.dat')
            plt.semilogy(data[ran,0],data[ran,1])
            print(prefix, np.max(np.abs(data[:,6]-ref2[:,6])))
            err[ind,ins] = np.sqrt(np.sum((data[:,1]-ref2[:,1])**2))#np.max(np.abs(data[:,1]-ref2[:,1]))
            ins = ins+1
            #check_conservation(data)
       # plt.semilogy(ref[ran,0],ref[ran,3])#,ref2[ran,0],ref2[ran,3])
        #plt.semilogy(ref2[ran,0],ref2[ran,6])
        plt.legend(sdeg_vec)
        plt.title(str(dfem)+' '+propagator)
        plt.show()
    ind = ind +1