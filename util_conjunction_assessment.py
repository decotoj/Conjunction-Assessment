# Conjunction Assessment 
# Jake Decoto (decotoj@gmail.com)
# Original: September 22, 2021
# Last Update: October 14, 2022

from datetime import timedelta
import numpy as np
import time
from collections import defaultdict
from numba import njit

# try:
#     import catapult.util_spacetrack_tles as strack
# except:
import util_spacetrack_tles as strack # standalone version for demo (not in catapult library)

@njit
def helper(n, x, y, z, rt, a, b):
    for j in range(len(x)):
        q = np.searchsorted(x,x[j]+rt, 'left') 
        s = [[n[j], n[k]] for k in range(j+1,q) if abs(y[j]-y[k]) < rt and abs(z[j]-z[k])<rt]  
        for ss in s:
            if ss[0] != ss[1]:
                a = np.append(a, min(ss))
                b = np.append(b, max(ss))

    return a[1:], b[1:]

def conjunction_run(cat, rt, dt0, steps, stepSize, maxEphemBlockSize=3600):

    conj = defaultdict(lambda: defaultdict(lambda:100000))

    # Step Through Blocks of Epochs
    pairs = set()
    i1 = 0
    i2 = min(steps, maxEphemBlockSize)
    while i1 < steps:

        print(f'Step {i1} to {i2} of {steps}')

        # epochs
        epochs = [dt0 + timedelta(seconds=i*stepSize) for i in range(i1,i2)]

        # Build TEME Frame Ephemeris
        t1 = time.time()
        teme, scc = strack.buildTLEEphem(cat, epochs)
        scc = list(scc)
        print(f'A: Build TEME Ephem: Run Time {round(time.time()-t1,2)} seconds')

        # Execute Component Sort Filter
        t1 = time.time()
        for i in range(len(epochs)):
            x = teme[:,i,:][:,0]
            n = x.argsort()
            x = x[n]
            y = teme[:,i,:][n,1]
            z = teme[:,i,:][n,2]

            a,b = helper(n, x, y, z, rt, np.array([0]), np.array([0]))
            
            pairs = set.union(pairs, set([(a[k],b[k]) for k in range(len(a))] ))
        
        i1 = i2
        i2 = min(steps, i2+maxEphemBlockSize)

        print(f'B: Component Sort Filter: Run Time {round(time.time()-t1,2)} seconds')

        # Step Through Valid Pair Indices
        t1 = time.time()
        for p in pairs:
            r = np.linalg.norm(np.subtract(teme[p[0],:,:], teme[p[1],:,:]), axis=1)

            if conj[(scc[p[0]], scc[p[1]])]['Range (km)'] > min(r):
                conj[(scc[p[0]], scc[p[1]])] = {'Epoch': epochs[np.argmin(r)], 'Range (km)': min(r)} 

        print(f'C: Process Conjunction Candidates {round(time.time()-t1,2)} seconds')

        del teme, scc

    return conj