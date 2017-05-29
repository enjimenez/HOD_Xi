import numpy as np
import os 


# ============= TOBIN AND SORT ================
path = "/home/esteban/Escritorio/Practica/Data/G14_millenium/"
os.chdir( path )
retval = os.getcwd()
print "Directory changed successfully %s" % retval

tt = np.argsort(np.loadtxt('dhaloid.csv'))
filenames = glob.glob('*csv')
for f in filenames:
    data = np.loadtxt(f)
    data = data[tt]
    np.save(f[:-4], data)
    
path = "/home/esteban/Escritorio/Practica/Data/L16_millenium/"
os.chdir( path )  
retval = os.getcwd()
print "Directory changed successfully %s" % retval
    
tt = np.argsort(np.loadtxt('dhaloid.csv'))
filenames = glob.glob('*csv')
for f in filenames:
    data = np.loadtxt(f)
    data = data[tt]
    np.save(f[:-4], data)    

# ============ ORDER BY TYPE GAL ===============

import glob

filenames = glob.glob('../Data/G13_millennium/*.npy')
fofid = np.load('../Data/G13_millennium/fofid.npy')
type_gal = np.load('../Data/G13_millennium/type.npy')
del fofid, type_gal
temp, idx = np.unique(fofid, return_index = True)
zeros = np.where(type_gal == 0)[0]
del temp
for f in filenames:
    print f
    data = np.load(f)
    for ix1, ix2 in zip(idx, zeros):
        data[ix1], data[ix2] = data[ix2], data[ix1]
    np.save(f[:-4], data)    
    
    
 # =========== ORDER BY MASS_INDEX ===============
 import glob
filenames = glob.glob('../Data/G13_millennium/*.npy')    
tt = np.argsort(mass_index)
del log_centralmvir, mass_index
for f in filenames:
    data = np.load(f)
    data = data[tt]
    np.save(f[:-4], data)   
    print "Done: ", f
    
    
DD = np.zeros(len(mass_index)).astype('int')
count,N = 0,0
for idx in range(-1,31):
    print idx
    indices = np.where((type_gal == 0) & (mass_index == idx))[0]
    for pos in indices:
        if pos == indices[-1]:
            print "ultimo pos"
        cond = True
        while cond:
            N += 1
            if pos+N == len(mass_index):
                break
            cond = (fofid[pos+N] == fofid[pos])
        if pos + N == len(mass_index): 
            print "sali del while y estoy en el for" 
            print "N del ultimo pos:",N,"pos: ", pos    
        DD[count:count+N] = np.arange(pos,pos+N) 
        count += N
        N = 0    

    
