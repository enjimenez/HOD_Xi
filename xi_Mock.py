def xi_Mock (method, scatter=False, **inputs):
    import numpy as np
    import SPack as SP
    input_path = '../Data/G13_millennium'
    
    d = int(inputs['density_cut'])
    selection  = inputs['selection']
    fvir = int(inputs['rmax'])
    if scatter: spread = inputs['spread']
    
    if method == 'Full-HOD': ext_m = 'Full'
    if method == '2HOD': ext_m = '2HOD'
    if method == '4HOD': ext_m = '4HOD'
            
    data = np.load('../Data/G13_millennium/Mock Catalogues/%s/cat_d%i_%irvir.npy' %(sel, d, fvir))  ##
    x = data[:,0].astype('float64')  
    y = data[:,1].astype('float64')
    z = data[:,2].astype('float64')
    
    mask = data[:,4] != -1 # GALAXIES THAT DID NOT INCLUDED IN THE SAMPLE
    x = x[mask]
    y = y[mask]
    z = z[mask]
    N = len(x)

    Dmin = -2.0 # Min Log distance in Mpc
    Dmax =  2.0 # Max Log distance in Mpc
    NBIN =  40  # Number of bins
    Bw   =  (Dmax - Dmin) / NBIN # Bin width 
        
    iNBin = 50
    ilim  = 6

    DD_array = np.zeros(NBIN).astype('int32')
    RR_array = np.zeros(NBIN).astype('float64')

    def RR(r, Ngal, dr):
        Rmax = 10**(r+dr)
        Rmin = 10**r
        Vol = (4./3)*3.14*(Rmax**3 - Rmin**3)
        Den = Ngal/(500.**3)
        rr = Vol*Den*Ngal
        return rr
    
    bins = np.array([Dmin + Bw*i for i in range(NBIN)])
    bins_p = np.array([Dmin + Bw*i for i in range(NBIN)]) + Bw * 0.5
    
    SP.xi2(DD_array, x, y, z, N, Dmin, Dmax, NBIN, iNBin, ilim)
    RR_array = np.array([RR(r, N, Bw) for r in bins])

    xi = DD_array/RR_array -1
    np.savetxt('../Data/G13_millennium/xi_mocks/%s/xi_%s_d%i_%irvir.txt' %(sel, ext_m, d, fvir), np.array([bins, xi, DD, RR]).T) 