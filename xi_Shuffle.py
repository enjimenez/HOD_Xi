def xi_Shuffle (**inputs):
    import numpy as np
    import SPack as SP
    input_path = '../Data/G13_millennium'
    
    d = int(inputs['density_cut'])
    selection  = inputs['selection']
    fvir = int(inputs['rmax'])
    
    pos = int(d * 5)
    data  = np.loadtxt(input_path + '/density_data.txt')
    mass_cut, sfr_cut = data[pos][1], data[pos][2]
    
    if selection == 'Mstell': m_data = np.load(input_path + '/Shuffle/stellarmass_shuffle.npy'); mask = m_data > mass_cut; sel = 'Stellar Mass'
    if selection == 'SFR': m_data = np.load(input_path + '/Shuffle/sfr_shuffle.npy'); mask = m_data > sfr_cut; sel = 'SFR'
                
    cat_data = np.load(input_path + '/Shuffle/%s/cat_shuffle_nfw_%irvir.npy' %(sel, fvir))  ##
    x = cat_data[:,0].astype('float64')  
    y = cat_data[:,1].astype('float64')
    z = cat_data[:,2].astype('float64')
            
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
    
    x = x[mask]
    y = y[mask]
    z = z[mask]
    N = len(x)

    bins = np.array([Dmin + Bw*i for i in range(NBIN)])
    bins_p = np.array([Dmin + Bw*i for i in range(NBIN)]) + Bw * 0.5
    
    SP.xi2(DD_array, x, y, z, N, Dmin, Dmax, NBIN, iNBin, ilim)
    RR_array = np.array([RR(r, N, Bw) for r in bins])

    xi = DD_array/RR_array -1
    np.savetxt(input_path + '/xi_shuffle/%s/xi_shuffle_d%i_%irvir.txt' %(sel, d, fvir), np.array([bins_p, xi, DD_array, RR_array]).T) 