def xi_SAM (**inputs):
    import numpy as np
    import SPack as SP
    input_path = '../Data/G13_millennium'
    
    d = int(inputs['density_cut'])
    selection  = inputs['selection']
        
    x = np.load('../Data/G13_millennium/x.npy')
    y = np.load('../Data/G13_millennium/y.npy')
    z = np.load('../Data/G13_millennium/z.npy')
    
    pos = int(d * 5)
    data  = np.loadtxt(input_path + '/density_data.txt')
    mass_cut, sfr_cut = data[pos][1], data[pos][2]
    
    if selection == 'Mstell': m_data = np.load(input_path + '/stellarmass.npy'); mask = m_data > mass_cut
    if selection == 'SFR': m_data = np.load(input_path + '/sfr.npy'); mask = m_data > sfr_cut
    
    Dmin = -2.0 # Min Log distance in Mpc
    Dmax =  2.0 # Max Log distance in Mpc
    NBIN =  40  # Number of bins
    Bw   =  (Dmax - Dmin) / NBIN # Bin width 
        
    iNBin = 50
    ilim  = 6
    bins_p = np.array([Dmin + Bw*i for i in range(NBIN)]) + Bw * 0.5

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
    N = len(x_m)

    bins = np.array([Dmin + Bw*i for i in range(NBIN)])
    SP.xi2(DD_array, x, y, z, N, Dmin, Dmax, NBIN, iNBin, ilim)
    RR_array = np.array([RR(r, N, Bw) for r in bins])

    xi = DD_array/RR_array -1
    np.savetxt('../Data/G13_millennium/TESTING.txt', np.array([bins_p, xi, DD_array, RR_array]).T)