def xi_or(x,y,z, type_cuts, density_cuts, input_path):
    import numpy as np
    import SPack as SP
    
    Dmin = -2.0 # Min Log distance in Mpc
    Dmax =  2.0 # Max Log distance in Mpc
    NBIN =  40  # Number of bins
    Bw   =  (Dmax - Dmin) / NBIN # Bin width 

    bins = np.array([Dmin + Bw*i for i in range(NBIN)])
    bins_p = np.array([Dmin + Bw*i for i in range(NBIN)]) + Bw * 0.5

    iNBin = 50
    ilim  = 6
    
    def RR(r, Ngal, dr):
        Rmax = 10**(r+dr)
        Rmin = 10**r
        Vol = (4./3)*3.14*(Rmax**3 - Rmin**3)
        Den = Ngal/(500.**3)
        rr = Vol*Den*Ngal
        return rr
    
    if (type(type_cut) == list) or (type(density_cut) == list):
        for Type in type_cut:
            for Den in density_cut:
                if Type == 'Mstell' and Den <= 4 and Den > 0:
                    pos = int(Den * 5)
                    data  = np.loadtxt(input_path + '/density_data.txt')
                    Log_Den, mass_cut = data[pos][0], data[pos][1]
                    
                    m_data = np.load(input_path + '/stellarmass.npy')
                    mask = m_data > mass_cut
                    
                    x_m = x[mask]
                    y_m = y[mask]
                    z_m = z[mask]
                    N = len(x_m)
                
                    DD = np.zeros(NBIN).astype('int32')
                    SP.xi(DD, x_m, y_m, z_m, N, Dmin, Dmax, NBIN, iNBin, ilim)
                    RR_array = np.array([RR(r, N, Bw) for r in bins])
                    return DD, RR
                    
                if Type == 'SFR' and Den <= 4 and Den > 0:
                    pos = int(Den * 5)
                    data  = np.loadtxt(input_path + '/density_data.txt')
                    Log_Den, sfr_cut = data[pos][0], data[pos][2]
                    
                    m_data = np.load(input_path + '/sfr.npy')
                    mask = m_data > sfr_cut
                    
                    x_m = x[mask]
                    y_m = y[mask]
                    z_m = z[mask]
                    N = len(x_m)
                    
                    DD = np.zeros(NBIN).astype('int32')
                    SP.xi2(DD, x_m, y_m, z_m, N, Dmin, Dmax, NBIN, iNBin, ilim)
                    RR = np.array([RR(r, N, Bw) for r in bins])
                    return DD, RR
    
    else:
        if (type_cut == 'Mstell') and (density_cut <= 4) and (density_cut > 0):
            
            pos = int(density_cut * 5)
            data  = np.loadtxt(input_path + '/density_data.txt')
            Log_Den, mass_cut = data[pos][0], data[pos][1]
            
            m_data = np.load(input_path + '/stellarmass.npy')
            mask = (m_data > mass_cut)

            x_m = x[mask]
            y_m = y[mask]
            z_m = z[mask]
            N = len(x_m)
            
            DD_array = np.zeros(NBIN).astype('int32')
            SP.xi2(DD_array, x_m, y_m, z_m, N, Dmin, Dmax, NBIN, iNBin, ilim)
            RR_array = np.array([RR(r, N, Bw) for r in bins])
            return DD_array, RR_array
    
        if type_cut == 'SFR' and density_cut <= 4 and density_cut > 0:
           
            pos = int(density_cut * 5)
            data  = np.loadtxt(input_path + '/density_data.txt')
            Log_Den, sfr_cut = data[pos][0], data[pos][2]
            
            m_data = np.load(input_path + '/sfr.npy')
            mask = m_data > sfr_cut
            
            x_m = x[mask]
            y_m = y[mask]
            z_m = z[mask]
            N = len(x_m)
            
            DD = np.zeros(NBIN).astype('int32')
            SP.xi2(DD, x_m, y_m, z_m, N, Dmin, Dmax, NBIN, iNBin, ilim)
            RR = np.array([RR(r, N, Bw) for r in bins])
            return DD, RR
            
        else:
            print "Error: select the correct parameters!"
    
def xi_ll(x,y,z, type_cuts, density_cuts, input_path):    
    import numpy as np
    import SPack as SP
    
    Dmin = -2.0 # Min Log distance in Mpc
    Dmax =  2.0 # Max Log distance in Mpc
    NBIN =  40  # Number of bins
    Bw   =  (Dmax - Dmin) / NBIN # Bin width 

    bins = np.array([Dmin + Bw*i for i in range(NBIN)])
    bins_p = np.array([Dmin + Bw*i for i in range(NBIN)]) + Bw * 0.5

    iNBin = 50
    ilim  = 6
    
    def RR(r, Ngal, dr):
        Rmax = 10**(r+dr)
        Rmin = 10**r
        Vol = (4./3)*3.14*(Rmax**3 - Rmin**3)
        Den = Ngal/(500.**3)
        rr = Vol*Den*Ngal
        return rr

    if (type(type_cuts) == list) or (type(density_cuts) == list):
        Len_densities, Len_types = 1, 1
        
        if type(density_cuts) == list: Len_densities = len(density_cuts)
        if type(type_cuts) == list: Len_types = len(type_cuts)
        
        DD = np.zeros((Len_types, Len_densities, NBIN)).astype('int32')
        RR_array = np.zeros((Len_types, Len_densities, NBIN)).astype('float64')
        
        if Len_densities == 1: density_cuts = [density_cuts]
        if Len_types == 1: type_cuts = [type_cuts]
        
        for i, Type in enumerate(type_cuts):
            for j, Den in enumerate(density_cuts):
                if Type == 'Mstell' and Den <= 4 and Den > 0:
                    pos = int(Den * 5)
                    data  = np.loadtxt(input_path + '/density_data.txt')
                    Log_Den, mass_cut = data[pos][0], data[pos][1]
                    
                    m_data = np.load(input_path + '/stellarmass.npy')
                    mask = m_data > mass_cut
                    
                    x_m = x[mask]
                    y_m = y[mask]
                    z_m = z[mask]
                    N = len(x_m)
                
                    #DD = np.zeros(NBIN).astype('int32')
                    SP.xi2(DD[i][j], x_m, y_m, z_m, N, Dmin, Dmax, NBIN, iNBin, ilim)
                    print i,j
                    RR_array[i][j] = np.array([RR(r, N, Bw) for r in bins])
                    print "Sali!"

                    #np.savetxt(output_path + '/%s_d%i.txt' %(type_cut, density_cut), np.array([bins_p, DD, RR])
                    
                if Type == 'SFR' and Den <= 4 and Den > 0:
                    pos = int(Den * 5)
                    data  = np.loadtxt(input_path + '/density_data.txt')
                    Log_Den, sfr_cut = data[pos][0], data[pos][2]
                    
                    m_data = np.load(input_path + '/sfr.npy')
                    mask = m_data > sfr_cut
                    
                    x_m = x[mask]
                    y_m = y[mask]
                    z_m = z[mask]
                    N = len(x_m)
                    
                    #DD = np.zeros(NBIN).astype('int32')
                    SP.xi2(DD[i][j], x_m, y_m, z_m, N, Dmin, Dmax, NBIN, iNBin, ilim)
                    print i,j
                    RR_array[i][j] = np.array([RR(r, N, Bw) for r in bins])

        return DD, RR_array
    else:
        if (type_cuts == 'Mstell') and (density_cuts <= 4) and (density_cuts > 0):
            
            pos = int(density_cuts * 5)
            data  = np.loadtxt(input_path + '/density_data.txt')
            Log_Den, mass_cut = data[pos][0], data[pos][1]
            
            m_data = np.load(input_path + '/stellarmass.npy')
            mask = (m_data > mass_cut)

            x_m = x[mask]
            y_m = y[mask]
            z_m = z[mask]
            N = len(x_m)
            
            DD = np.zeros(NBIN).astype('int32')
            SP.xi2(DD, x_m, y_m, z_m, N, Dmin, Dmax, NBIN, iNBin, ilim)
            RR = np.array([RR(r, N, Bw) for r in bins])
            return DD, RR
    
        if type_cuts == 'SFR' and density_cuts <= 4 and density_cuts > 0:
           
            pos = int(density_cuts * 5)
            data  = np.loadtxt(input_path + '/density_data.txt')
            Log_Den, sfr_cut = data[pos][0], data[pos][2]
            
            m_data = np.load(input_path + '/sfr.npy')
            mask = m_data > sfr_cut
            
            x_m = x[mask]
            y_m = y[mask]
            z_m = z[mask]
            N = len(x_m)
            
            DD = np.zeros(NBIN).astype('int32')
            SP.xi2(DD, x_m, y_m, z_m, N, Dmin, Dmax, NBIN, iNBin, ilim)
            RR = np.array([RR(r, N, Bw) for r in bins])
            return DD, RR
            
        else:
            print "Error: select the correct parameters!"