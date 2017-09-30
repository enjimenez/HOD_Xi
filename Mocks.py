def ToMakeMock (method, scatter= False, **inputs):
    import numpy as np
    import SPack as SP
    
    m = inputs['model']
    type_cut = inputs['selection']
    if type_cut = 'Mstell': type_cut = 'Stellar Mass'
    d = inputs['density_cut']
    
    if 'scatter' in params: ext_scat = '_poisson'
    else: ext_scat = ''
    if 'spread' in params: ext_sp = '_ls'
    else: ext_sp = ''

    # Add an status # TODO
    
    model = m + "_millennium"
    path = "../Data/" + model
    # Loading
    x = np.load( path + '/x.npy' )
    y = np.load( path + '/y.npy' )
    z = np.load( path + '/z.npy' )
    if ab == 'G13': halomass = np.log10( np.load( 'centralmvir.npy' ) + 1e-10 ) + 10
    else: halomass = np.log10( np.load( 'mdhalo.npy' ) + 1e-10 ) + 10
    type_gal = np.load( 'type.npy' )

   
    idx = np.where( type_gal == 0 )[0] # Positions of the centrals

    # Centrals Data
    x_cen = x[idx]
    y_cen = y[idx]
    z_cen = z[idx]
    halomass_cen = halomass[idx]
    del x, y, z, halomass
    
    hod_path = path + '/hod_sams/' + type_cut
    cat_path = path + '/Mock Catalogues/' + type_cut
    
    if 'spread' in inputs: type_spread = inputs['spread']
    else: type_spread = -1
    def KvaluE (k, type_spread):
        if type_spread ==  'ls':
            if k <= n:
                k -= int((n-k)/2.)
                if k < 0: k = 0
            else: k += int((k-n)/2.)     
                    
        elif type_spread == 'ms':
            if k <= n: 
                k += int((n-k)/2.)
            else: k -= int((k-n)/2.)
                    
        return k
    
    # ========================================= FULL HOD ======================================
    if method == 'Full-HOD':
        ext_m = ''
        params = np.zeros((NBIN - 1, 2), dtype = np.float32)
        bins, HOD, Cent, Sats = np.loadtxt( hod_path + '/hod_d%i_30b' %d, unpack = True )
        
        for i in range( NBIN-1 ):
            y1, y2 = HOD[i], HOD[i+1]
            x1, x2 = bins[i], bins[i+1]
            m = ( y2 - y1 ) / ( x2 - x1)
            n = y1 - (m * x1)
            params[i] = np.array([m, n])

        mask = (halomass_cen >= bins[0]) & (halomass_cen <= bins[-1])
        halomass_cen = halomass_cen[mask]
        x_cen, y_cen, z_cen = x_cen[mask], y_cen[mask], z_cen[mask]
        Len_Cat = len(x_cen) # halo number
        
        mass_index = ((( halomass_cen - bins[0] ) / ( bins[-1] - bins[0] )) * (NBIN-1)).astype('int32')
        n_gal = np.zeros(Len_Cat).astype('float64')
        
        M, N = params[:,0], params[:,1]
        M = M.copy(order='C')
        N = N.copy(order='C')
        
        # It make the interpolation for the float halomass, mass index is to find the parameters M, N
        # that will be used
        SP.Ngals(n_gal, halomass_cen, mass_index, M, N, Len_Cat) 
        
        nn_gal = np.zeros(Len_Cat)
                
        for i in xrange(Len_Cat):
            n = n_gal[i]
                
            if scatter = True:
                spreading = -1
                rr = np.random.random_sample()
                cdf = np.exp(-n)
                while rr > cdf:
                    k += 1
                    cdf = np.exp(-n)*np.sum([n**j / math.factorial(j) for j in range(k+1)])
                
                
                k = KvaluE(k, type_spread)
                nn_gal[i] = k
                    
            else:    
                dec = n % 1
                rr = np.random.random_sample()
                if rr <= dec: nn_gal[i] = int(n) + 1
                else: nn_gal[i] = int(n)                                        
            
        del n_gal
        cat_data = np.zeros((Len_Cat,5), dtype = np.float32)
    
        for i,n in enumerate(nn_gal):
            cat_data[i] = np.array([x_cen[i], y_cen[i], z_cen[i], halomass_cen[i], n])
            
        np.save('cat_d%i%s%s%s' %(d, ext_m, ext_scat, ext_sp), cat_data)  
    
    # ========================================= 2 HOD ========================================= 
    if method == '2-HOD':
        ext_m = '_2HOD'
        bins, HOD, Cent, Sats = np.loadtxt( hod_path + '/hod_d%i_30b' %d, unpack = True )
        m_all, n_all = np.zeros(NBIN - 1).astype('float32'), np.zeros(NBIN - 1).astype('float32')
        m_cen, n_cen = np.zeros(NBIN - 1).astype('float32'), np.zeros(NBIN - 1).astype('float32')
        m_sat, n_sat = np.zeros(NBIN - 1).astype('float32'), np.zeros(NBIN - 1).astype('float32')
        
        for i in range( NBIN-1 ):
            y1_all, y2_all = HOD[i], HOD[i+1]
            y1_cen, y2_cen = Cent[i], Cent[i+1]
            y1_sat, y2_sat = Sats[i], Sats[i+1]
            x1, x2 = bins[i], bins[i+1]
        
            m_all[i] = ( y2_all - y1_all ) / ( x2 - x1)
            m_cen[i] = ( y2_cen - y1_cen ) / ( x2 - x1)
            m_sat[i] = ( y2_sat - y1_sat ) / ( x2 - x1)
    
            n_all[i] = y1_all - (m_all[i] * x1)
            n_cen[i] = y1_cen - (m_cen[i] * x1)
            n_sat[i] = y1_sat - (m_sat[i] * x1)
            
        mask = (halomass_cen >= bins[0]) & (halomass_cen <= bins[-1])
        halomass_cen = halomass_cen[mask]
        x_cen, y_cen, z_cen = x_cen[mask], y_cen[mask], z_cen[mask]
        Len_Haloes = len(x_cen)    
        
        mass_index = ((( halomass_cen - bins[0] ) / ( bins[-1] - bins[0] )) * (NBIN-1)).astype('int32')
        n_gal_all = np.zeros(Len_Haloes).astype('float64')
        n_gal_cen = np.zeros(Len_Haloes).astype('float64')
        n_gal_sat = np.zeros(Len_Haloes).astype('float64')
    
        SP.Ngals2(n_gal_all, n_gal_cen, n_gal_sat, halomass_cen, mass_index, m_all, n_all, m_cen, n_cen,
                     m_sat, n_sat, Len_Haloes)
        
        for i in xrange(Len_Cat):
            n1 = n_gal_all[i]
            n2 = n_gal_cen[i]
            n3 = n_gal_sat[i]
            
            rr1 = np.random.random_sample()
            rr2 = np.random.random_sample()
            rr3 = np.random.random_sample()
                
            if scatter = True:
                k1,k2,k3 = 0,0,0
                cdf1,cdf2,cdf3 = np.exp(-n1), np.exp(-n2), np.exp(-n3)
                
                while rr1 > cdf1:
                    k1 += 1
                    cdf1 = 1+np.exp(-n1)*np.sum([n1**j / math.factorial(j) for j in range(k1+1)])
                while rr2 > cdf2:
                    k2 += 1
                    cdf2 = 1+np.exp(-n2)*np.sum([n2**j / math.factorial(j) for j in range(k2+1)])    
                while rr3 > cdf3:
                    k3 += 1
                    cdf3 = 1+np.exp(-n3)*np.sum([n3**j / math.factorial(j) for j in range(k3+1)])   
                
                k1 = KvaluE(k1, type_spread)
                k2 = KvaluE(k2, type_spread)
                k3 = KvaluE(k3, type_spread)
               
                n_gal_all[i] = k1
                n_gal_cen[i] = k2
                n_gal_sat[i] = k3
            
            else:    
                dec1 = n1 % 1
                dec2 = n2 % 1
                dec3 = n3 % 1
                
                if rr1 <= dec1: n_gal_all[i] = int(n1) + 1
                else: n_gal_all[i] = int(n1)      
                if rr2 <= dec2: n_gal_cen[i] = int(n2) + 1
                else: n_gal_cen[i] = int(n2)      
                if rr3 <= dec3: n_gal_sat[i] = int(n3) + 1
                else: n_gal_sat[i] = int(n3)  
            
        cat_data = np.zeros((Len_Cat,7), dtype = np.float32)
        
        for i in range(Len_Haloes):
                cat_data[i] = np.array([x_cen[i], y_cen[i], z_cen[i], halomass_cen[i],
                                        n_gal_all[i], n_gal_cen[i], n_gal_sat[i]])
        
        np.save('cat_d%i%s%s%s' %(d, ext_m, ext_scat, ext_sp), cat_data)   
    
    # ========================================= 4 HOD =========================================
    if method == '4-HOD':
        ext_m = '_4HOD'
        
        bins, N_halos, N_cen, N_nocen, N_cen_sat, N_cen_nosat,
        N_sat_cen, N_sat_nocen =  np.loadtxt( hod_path + '/hod4_d%i_30b' %d, unpack = True )
    
        N_halos += 1e-10
        N_cen += 1e-10
        N_nocen += 1e-10
        N_cen_sat += 1e-10
        N_cen_nosat += 1e-10
        N_sat_cen += 1e-10
        N_sat_nocen += 1e-10
        
        
        fofid = np.load( path + '/fofid.npy' )
        
        for i in range( NBIN-1 ):
            y1_sat_cen, y2_sat_cen = N_sat_cen[i] / N_cen_sat[i], N_sat_cen[i+1] / N_cen_sat[i+1] 
            y1_sat_nocen, y2_sat_nocen = N_sat_nocen[i] / N_nocen_sat[i], N_sat_nocen[i+1] / N_nocen_sat[i+1]
            x1, x2 = bins[i], bins[i+1]
        
            m_sat_cen[i] = ( y2_sat_cen - y1_sat_cen ) / ( x2 - x1)
            m_sat_nocen[i] = ( y2_sat_nocen - y1_sat_nocen ) / ( x2 - x1)
    
    
            n_sat_cen[i] = y1_sat_cen - (m_sat_cen[i] * x1)
            n_sat_nocen[i] = y1_sat_nocen - (m_sat_nocen[i] * x1)
    
        mask = (halomass_cen >= bins[0]) & (halomass_cen <= bins[-1])
        halomass_cen = halomass_cen[mask]
        mass_index = ((( halomass_cen - bins[0] ) / ( bins[-1] - bins[0] )) * (NBIN-1)).astype('int32')
        x_cen, y_cen, z_cen = x_cen[mask], y_cen[mask], z_cen[mask]
        Len_Haloes = len(x_cen)    
        
        ff, counts_array = np.unique(fofid[mask], return_counts=True)
        del ff
        
        j = 0
        cat_data = np.zeros((Len_Haloes,5), dtype = np.float32)

        for i in range(Len_Haloes):  
            
            index = mass_index[j] # Go through all haloes
            counts = counts_array[i]
            cut_cen = N_cen[index] / N_halos[index] # tiene o no central?
            rr1 = np.random.rand(1)[0]
            
            if rr1 < cut_cen:
                cut_cen_sat = N_cen_sat[index] / N_cen[index]
                rr2 = np.random.rand(1)[0]
                
                if rr2 < cut_cen_sat:
                    m = m_sat_cen[index]
                    n = n_sat_cen[index] 
                    Nsat = (m * halomass_cen[i]) + n;
                    
                    if scatter = True:
                        rr = np.random.random_sample()
                        cdf = np.exp(-Nsat)
                        k = 0
                        while rr > cdf:
                            k += 1
                            cdf = np.exp(-Nsat)*np.sum([n**p / math.factorial(p) for p in range(k+1)])
                    
                   
        
                        Nsat = KvaluE(k, type_spread)
        
                    else:    
                        dec = Nsat % 1
                        rr = np.random.random_sample()
                        if rr <= dec: Nsat = int(Nsat) + 1
                        else: Nsat[i] = int(Nsat)                                        
            
    
                    NGAL = Nsat + 1
                    cat_data[i] = np.array([x[j], y[j], z[j], index, NGAL])
                    
                else:
                    cat_data[i] = np.array([x[j], y[j], z[j], index, 1])
                    
            else:
                cut_nocen_sat = N_nocen_sat[index] / N_nocen[index]
                rr3 = np.random.rand(1)[0]
                
                if rr3 < cut_nocen_sat:
                    m = m_sat_nocen[index]
                    n = n_sat_nocen[index] 
                    Nsat = (m * halomass_cen[i]) + n;
                    
                    if scatter = True:
                        rr = np.random.random_sample()
                        cdf = np.exp(-Nsat)
                        k = 0
                        while rr > cdf:
                            k += 1
                            cdf = np.exp(-Nsat)*np.sum([n**p / math.factorial(p) for p in range(k+1)])
                    
                        Nsat = KvaluE(k, type_spread)
        
                    else:    
                        dec = Nsat % 1
                        rr = np.random.random_sample()
                        if rr <= dec: Nsat = int(Nsat) + 1
                        else: Nsat[i] = int(Nsat)           
                        
                    cat_data[i] = np.array([x[j], y[j], z[j], index, Nsat])
                    
                else:
                    cat_data[i] = np.array([x[j], y[j], z[j], index, 0])
                                    
            j += counts
            #if j % display_step == 0:
            #    print "Step: %i Time: %3fs" % (j, (time.time() - start_time))
                
        np.save('cat_d%i%s%s%s' %(d, ext_m, ext_scat, ext_sp), cat_data) 
