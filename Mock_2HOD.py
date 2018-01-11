def ToMakeMock (scatter= False, **inputs):
    import numpy as np
    import SPack as SP
    from scipy.interpolate import interp1d
    import math
    
    m = inputs['model']
    selection = inputs['selection']
    if selection == 'Mstell': sel = 'Stellar Mass'
    if selection == 'SFR': sel = 'SFR'
    
    d = inputs['density_cut']
    fvir = int(inputs['rmax'])
    rhos = np.array([-3.5,-3.0,-2.5,-2.0,-1.5])
    s_density = 10**rhos[d-1]
    
    if scatter == True: ext_scat = '_poisson'; s_scatter = 'Poisson'
    else: ext_scat = ''; s_scatter = 'No Scatter'
    if 'spread' in inputs: ext_sp = '_' + inputs['spread']  #TODO
    else: ext_sp = ''

    # Add an status
    print("--- Building Mock Catalogue ---")
    print "Density: %.5f" %s_density
    print "Method: 2-HOD"
    print "Scatter: %s " %s_scatter
    print "Rmax: %i r200" %fvir
    
    model = m + "_millennium"
    path = "../Data/" + model
    # Loading
    x = np.load( path + '/x.npy' )
    y = np.load( path + '/y.npy' )
    z = np.load( path + '/z.npy' )
    if m == 'G13': halomass = np.log10( np.load(path + '/centralmvir.npy' ) + 1e-10 ) + 10
    else: halomass = np.log10( np.load(path + '/mdhalo.npy' ) + 1e-10 ) + 10
    type_gal = np.load(path + '/type.npy' )

   
    idx = np.where( type_gal == 0 )[0] # Positions of the centrals

    # Centrals Data
    x_cen = x[idx]
    y_cen = y[idx]
    z_cen = z[idx]
    halomass_cen = halomass[idx]
    del x, y, z, halomass
    
    hod_path = path + '/hod_sams/' + sel
    cat_path = path + '/Mock Catalogues/' + sel
    
    if 'spread' in inputs: type_spread = inputs['spread']
    else: type_spread = -1
    def KvaluE (k, type_spread, n):
        if type_spread ==  "ms":
            if k <= n:
                k = k - int((n-k)/2.)
                if k < 0: k = 0
            else: k = k + int((k-n)/2.)     
                    
        if type_spread == "ls":
            if k <= n: 
                k = k + int((n-k)/2.)
            else: k = k - int((k-n)/2.)
                    
        return k
    
    h = 0.704
    Lbox = 500
    omega_m = 0.272
    omega_t = 1
    
    def r200 (logM):
        t = (10**logM)**(1./3)  * (omega_m/omega_t)**(-1./3)
        return  (1.63e-2) * t * 1e3 # pc/h
        
    def conc (log_Mh,z):
        c1 = 2e13
        t1 = (10.**(log_Mh) /c1)**(-0.13)
        return (11./(1+z)) * t1

    def delta_c (c):
        return (200./3)* (c**3) /(np.log(1+c) -c/(1+c))

    def g(x):
        return np.log(1+float(x)) - x/(1+float(x))   
    
    # ========================================= 2 HOD ========================================= 
    ext_m = '_2HOD'
    NBIN = 60
    bins, HOD, Cent, Sats = np.loadtxt( hod_path + '/hod_d%i_60b.txt' %d, unpack = True )
    
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
        
    mask = (halomass_cen > bins[0]) & (halomass_cen < bins[-1])
    halomass_cen = halomass_cen[mask]
    x_cen, y_cen, z_cen = x_cen[mask], y_cen[mask], z_cen[mask]
    Len_Cat = len(x_cen)    
    
    #np.save('../Data/G13_millennium/Mock Catalogues/halomass_cen', halomass_cen) #TODO TO SAVE HALOMASSES 
    
    mass_index = ((( halomass_cen - bins[0] ) / ( bins[-1] - bins[0] )) * (NBIN-1)).astype('int32')
    n_gal_all = np.zeros(Len_Cat).astype('float64') # I DON'T NEED IT, ADAPT THE NGALS 2 CODE TODO
    n_gal_cen = np.zeros(Len_Cat).astype('float64')
    n_gal_sat = np.zeros(Len_Cat).astype('float64')

    SP.Ngals2(n_gal_all, n_gal_cen, n_gal_sat, halomass_cen, mass_index, m_all, n_all, m_cen, n_cen, m_sat, n_sat, Len_Cat)
    
    for i in xrange(Len_Cat):
        # CENTRAL
        n2 = n_gal_cen[i]
        dec2 = n2 % 1
        rr2 = np.random.rand()
        if rr2 <= dec2: n_gal_cen[i] = int(n2) + 1
        else: n_gal_cen[i] = int(n2)
        
        # SATTELITE
        n3 = n_gal_sat[i]
        rr3 = np.random.rand()
        if scatter == True:
            k3 = 0
            cdf3 = np.exp(-n3)
            
            while rr3 > cdf3:
                k3 += 1
                cdf3 = np.exp(-n3)*np.sum([(n3**j)/math.factorial(j) for j in range(0,k3+1)]) 
            
            if type_spread == "ms" or type_spread == "ls": k3 = KvaluE(k3, type_spread, n3)
            n_gal_sat[i] = int(k3)

        else:    
            dec3 = n3 % 1   
            if rr3 <= dec3: n_gal_sat[i] = int(n3) + 1
            else: n_gal_sat[i] = int(n3)  
            
    del n_gal_all # I DON'T NEED IT!
    
    
    # ===================================================== 1 HALO-TERM =======================================================
    
    # === ONLY APPLY FOR CONSTANT CONCENTRATION ===
    c = 13.981 #conc for a 10^12.5 Msun halo mass
    x0 = np.linspace(0,1,1000)
    y0 = np.zeros(1000)
    for l in range(len(x0)):
        r_norm0 = 1e-10
        cmf0 = g(c*r_norm0)/g(fvir*c)
        div = g(fvir*c)
        while x0[l] > cmf0:
            r_norm0 += 0.001
            cmf0 = g(c*r_norm0)/div
        y0[l] = r_norm0 
    y0[-1] = fvir
    spl = interp1d(x0, y0)
    # ==============================================
    
    zeros = np.where(n_gal_cen == 0)[0]
    NGAL = int(np.sum(n_gal_cen) + np.sum(n_gal_sat) + len(zeros))
    cat_data = np.zeros((NGAL,5), dtype = np.float32)
    
    #print NGAL
    #print np.sum(n_gal_cen)
    #print np.sum(n_gal_sat)
    
    k = 0
    for i in range(len(n_gal_cen)):
        if n_gal_cen[i] == 0:
            cat_data[k] = np.array([-1, -1, -1, halomass_cen[i], -1]) # HALOES THAT DO NOT ENTER TO THE SAMPLE
            k += 1
        
        else:
            cat_data[k] = np.array([x_cen[i], y_cen[i], z_cen[i], halomass_cen[i], 0])
            k += 1
        
        if n_gal_sat[i] != 0:
            NSAT = int(n_gal_sat[i])
            rr = np.random.random_sample(NSAT)
            r_norm = spl(rr)
            
            r = r_norm * r200(halomass_cen[i])/1e6
            phi = np.random.random_sample(NSAT) * 2*np.pi
            theta = np.random.random_sample(NSAT) * np.pi
            
            X = r*np.sin(theta)*np.cos(phi)
            Y = r*np.sin(theta)*np.sin(phi)
            Z = r*np.cos(theta)

            x_sat = x_cen[i] + X
            y_sat = y_cen[i] + Y
            z_sat = z_cen[i] + Z
            
            ones_array = np.ones(NSAT)
            halomass_array = np.ones(NSAT)*halomass_cen[i]

            fn = k + NSAT
            cat_data[k:fn,:] = np.array([x_sat, y_sat, z_sat, halomass_array, ones_array]).T
            k += NSAT
        
        #if k == NGAL: break
            
    
    # ==============================================================================================================================
    
    for i in range(NGAL):
        
        if cat_data[i][4] == -1:
            continue
        else:
            x,y,z = cat_data[i,:3]
            
            if x < 0: newx = x + Lbox; cat_data[i][0] = newx
            if x > Lbox: newx = x - Lbox; cat_data[i][0] = newx
            
            if y < 0: newy = y + Lbox; cat_data[i][1] = newy
            if y > Lbox: newy = y - Lbox; cat_data[i][1] = newy
            
            if z < 0: newz = z + Lbox; cat_data[i][2] = newz
            if z > Lbox: newz = z - Lbox; cat_data[i][2] = newz
    
    np.save('../Data/G13_millennium/Mock Catalogues/%s/2-HOD/cat_d%i%s%s_%irvir' %(sel, d, ext_scat, ext_sp, fvir), cat_data)   