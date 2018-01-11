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
    print "Method: 4-HOD"
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
    
    # =========================================== 4 HOD =================================================
    ext_m = '_4HOD'
    NBIN = 60
    bins, N_halos, N_cen, N_nocen, N_nocen_sat, N_cen_sat, N_cen_nosat, N_sat_cen, N_sat_nocen =  np.loadtxt( hod_path + '/hod_d%i_4HOD_60b.txt' %d, unpack = True)
    
    N_halos += 1e-10
    N_cen += 1e-10
    N_nocen += 1e-10
    N_nocen_sat += 1e-10
    N_cen_sat += 1e-10
    N_cen_nosat += 1e-10
    N_sat_cen += 1e-10
    N_sat_nocen += 1e-10
    
    Cent = np.log10( N_cen/N_halos)
    Cent_Sat = np.log10( N_sat_cen/N_cen_sat )
    Sat_noCent = np.log10( N_sat_nocen/N_nocen_sat )
    Cent_Sat2 = np.log10( N_cen_sat/N_cen )
    noCent_Sat = np.log10( N_nocen_sat/N_nocen)
    
    m_sat_cen, n_sat_cen = np.zeros(NBIN - 1).astype('float32'), np.zeros(NBIN - 1).astype('float32')
    m_sat_nocen, n_sat_nocen = np.zeros(NBIN - 1).astype('float32'), np.zeros(NBIN - 1).astype('float32')
    m_cen, n_cen = np.zeros(NBIN-1).astype('float32'), np.zeros(NBIN-1).astype('float32')
    m_cen_sat, n_cen_sat = np.zeros(NBIN-1).astype('float32'), np.zeros(NBIN-1).astype('float32')
    m_nocen_sat, n_nocen_sat = np.zeros(NBIN-1).astype('float32'), np.zeros(NBIN-1).astype('float32')
    
    
    for i in range( NBIN-1 ):
        y1_cen, y2_cen = Cent[i], Cent[i+1]
        y1_sat_cen, y2_sat_cen = Cent_Sat[i], Cent_Sat[i+1]
        y1_sat_nocen, y2_sat_nocen = Sat_noCent[i], Sat_noCent[i+1]
        y1_cen_sat, y2_cen_sat = Cent_Sat2[i], Cent_Sat2[i+1]
        y1_nocen_sat, y2_nocen_sat = noCent_Sat[i], noCent_Sat[i+1]
        
        x1, x2 = bins[i], bins[i+1]
        
        m_cen[i] = (y2_cen - y1_cen) / (x2-x1)
        n_cen[i] = y1_cen - (m_cen[i] * x1)
        
        m_cen_sat[i] = ( y2_cen_sat - y1_cen_sat ) / (x2 - x1)
        n_cen_sat[i] = y1_cen_sat - (m_cen_sat[i] * x1)
        
        m_nocen_sat[i] = ( y2_nocen_sat - y1_nocen_sat ) / (x2 - x1)
        n_nocen_sat[i] = y1_nocen_sat - (m_nocen_sat[i] * x1)
        
        m_sat_cen[i] = ( y2_sat_cen - y1_sat_cen ) / (x2 - x1)
        n_sat_cen[i] = y1_sat_cen - (m_sat_cen[i] * x1)
        
        m_sat_nocen[i] = ( y2_sat_nocen - y1_sat_nocen ) / (x2 - x1)
        n_sat_nocen[i] = y1_sat_nocen - (m_sat_nocen[i] * x1)

    mask = (halomass_cen > bins[0]) & (halomass_cen < bins[-1])
    halomass_cen = halomass_cen[mask]
    x_cen, y_cen, z_cen = x_cen[mask], y_cen[mask], z_cen[mask]
    Len_Cat = len(x_cen)
    
    mass_index = ((( halomass_cen - bins[0] ) / ( bins[-1] - bins[0] )) * (NBIN-1)).astype('int32')
    n_gal_cen = np.zeros(Len_Cat).astype('float64')
    n_gal_sat = np.zeros(Len_Cat).astype('float64')

    for i in range(Len_Cat):  
        index = mass_index[i] # Go through all haloes
        hmass = halomass_cen[i]
        cut_cen = 10**(m_cen[index]*hmass + n_cen[index])
        #cut_cen = N_cen[index] / N_halos[index] # tiene o no central?
        rr1 = np.random.rand()
        
        # CENTRAL
        if rr1 < cut_cen:
            n_gal_cen[i] = 1
            cut_cen_sat = 10**(m_cen_sat[index]*hmass + n_cen_sat[index])
            rr2 = np.random.rand()
            
            # SATELLITES
            if rr2 < cut_cen_sat:
                m = m_sat_cen[index]
                n = n_sat_cen[index] 
                Nsat = 10**((m * hmass) + n)    
                rr3 = np.random.rand()
            
                if scatter == True:
                    k3 = 0
                    cdf = np.exp(-Nsat)
                    while rr3 > cdf:
                        k3 += 1
                        cdf = np.exp(-Nsat)*np.sum([(Nsat**j)/math.factorial(j) for j in range(0,k3+1)]) 
                    
                    if type_spread == "ms" or type_spread == "ls": k3 = KvaluE(k3, type_spread, Nsat)
                    n_gal_sat[i] = int(k3)
                
                else:    
                    dec = Nsat % 1   
                    if rr3 <= dec: n_gal_sat[i] = int(Nsat) + 1
                    else: n_gal_sat[i] = int(Nsat)  
                                                
        else:
            cut_nocen_sat = 10**(m_nocen_sat[index]*hmass + n_nocen_sat[index])
            rr2 = np.random.rand()
            
            if rr2 < cut_nocen_sat:
                m = m_sat_nocen[index]
                n = n_sat_nocen[index] 
                Nsat = 10**((m * hmass) + n)
                rr3 = np.random.rand()
                
                if scatter == True:
                    k3 = 0
                    cdf = np.exp(-Nsat)
                    while rr3 > cdf:
                        k3 += 1
                        cdf = np.exp(-Nsat)*np.sum([(Nsat**j)/math.factorial(j) for j in range(0,k3+1)]) 
                    
                    if type_spread == "ms" or type_spread == "ls": k3 = KvaluE(k3, type_spread, Nsat)
                    n_gal_sat[i] = int(k3)
    
                else:    
                    dec = Nsat % 1   
                    if rr3 <= dec: n_gal_sat[i] = int(Nsat) + 1
                    else: n_gal_sat[i] = int(Nsat)           
                            


    # ===================================================== 1 HALO-TERM ==================================================

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
    
    #import time
    #start_time = time.time()
    #display_step = Len_Cat/100
    k = 0
    for i in range(Len_Cat):
        if n_gal_cen[i] == 0:
            cat_data[k] = np.array([-1, -1, -1, halomass_cen[i], -1]) # CENTRAL THAT DO NOT ENTER TO THE SAMPLE
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
            
        #if i % display_step == 0:
        #    print "Step: %i Time: %3fs" % (i, (time.time() - start_time))    

    # ===================================================================================================================

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

    np.save('../Data/G13_millennium/Mock Catalogues/%s/4-HOD/cat_d%i%s%s_%irvir' %(sel, d, ext_scat, ext_sp, fvir), cat_data)   
