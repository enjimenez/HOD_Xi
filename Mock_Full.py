def ToMakeMock (scatter= False, **inputs):
    import numpy as np
    import SPack as SP
    import math #TODO
    from scipy.interpolate import interp1d
    
    m = inputs['model']
    type_cut = inputs['selection']
    if type_cut == 'Mstell': type_cut = 'Stellar Mass' #TODO and SFR?
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
    print "Method: Full-HOD" 
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
    
    hod_path = path + '/hod_sams/' + type_cut
    cat_path = path + '/Mock Catalogues/' + type_cut
    
    if 'spread' in inputs: type_spread = inputs['spread']
    else: type_spread = -1
    def KvaluE (k, type_spread):
        if type_spread ==  'ms':
            if k <= n:
                k -= int((n-k)/2.)
                if k < 0: k = 0
            else: k += int((k-n)/2.)     
                    
        elif type_spread == 'ls':
            if k <= n: 
                k += int((n-k)/2.)
            else: k -= int((k-n)/2.)
                    
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
    
    # ========================================= FULL HOD ======================================
    NBIN = 60
    params = np.zeros((NBIN - 1, 2), dtype = np.float32)
    bins, HOD, Cent, Sats = np.loadtxt( hod_path + '/hod_d%i_60b.txt' %d, unpack = True )
    
    for i in range( NBIN-1 ):
        y1, y2 = HOD[i], HOD[i+1]
        x1, x2 = bins[i], bins[i+1]
        m = ( y2 - y1 ) / ( x2 - x1)
        n = y1 - (m * x1)
        params[i] = np.array([m, n])

    mask = (halomass_cen > bins[0]) & (halomass_cen < bins[-1])
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
        n3 = n_gal[i]
            
        if scatter == True:
            rr3 = np.random.random_sample()
            cdf3 = np.exp(-n3)
            k3 = 0 
            while rr3 > cdf3:
                k3 += 1
                cdf3 = np.exp(-n3)*np.sum([(n3**j)/math.factorial(j) for j in range(0,k3+1)]) 
            
            if type_spread == "ms" or type_spread == "ls": k3 = KvaluE(k3, type_spread, n3)
            nn_gal[i] = int(k3)
            
                
        else:    
            dec = n3 % 1
            rr = np.random.random_sample()
            if rr <= dec: nn_gal[i] = int(n3) + 1
            else: nn_gal[i] = int(n3)                                        
        
    del n_gal
    #cat_data = np.zeros((Len_Cat,5), dtype = np.float32)
    nzeros = np.where(nn_gal == 0)[0]
    Ntot = int(np.sum(nn_gal) + len(nzeros))
    cat_data = np.zeros((Ntot,5),dtype = np.float32)
    k = 0
    
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

    for i,n in enumerate(nn_gal):
        if n == 0: 
            cat_data[k] = np.array([-1, -1, -1, halomass_cen[i], -1]) # HALOES THAT DO NOT ENTER TO THE SAMPLE
            k += 1
            continue
        
        cat_data[k] = np.array([x_cen[i], y_cen[i], z_cen[i], halomass_cen[i], 0])
        k += 1
        if n == 1: continue
        
        # TODO Generate n-1 random numbers instead one (Use Shuffle code as example)
        j = 0
        while j < n-1:
            rr = np.random.rand()
            r_norm = spl(rr)
            r = r_norm * r200(halomass_cen[i])/1e6
            phi = np.random.rand() * 2*np.pi
            theta = np.random.rand() * np.pi

            x = r*np.sin(theta)*np.cos(phi)
            y = r*np.sin(theta)*np.sin(phi)
            z = r*np.cos(theta)
            
            x_sat = x_cen[i] + x
            if x_sat < 0: x_sat = x_sat + Lbox
            if x_sat > Lbox: x_sat = x_sat - Lbox

            y_sat = y_cen[i] + y
            if y_sat < 0: y_sat = y_sat + Lbox
            if y_sat > Lbox: y_sat = y_sat - Lbox
            
            z_sat = z_cen[i] + z
            if z_sat < 0: z_sat = z_sat + Lbox
            if z_sat > Lbox: z_sat = z_sat - Lbox 
        
            cat_data[k] = np.array([x_sat, y_sat, z_sat, halomass_cen[i], 1])
            k += 1
            j += 1
    
    np.save(cat_path + '/Full-HOD/cat_d%i%s%s_%irvir' %(d, ext_scat, ext_sp, fvir), cat_data)  
    # ext_m = 2HOD/_4HOD
    # ext_scat = poisson
    # ext_sp = ms/ls
    
    