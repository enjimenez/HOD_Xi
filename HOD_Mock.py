def HOD_Mock (**inputs):
    import numpy as np
    
    data_path = '../Data/G13_millennium'
    
    d = inputs['density_cut']
    selection = inputs['selection']
    fvir = int(inputs['rmax'])
    
    if selection == 'Mstell': sel = 'Stellar Mass'
    if selection == 'SFR': sel = 'SFR'
    
    bins = np.loadtxt('../Data/G13_millennium/X_axis_60b.txt')
    Mock = np.load('../Data/G13_millennium/Mock Catalogues/%s/cat_d%i_nfw_%irvir.npy' %(sel,d,fvir)) #THE CATALOGUE!
    log_centralmvir = Mock[:,3]
    type_gal = Mock[:,4]
    del Mock

    logm_min = 10.
    logm_max = 15
    NBIN = 60
    mass_labels = ((log_centralmvir - logm_min)/(logm_max - logm_min) * NBIN).astype(int)

    HOD  = np.zeros(NBIN)
    Cent = np.zeros(NBIN)
    Sats = np.zeros(NBIN)

    summ = 0
    for i in range(NBIN):
        
        mass_mask = mass_labels == i
        #idx = np.where(mass_mask)[0]  # numero de halos

        idx = np.where(mass_labels == i)[0]
        mask_type = type_gal[idx]
        idx2 = np.where(type_gal[idx] !=-1)[0]
        Halos = np.where((type_gal[idx] == 0) | (type_gal[idx] == -1))[0] 

        mask_Cent = np.where(mask_type == 0)[0]
        mask_Sats = np.where(mask_type == 1)[0]

        Len_Full = len(idx2) 
        Len_Cent = len(mask_Cent) 
        Len_Sats = len(mask_Sats) 
        Len_Haloes = float(len(Halos)) + 1e-10

        HOD[i]      = Len_Full / Len_Haloes 
        Cent[i]     = Len_Cent / Len_Haloes
        Sats[i]     = Len_Sats / Len_Haloes

    HOD = np.log10(HOD + 1e-10)
    Cent = np.log10(Cent + 1e-10)
    Sats = np.log10(Sats + 1e-10)
    
    np.savetxt(data_path + '/hod_mocks/%s/hod_d%i_60b_%irvir.txt' %(sel, d, fvir), np.array([bins, HOD, Cent, Sats]).T)