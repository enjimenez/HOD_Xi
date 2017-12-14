def HOD_Shuffle (**inputs):
    
    import numpy as np
    
    data_path = '../Data/G13_millennium'
    
    d = inputs['density_cut']
    selection = inputs['selection']
    fvir = int(inputs['rmax'])
    
    if selection == 'Mstell': stellarmass_shuffle = np.load(data_path + '/Shuffle/stellarmass_shuffle.npy'); sel = 'Stellar Mass'
    if selection == 'SFR': sfr_shuffle = np.load(data_path + '/Shuffle/sfr_shuffle.npy'); sel = 'SFR'
    
    bins = np.loadtxt(data_path + '/X_axis_60b.txt')
    Shuff_cat = np.load(data_path + '/Shuffle/%s/cat_shuffle_nfw_%irvir.npy' %(sel, fvir)) # THE CATALOGUE!
    
    pos = int(d * 5)
    log_rho, mass_cuts, sfr_cuts = np.loadtxt('../Data/G13_millennium/density_data.txt', unpack = True)
    rho, mass_cut, sfr_cut = log_rho[pos], mass_cuts[pos], sfr_cuts[pos]

    log_centralmvir = Shuff_cat[:,3]
    type_gal = Shuff_cat[:,4]
    del Shuff_cat

    logm_min = 10.
    logm_max = 15.
    NBIN = 60
    mass_labels = ((log_centralmvir - logm_min)/(logm_max - logm_min) * NBIN).astype(int)

    HOD  = np.zeros(NBIN)
    Cent = np.zeros(NBIN)
    Sats = np.zeros(NBIN)
    
    i = 0
    for i in range(NBIN):
        
        idx = np.where(mass_labels == i)[0]
        Halos = np.where(type_gal[idx] == 0)[0] 

        if selection == 'Mstell': mask_HOD = stellarmass_shuffle[idx] > mass_cut
        if selection == 'SFR': mask_HOD = sfr_shuffle[idx] > sfr_cut
        HOD_ID = np.where(mask_HOD)[0]
        
        mask_gals = type_gal[idx][HOD_ID]
        mask_Cent = np.where(mask_gals == 0)[0]
        mask_Sats = np.where(mask_gals == 1)[0]

        Len_Full = len(HOD_ID) 
        Len_Cent = len(mask_Cent) 
        Len_Sats = len(mask_Sats) 
        Len_Haloes = float(len(Halos)) + 1e-10
    
        HOD[i]      = Len_Full / Len_Haloes 
        Cent[i]     = Len_Cent / Len_Haloes
        Sats[i]     = Len_Sats / Len_Haloes
        
    HOD = np.log10(HOD + 1e-10)
    Cent = np.log10(Cent + 1e-10)
    Sats = np.log10(Sats + 1e-10)
    
    np.savetxt(data_path + '/hod_shuffle/%s/hod_d%i_60b_%irvir.txt' %(sel, d, fvir), np.array([bins, HOD, Cent, Sats]).T)