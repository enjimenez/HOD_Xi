def HOD_SAM (**inputs):
    import numpy as np

    data_path = '../Data/G13_millennium'
    
    d = inputs['density_cut']
    selection = inputs['selection']

    log_centralmvir = np.log10(np.load('../Data/G13_millennium/centralmvir.npy') + 1e-10) + 10
    type_gal = np.load('../Data/G13_millennium/type.npy')

    logm_min = 10.
    logm_max = 15.
    NBIN = 60
    Bw = (logm_max - logm_min)/NBIN
    mass_labels = ((log_centralmvir - logm_min)/(logm_max - logm_min) * NBIN).astype(int)
    bins = np.loadtxt(data_path + '/X_axis_60b.txt')  #TODO
    densities_file = np.loadtxt(data_path + '/density_data.txt')

    del log_centralmvir

    if selection == 'Mstell': stellarmass = np.load(data_path + '/stellarmass.npy'); sel = 'Stellar Mass'

    if selection == 'SFR': sfr = np.load(data_path + '/sfr.npy'); sel = 'SFR'

    pos = int(d * 5) # position of the density cut in the density_data file 
    log_rho = densities_file[pos][0]
    if selection == 'Mstell': cut = densities_file[pos][1]
    if selection == 'SFR': cut = densities_file[pos][2]

    # To make HODs
    HOD  = np.zeros(NBIN)
    Cent = np.zeros(NBIN)
    Sats = np.zeros(NBIN)

    for i in range(NBIN):

        idx = np.where(mass_labels == i)[0]
        Halos = np.where(type_gal[idx] == 0)[0] 

        if selection == 'Mstell': mask_HOD = stellarmass[idx] > cut
        if selection == 'SFR': mask_HOD = sfr[idx] > cut
        HOD_ID = np.where(mask_HOD)[0]

        mask_Cent = np.where(type_gal[idx][HOD_ID] == 0)[0]
        mask_Sats = np.where(type_gal[idx][HOD_ID] != 0)[0]

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

    np.savetxt(data_path + '/hod_sams/%s/hod_d%i_60b.txt' %(sel, d), np.array([bins, HOD, Cent, Sats]).T)