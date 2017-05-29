# TODO ADD THE PREDICTIONS FOR CENTRALS AND SATS

def HOD (NBIN, type_cuts, number_densities, model, SAM = True, Mock = False, Shuffle = False ,**kwargs):
    import numpy as np
    data_path = "/home/esteban/Escritorio/Practica/Data/" + model
    
    for key in kwargs:
        if key == 'log_halomass': log_centralmvir = kwargs[key]
        if key == 'type_gal': type_array = kwargs[key]
        
    if SAM or Shuffle:
        logm_min = 10.
        logm_max = 15.
        Bw = (logm_max - logm_min)/NBIN
        mass_labels = ((log_centralmvir - logm_min)/(logm_max - logm_min) * NBIN).astype(int)
        bins = np.array([logm_min + Bw * i for i in range(NBIN)]) + Bw * 0.5
        densities_file = np.loadtxt(data_path + '/density_data.txt')
       
        if Shuffle:
            mask = (mass_labels >= 0) & (mass_labels < NBIN)
            type_array = type_array[mask]
            mass_labels = mass_labels[mask]
        del log_centralmvir
        
        if type(type_cuts) != list: type_cuts = [type_cuts]
        if type(number_densities) != list: number_densities = [number_densities]
        
        for label in type_cuts:
            
            if label == 'Mstell':                                               # ver bien donde poner esto
                if Shuffle: stellarmass = np.load(data_path + '/stellarmass_shuffle.npy')
                else: stellarmass = np.load(data_path + '/stellarmass.npy')
                tc = 'Stellar Mass/'
            if label == 'SFR':
                if Shuffle: sfr = np.load(data_path + '/sfr_shuffle.npy')
                else: sfr = np.load(data_path + '/sfr.npy')
                tc = 'SFR/'
            
            for d in number_densities:
                
                pos = int(d * 5) # position of the density cut in the density_data file 
                log_rho = densities_file[pos][0]
                if label == 'Mstell': cut = densities_file[pos][1]
                if label == 'SFR': cut = densities_file[pos][2]

                # To make HODs
                HOD  = np.zeros(NBIN)
                Cent = np.zeros(NBIN)
                Sats = np.zeros(NBIN)

                for i in range(NBIN):

                    idx = np.where(mass_labels == i)[0]
                    Halos = np.where(type_array[idx] == 0)[0] 

                    if label == 'Mstell': mask_HOD = stellarmass[idx] > cut
                    if label == 'SFR': mask_HOD = sfr[idx] > cut
                    HOD_ID = np.where(mask_HOD)[0]

                    mask_Cent = np.where(type_array[idx][HOD_ID] == 0)[0]
                    mask_Sats = np.where(type_array[idx][HOD_ID] != 0)[0]

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

                if Shuffle:
                    hod_path = data_path + '/hod_shuffle/' + tc 
                else: hod_path = data_path + '/hod_sams/' + tc 
                np.savetxt(hod_path + 'hod_d%i_%ib' %(d,NBIN), np.array([bins, HOD, Cent, Sats]).T)
                
    if Mock:
        logm_min = 10
        logm_max = 14.8333333
        NBIN = 29
        Bw = (logm_max - logm_min)/NBIN
        bins = np.array([logm_min + Bw*i for i in range(NBIN)]) + Bw * 0.5
        
        if type(type_cuts) != list: type_cuts = [type_cuts]
        if type(number_densities) != list: number_densities = [number_densities]
        
        for label in type_cuts:
            
            if label == 'Mstell': tc = 'Stellar Mass/'
            if label == 'SFR': tc = 'SFR/'
        
            cat_path = data_path + '/Mock Catalogues/'+ tc
    
            for d in number_densities:

                data = np.load(cat_path + 'cat_d%i_.npy' %d)
                halomass, n_all, n_cen, n_sat = data[:,3], data[:,4], data[:,5], data[:,6]

                mass_labels = ((halomass - logm_min)/(logm_max - logm_min) * NBIN).astype(int)
                
                HOD = np.zeros(NBIN)
                Cent = np.zeros(NBIN)
                Sats = np.zeros(NBIN)
                
                for i in range(NBIN):
                    
                    mass_mask = mass_labels == i
                    idx1 = np.where(mass_mask)[0]  # numero de halos
                    idx2 = np.where(mass_mask & (n_all != 0))[0]
                    idx3 = np.where(mass_mask & (n_cen != 0))[0]
                    idx4 = np.where(mass_mask & (n_sat != 0))[0]
                     
                    N_mask_all = n_all[idx2]
                    N_mask_cen = n_cen[idx3]
                    N_mask_sat = n_sat[idx4]
                    
                    Len_Haloes = len(idx1) + 1e-10
                    Len_Full = np.sum(N_mask_all)
                    Len_Cent = np.sum(N_mask_cen) 
                    Len_Sats = np.sum(N_mask_sat) 
                    HOD[i]   = Len_Full / Len_Haloes 
                    Cent[i] =  Len_Cent / Len_Haloes 
                    Sats[i] =  Len_Sats / Len_Haloes 

                HOD = np.log10(HOD + 1e-10)
                Cent = np.log10(Cent + 1e-10)
                Sats = np.log10(Sats + 1e-10)
                
                hod_path = data_path + '/hod_mocks/' + tc 
                np.savetxt(hod_path + 'hod_d%i_%ib_02' %(d,NBIN), np.array([bins, HOD, Cent, Sats]).T)
            