def HOD_std (densitycut, typecut):
    import numpy as np
    import time
    import matplotlib.pyplot as plt
    import os
    from matplotlib import rc
    rc('text', usetex=True)
    rc('font', family='serif')
    
    if typecut == 'Mstell': tp = 'Stellar Mass'
    if typecut == 'SFR': tp = 'SFR'
    
    rhos = np.array([-3.5, -3.0, -2.5, -2.0, -1.5])
    
    path = "/home/esteban/Escritorio/Esteban HOD/Data/G13_millennium" # Guo2013 
    #path = "/home/esteban/Escritorio/Practica/Data/G14_millennium" # Gonzalez2014 
    #path = "/home/esteban/Escritorio/Practica/Data/L16_millennium" # Lacey2016 

    # Loading database
    #x = np.load(path +'/x.npy')
    fofid = np.load(path +'/fofid.npy')
    if typecut == 'Mstell': stellarmass =  np.load(path +'/stellarmass.npy')
    log_centralmvir  = np.log10(np.load(path +'/centralmvir.npy') + 1e-10) + 10
    type_gal = np.load(path +'/type.npy')
    if typecut == 'SFR': sfr = np.load(path + '/sfr.npy')
    
    #log_centralmvir = centralmvir
    #del centralmvir
    logm_min, logm_max = 10, 15.
    NBIN = 60

    labels = ((log_centralmvir - logm_min)/(logm_max - logm_min) * NBIN).astype('int')
    del log_centralmvir

    mask = (labels >= NBIN) | (labels < 0)
    fofid = fofid[~mask]
    if typecut == 'Mstell': stellarmass = stellarmass[~mask]
    #log_centralmvir = log_centralmvir[~mask]
    type_gal = type_gal[~mask]
    if typecut == 'SFR': sfr = sfr[~mask]
    labels = labels[~mask]
    #x = x[~mask]

    bin_width = (logm_max - logm_min)/NBIN
    bins = np.array([logm_min + bin_width*i for i in range(NBIN)]) + bin_width *.5
    
    # Density selection
    # =================
    #   Fila    log(rho)
    #    25      -1.5
    #    20      -2.0
    #    15      -2.5
    #    10      -3.0
    #     5      -3.5
    # =================
    pos  = int(densitycut*5)
    log_rho, mass_cuts, sfr_cuts = np.loadtxt(path + '/density_data.txt', unpack = True)
    rho, mass_cut, sfr_cut = log_rho[pos], mass_cuts[pos], sfr_cuts[pos]

    # Centrals data
    idx = np.where(type_gal == 0)[0]
    fofid_ar = fofid[idx]
    #logcentralmvir_ar = log_centralmvir[idx]
    labels_ar = labels[idx]
    del idx
    #x_ar = x[idx]

    # Number of central = number of halos
    len_data = len(fofid_ar)

    # To save halos data
    All_ar = np.zeros(len_data)
    Cent_ar = np.zeros(len_data)
    Sats_ar = np.zeros(len_data)

    i,j = 0,0
    Nsats, Ncent, Nall = 0,0,0
    len_fulldata = len(fofid)

    # To make HODs
    Len_Sats = np.zeros(NBIN)
    Len_Cent = np.zeros(NBIN)
    Len_All = np.zeros(NBIN)
    Halos = np.zeros(NBIN)

    display_step = len_data/10
    
    if typecut == 'Mstell':
        while i < len_data:
            index = labels_ar[i]
            cond = (fofid[j] == fofid_ar[i])
            while cond:
                if stellarmass[j] > mass_cut: 
                    Len_All[index] += 1
                    Nall += 1
                    if type_gal[j] != 0: 
                        Len_Sats[index] += 1
                        Nsats += 1
                    else:
                        Len_Cent[index] += 1 
                        Ncent += 1
                j += 1
                if j == len_fulldata: break
                cond = (fofid[j] == fofid_ar[i])
            
            Halos[index] += 1
            # Atencion!! tambien se guardan halos sin galaxias (no cumplen el criterio)
            Cent_ar[i] = Ncent
            Sats_ar[i] = Nsats
            All_ar[i] = Nall
            Ncent, Nsats, Nall = 0,0,0
            i += 1
        
    print "termine"
    #return All_ar, Cent_ar, Sats_ar, Len_All, Len_Cent, Len_Sats, labels_ar

    print "me pase a otro lado"    
        
    if typecut == 'SFR':
        while i < len_data:
            index = labels_ar[i]
            cond = (fofid[j] == fofid_ar[i])
            while cond:
                if sfr[j] > sfr_cut: 
                    Len_All[index] += 1
                    Nall += 1
                    if type_gal[j] != 0: 
                        Len_Sats[index] += 1
                        Nsats += 1
                    else:
                        Len_Cent[index] += 1 
                        Ncent += 1
                j += 1
                if j == len_fulldata: break
                cond = (fofid[j] == fofid_ar[i])
            
            Halos[index] += 1
            # Atencion!! tambien se guardan halos sin galaxias (no cumplen el criterio)
            Cent_ar[i] = Ncent
            Sats_ar[i] = Nsats
            All_ar[i] = Nall
            Ncent, Nsats, Nall = 0,0,0
            i += 1
    
    HOD, Cent, Sats = np.zeros(NBIN), np.zeros(NBIN), np.zeros(NBIN)
    std_hod, std_cen, std_sat = np.zeros(NBIN), np.zeros(NBIN), np.zeros(NBIN)
    alpha_gal = np.zeros(NBIN)
    alpha_sat = np.zeros(NBIN)
    
    del fofid, fofid_ar, labels
    Nhalos = np.zeros(NBIN)
    
    for i in range(NBIN):
        idx = np.where(labels_ar == i)[0]
        Len_Haloes = len(idx)
        
        Ngal_halo = All_ar[idx]
        Ncen_halo = Cent_ar[idx]
        Nsat_halo = Sats_ar[idx]

        Len_Full = np.sum(Ngal_halo)
        Len_Cent = np.sum(Ncen_halo)
        Len_Sats = np.sum(Nsat_halo)

        HOD[i]      = Len_Full / Len_Haloes 
        Cent[i]     = Len_Cent / Len_Haloes
        Sats[i]     = Len_Sats / Len_Haloes
        all_pairs_halo = np.array([Ngal_halo[k]*(Ngal_halo[k]-1) for k in range(Len_Haloes)])
        sat_pairs_halo = np.array([Nsat_halo[k]*(Nsat_halo[k]-1) for k in range(Len_Haloes)])
        
        alpha_gal[i] = np.sqrt(np.mean(all_pairs_halo + 1e-10))/(HOD[i] + 1e-10)
        alpha_sat[i] = np.sqrt(np.mean(sat_pairs_halo + 1e-10))/(Sats[i] + 1e-10)
        
        #sum_all = np.sum(np.array([(Ngal_halo[k] - HOD[i])**2 for k in range(Len_Haloes)]))
        #sum_cen = np.sum(np.array([(Ncen_halo[k] - Cent[i])**2 for k in range(Len_Haloes)]))
        #sum_sat = np.sum(np.array([(Nsat_halo[k] - Sats[i])**2 for k in range(Len_Haloes)]))
        
        # STANDARD SCATTER
        #std_hod[i] = np.sqrt(sum_all/(Len_Haloes -1))
        #std_cen[i] = np.sqrt(sum_cen/(Len_Haloes -1))
        #std_sat[i] = np.sqrt(sum_sat/(Len_Haloes -1))
        Nhalos[i] = Len_Haloes
         
    #HOD += 1e-10
    #Cent += 1e-10
    #Sats += 1e-10
    #log_HOD = np.log10(HOD)
    #log_Cent = np.log10(Cent)
    #log_Sats = np.log10(Sats)
    
    #poisson = np.sqrt(Nhalos) 
    #poisson = np.sqrt(Sats)
    #log_err_max =  np.log10(Sats + poisson)
    #log_err_min =  np.log10(Sats - poisson)
    
    #cent_std = np.sqrt(Cent)
    #log_err1 =  np.log10(Cent + cent_std)
    #log_err2 =  np.log10(Cent - cent_std)

    
    plot_path = "/home/esteban/Escritorio/Esteban HOD/Plots/G13_millennium/hod_sams/%s" %tp
    f = plt.figure(figsize =(10,7))
    plt.plot(bins, alpha_gal, 'k.', label = r"$\rm All$")
    plt.plot(bins, alpha_sat, 'b.',label = r"$\rm Sats$")
    plt.ylabel(r"$\left<N(N-1)\right>^{1/2}/ \left<N \right>$", fontsize = 24)
    #plt.plot(bins, log_HOD, 'k-', lw = 1, label = r"$\rm All$") # original lw = .8
    #plt.plot(bins, log_Cent, 'r-', lw = 1, label = r"$\rm Cent$")
    #plt.plot(bins, log_Sats, 'b-', lw = 1, label = r"$\rm Sats$")
    plt.axhline(y=1, color = 'black', ls = '--', lw = '.8')
    #plt.plot([], [], 'b--', lw = 2, label = r"\sqrt(\sigma^2)" )
    #plt.plot([], [], 'b:', lw = 2, label = r"$\sqrt(N)$" )
    #plt.fill_between(bins, log_HOD-log_HOD_std, log_HOD + log_HOD_std, facecolor='black', alpha=0.1)
    #plt.fill_between(bins, log_Cent-log_Cent_std, log_Cent + log_Cent_std, facecolor='red', alpha=0.1)
    #plt.fill_between(bins, log_Sats-log_Sats_std, log_Sats + log_Sats_std, facecolor='blue', alpha=0.2)
    
    #plt.plot(bins, np.log10(Cent - std_cen + 1e-10), 'r--', lw=0.5)
    #plt.plot(bins, np.log10(Cent + std_cen + 1e-10), 'r--', lw=0.5)

    #plt.plot(bins, np.log10(Sats - std_sat + 1e-10), 'b--', lw=0.5)
    #plt.plot(bins, np.log10(Sats + std_sat + 1e-10), 'b--', lw=0.5)
    
    
    #for i in range(NBIN):
    #    print log_Sats[i], log_err_max[i], log_err_min[i]
    # 
    #a = log_Sats
   # b = log_err_max
    #c = log_err_min
    
    #return bins, a,b,c
    
    #for i in range(NBIN):
    #    print log_Sats[i], log_err_max[i], log_err_min[i]
    
    
    #plt.plot(bins, log_err_max, 'b:', lw=0.5)
    #plt.plot(bins, log_err_min, 'b:', lw=0.5)

    plt.title(r"$n = 10^{%.1f} /\rm h^{-3}\ Mpc^3$" %rhos[densitycut-1], fontsize = 24 ) #TODO
    plt.xlabel(r"$\log(\rm M_h / h^{-1}M_{\odot})$", fontsize = 24)
    #plt.ylabel(r"$\log(\rm N)$", fontsize = 24)
    #plt.ylim(-2,3)
    plt.ylim(0,1.2)
    plt.xlim(10,15)
    plt.xticks(fontsize = 24)
    plt.yticks(fontsize = 24)
    plt.minorticks_on()
    plt.legend(loc='lower right',title = r"\Large{\textbf{G13}}", frameon = False, prop={'size':20}) #TODO
    plt.show()
    f.savefig(plot_path + '/alpha.pdf')