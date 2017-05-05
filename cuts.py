import numpy as np
Lbox = 500.

def density_data (sfr_array, smass_array):

    densities_array = 10**(np.arange(-4,-1.4,.1))

    len_data = len(sfr_array)
    idx_sfr = np.argsort(sfr_array)
    idx_smass = np.argsort(smass_array)
    
    sort_sfr = sfr_array[idx_sfr]             # Sorted and reversed
    sort_smass = smass_array[idx_smass]
    
    sort_sfr = sort_sfr[::-1]
    sort_smass = sort_smass[::-1]             # Ordenado de Mayor a Menor
    
    data = np.zeros((len(densities_array),3))
    len_samples = np.zeros(len(densities_array))
    
    for i,d in enumerate(densities_array):
        
        len_sample = int(d * Lbox**3)
        #len_sample = int(d * len_data) 
        #len_samples[i] = len_sample
        
        min_sfr = sort_sfr[len_sample]
        min_mass = sort_smass[len_sample]
        #sfr_sample = sort_sfr[:len_sample]
        #smass_sample = sort_smass[:len_sample]
        
        #min_sfr = sfr_sample[-1]
        #min_mass = smass_sample[-1]

        data[i] = np.log10(d), min_mass, min_sfr
        
    np.savetxt('density_data.txt', data) 
    print " ---txt data has been saved---"
    
