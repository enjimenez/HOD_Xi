def HOD (model, density, type_cut, Mock = False, cSAM = False, cShuffle = False):
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib import rc
    rc('text', usetex=True)
    rc('font', family='serif')
    
    ab = model
    model = ab + '_millennium'
    rhos = np.array([-3.5, -3.0, -2.5, -2.0, -1.5])
    
    if type_cut == 'Mstell': tp = 'Stellar Mass'
    if type_cut == 'SFR': tp = 'SFR'
    
    if Mock:
        
        data_path_mock = '../Data/' + model + '/hod_mocks/' + tp + '/hod_d%i_29b' % density
        data_path_sams = '../Data/' + model + '/hod_sams/' + tp + '/hod_d%i_30b' % density    
        data_path_shuffle = '../Data/' + model + '/hod_shuffle/' + tp + '/hod_d%i_30b' % density    
        
        plot_path = '../Plots/' + model + '/HOD/' + tp  # TODO /hod_cat
        
        f = plt.figure()
        bins, HOD = np.loadtxt(data_path_mock, unpack = True)
        plt.plot(bins, HOD, 'k--', markersize = 1.5, label = r"All / Mock")
        
        if cSAM or cShuffle:
            if cSAM:
                bins_sam, hod_sam = np.loadtxt(data_path_sams, usecols = (0,1), unpack = True)
                plt.plot(bins_sam, hod_sam, 'b--', markersize = 1.5, label = r"All / SAM")
            if cShuffle:
                bins_shuffle, hod_shuffle = np.loadtxt(data_path_shuffle, usecols = (0,1), unpack = True)
                plt.plot(bins_shuffle, hod_shuffle, 'r--', markersize = 1.5, label = r"All / Shuffle")
                
        plt.title(r"$n = 10^{%.1f} /h^{-3} Mpc^3$" %rhos[density-1])
        plt.xlabel(r"$log(M_h / h^{-1}M_{\odot})$", fontsize = 10)
        plt.ylabel(r"$log(N)$", fontsize = 10)
        plt.ylim(-2,3)
        plt.xlim(10,15)
        plt.minorticks_on()
        plt.legend(loc='upper left',title = r"\textbf{%s}" %ab, frameon = False, prop={'size':10})
        f.savefig(plot_path + '/hod_d%i.pdf' %density)
        
    else:
        print "Soon... :)"  #TODO
        
def xi (model, density, type_cut, Mock = False, Shuffle = False, Ratio = False, **kwargs):
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib import rc
    from matplotlib import gridspec
    rc('text', usetex=True)
    rc('font', family='serif')
    
    ab = model
    model = ab + '_millennium'
    rhos = np.array([-3.5, -3.0, -2.5, -2.0, -1.5])   #TODO Allow density array be unsorted
    
    if (Shuffle) or (Mock):
        
        if Shuffle and Mock:
            if type(type_cut) != list: type_cut = [type_cut]
            if type(density) != list: density = [density]
        
            colors =  ['g-', 'r-', 'b-', 'k-']
            
            if len(kwargs) > 1:
                for key in kwargs:
                    if key == xlim: xlim = kwargs[key]
                    if key == ry_range: ry_range = kwargs[key]
            else:
                print "entre al else"
                xlim = [0.3,1.7]
                ry_range = [0.8,1.2]
                
            for label in type_cut:
                
                if label == 'Mstell': tp = 'Stellar Mass'
                if label == 'SFR': tp = 'SFR'                    
                plot_path = '../Plots/' + model + '/xi_sam/Assembly Bias/'+ tp
                
                for d in density:
                    
                    #f = plt.figure()  
                    
                    data_path_sam_shuffle = '../Data/' + model + '/xi_sam/' + '/Shuffle/' + tp + '/xi_d%i.txt' % d
                    data_path_mock = '../Data/' + model + '/xi_mocks/' + tp + '/xi_cat_0%i.txt' % d
                    data_path_sam = '../Data/' + model + '/xi_sam/' + tp + '/xi_d%i_fixed.txt' % d
                    
                    bins_shuffle, xi_shuffle = np.loadtxt(data_path_sam_shuffle, usecols = (0,1), unpack = True)
                    bins_sam,     xi_sam = np.loadtxt(data_path_sam, usecols = (0,1), unpack = True)
                    bins_mock,    xi_mock = np.loadtxt(data_path_mock, usecols = (0,1), unpack = True)
                    
                    logxi_sam, logxi_mock, logxi_shuff = np.log10(xi_sam), np.log10(xi_mock), np.log10(xi_shuffle)
                    
                    mask = (bins_sam >= xlim[0]) & (bins_sam < xlim[1])
                    xi_sam_mask = xi_sam[mask]
                    xi_mock_mask = xi_mock[mask]
                    xi_shuffle_mask = xi_shuffle[mask]
                    
                    ratio1 = xi_sam_mask / xi_mock_mask
                    ratio2 = xi_sam_mask / xi_shuffle_mask
                
                    fig = plt.figure(figsize=(8, 6)) 
                    gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1]) 
                    gs.update(left = 0.1, right = 0.98, hspace = 0)
                    ax0 = plt.subplot(gs[0])
                    ax0.set_xticklabels([])
                    yticks = ax0.yaxis.get_major_ticks()
                    yticks[0].label1.set_visible(False)
                    ax0.plot(bins_sam, logxi_sam, 'k-', label = r"$\xi$ from SAM")
                    ax0.plot(bins_mock, logxi_mock, 'k--', label = r"$\xi$ from Mock")
                    ax0.plot(bins_mock, logxi_shuff, 'k:', label = r"$\xi$ from Shuffle")
                    ax0.set_title(r"$n = 10^{%.1f} /h^{-3} Mpc^3$" %rhos[d-1], fontsize = 12)
                    ax0.set_ylabel(r"$log(\xi)$", fontsize = 12)
                    ax0.legend(loc='best', frameon=False)
                    ax0.set_xlim(-0.5,2)
                    ax0.set_ylim(-2,3)

                    ax1 = plt.subplot(gs[1])
                    ax1.axhline(y=1, color='black', linestyle='--', lw=.4)
                    ax1.plot(bins_sam[mask], ratio1, 'r-', label = r"$\xi_{SAM} / \xi_{Mock}$")
                    ax1.plot(bins_shuffle[mask], ratio2, 'b-', label = r"$\xi_{SAM} / \xi_{Shuff}$")
                    #ax1.set_ylabel(r"$\xi_{SAM} / \xi_{Mock}$,$\xi_{SAM} / \xi_{Shuff}$")
                    ax1.set_xlim(-0.5,2)
                    ax1.set_ylim(ry_range[0],ry_range[1])
                    ax1.legend(loc='best', frameon=False, prop={'size':8})
                    ax1.set_xlabel("$log(d/h^{-1}Mpc)$", fontsize = 12)
                    
                    fig.savefig(plot_path + '/xi_cat_d%i.pdf' % d)
                    
                    #ratio = xi_org/xi_shuffle
                    #plt.plot(bins, ratio, 'k-')
                    #plt.axhline(y=1, color = 'black', linestyle = '--', lw = .4)
                    #plt.title(r"$n = 10^{%.1f} /h^{-3} Mpc^3$" %rhos[density-1])
                    ##plt.annotate(r"\large{\textbf{%s}}" %ab, xy=(0.06,0.1), xycoords='axes fraction' )
                    #plt.xlabel("$log(d/h^{-1}Mpc)$", fontsize = 14)
                    #plt.ylabel(r"$b = \left(\xi_{org} / \xi_{shuff}\right)^{1/2}$", fontsize = 14)
                    ##plt.ylim(-1,6)
                    #plt.xlim(-2,1.5)
                    #plt.minorticks_on
                    #f.savefig(plot_path + '/xi_d%i.pdf' %d)
        
        else:
            if len(kwargs) > 1:
                for key in kwargs:
                    if key == xlim: xlim = kwargs[key]
                    if key == ry_range: ry_range = kwargs[key]
            else:
                xlim = [0.3,1.7]
                ry_range = [0.8,1.2]
            
            if type(type_cut) != list: type_cut = [type_cut]
            if type(density) != list: density = [density]
            
            for label in type_cut:
                if label == 'Mstell': tp = 'Stellar Mass'
                if label == 'SFR': tp = 'SFR'                    
                plot_path = '../Plots/' + model + '/xi_cat/' + tp
                
                for d in density:
                    data_path_mock1 = '../Data/' + model + '/xi_mocks/' + tp + '/xi_cat_0%i.txt' % d  #TODO CAMBIA CON EL METODO
                    data_path_mock2 = '../Data/' + model + '/xi_mocks/' + tp + '/xi_cat_0%i_2HOD.txt' % d
                    data_path_sams = '../Data/' + model + '/xi_sam/' + tp + '/xi_d%i_fixed.txt' % d
                    plot_path = '../Plots/' + model + '/xi_cat/' + tp
                    
                    bins_mock1, xi_mock1 = np.loadtxt(data_path_mock1, usecols = (0,1), unpack = True)
                    logxi_mock1 = np.log10(xi_mock1)
                    
                    bins_mock2, xi_mock2 = np.loadtxt(data_path_mock2, usecols = (0,1), unpack = True)
                    logxi_mock2 = np.log10(xi_mock2)
                        
                    bins_sam, xi_sam = np.loadtxt(data_path_sams, usecols = (0,1), unpack = True)
                    logxi_sam = np.log10(xi_sam)
                        
                    if Ratio:
                        mask = (bins_sam >= xlim[0]) & (bins_sam < xlim[1])
                        xi_sam_mask = xi_sam[mask]
                        xi_mock_mask1 = xi_mock1[mask]
                        xi_mock_mask2 = xi_mock2[mask]
                        ratio1 = xi_sam_mask / xi_mock_mask1
                        ratio2 = xi_sam_mask / xi_mock_mask2
                    
                        fig = plt.figure(figsize=(8, 6)) 
                        gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1]) 
                        gs.update(left = 0.1, right = 0.98, hspace = 0)
                        ax0 = plt.subplot(gs[0])
                        ax0.set_xticklabels([])
                        yticks = ax0.yaxis.get_major_ticks()
                        yticks[0].label1.set_visible(False)
                        ax0.plot(bins_sam, logxi_sam, 'k-', label = r"$\xi$ from SAM")
                        ax0.plot(bins_mock1, logxi_mock1, 'r--', label = r"$\xi$ from Mock (1 HOD)")
                        ax0.plot(bins_mock2, logxi_mock2, 'b--', label = r"$\xi$ from Mock (2 HOD)")
                        ax0.set_title(r"$n = 10^{%.1f} /h^{-3} Mpc^3$" %rhos[d-1], fontsize = 12)
                        ax0.set_ylabel(r"$log(\xi)$", fontsize = 12)
                        ax0.legend(loc='best', frameon=False)
                        ax0.set_xlim(-0.5,2)
                        ax0.set_ylim(-2,3)

                        ax1 = plt.subplot(gs[1])
                        ax1.axhline(y=1, color='black', linestyle='--', lw=.4)
                        ax1.plot(bins_sam[mask], ratio1, 'r-')
                        ax1.plot(bins_sam[mask], ratio2, 'b-')
                        ax1.set_ylabel(r"$\xi_{SAM} / \xi_{Mock}$")
                        ax1.set_xlim(-0.5,2)
                        ax1.set_ylim(ry_range[0],ry_range[1])
                        ax1.set_xlabel("$log(d/h^{-1}Mpc)$", fontsize = 12)
                        
                        fig.savefig(plot_path + '/xi_cat_d%i.pdf' % d)   #TODO CAMBIA DE ACUERDO AL METODO DE POBLACION USADO
                        
                    else:
                        f = plt.figure()    
                        plt.plot(bins_mock, logxi_mock, 'k--', label = r"$\xi$ from Mock")
                        plt.plot(bins_sam, logxi_sam, 'k-', label = r"$\xi$ from SAM" )
                        plt.title(r"$n = 10^{%.1f} /h^{-3} Mpc^3$" %rhos[density-1], fontsize = 12)
                        plt.xlabel("$log(d/h^{-1}Mpc)$", fontsize = 12)
                        plt.ylabel(r"$log(\xi)$", fontsize = 12)
                        plt.ylim(-2,3)
                        plt.xlim(-2,2)
                        plt.legend(loc='best', frameon=False)
                        f.savefig(plot_path + '/xi_cat_d%i.pdf' % density)
        
    else:
        if type(type_cut) != list: type_cut = [type_cut]
        if type(density) != list: density = [density]
        
        colors =  ['g-', 'r-', 'b-', 'k-']
            
        for label in type_cut:
            
            if label == 'Mstell': tp = 'Stellar Mass'
            if label == 'SFR': tp = 'SFR'
            plot_path = '../Plots/' + model + '/xi_sam/' + tp
            f = plt.figure()  
            
            for d in density:
                
                data_path_sam = '../Data/' + model + '/xi_sam/' + tp + '/xi_d%i_fixed.txt' % d
                bins, xi = np.loadtxt(data_path_sam, usecols = (0,1), unpack = True)
                logxi = np.log10(xi)
                plt.plot(bins, logxi, colors[d-1], label = r"$n = 10^{%.1f} /h^{-3} Mpc^3$" %rhos[d-1])
            
            plt.legend(loc = 'best',title = r"\large{%s}" %tp, frameon = False, prop={'size':12})
            plt.annotate(r"\large{\textbf{%s}}" %ab, xy=(0.06,0.1), xycoords='axes fraction' )
            plt.xlabel("$log(d/h^{-1}Mpc)$", fontsize = 14)
            plt.ylabel(r"$log(\xi)$", fontsize = 14)
            plt.ylim(-1,6)
            plt.xlim(-2,1.5)
            plt.minorticks_on
            f.savefig(plot_path + '/xi_%s.pdf' % label)
            
            
                    
                    
                
                
                
                
            
      
        
        
            
            
            
     
     
        
