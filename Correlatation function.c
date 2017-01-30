# include <stdio.h>
# include <math.h>
# include <time.h>

# define ROWS 3269 // Number of galaxies  M > 10^9 Msun -> #14193, M > 10^10 Msun -> #3269
# define Lbox 62.5 // Size of milimillenium simulation
# define NBIN 30   // Number of bins

double logdist (float r1[], float r2[]){
    
    float Lx, Ly, Lz;
    double ret;
    
    Lx = fabs(r2[0] - r1[0]);
    Ly = fabs(r2[1] - r1[1]);    
    Lz = fabs(r2[2] - r1[2]);
    if (Lx >= Lbox/2){
        Lx = Lbox - Lx;
    }
    if (Ly >= Lbox/2){
        Ly = Lbox - Ly;
    }
    if (Lz >= Lbox/2){
        Lz = Lbox - Lz;
    }

    ret = .5 * log10(Lx*Lx + Ly*Ly + Lz*Lz);
    return ret;
}

int ToClass (float d,  float Xmin, float Xmax){
    
    int label;
    label = ((d-Xmin)/(Xmax - Xmin)) * NBIN;
    return label;
}

double RR (float r, int n, float dr){
    
    double V_c, V_t, ret;
    V_c = 4/3. * M_PI *(pow(pow(10,(r+dr)), 3) - pow(pow(10,r), 3));
    V_t = pow(62.5, 3);
    ret = (V_c/V_t) * pow(n,2);
    return ret;
}

int main(void){
    
    clock_t begin = clock();
    
    // Read Data
    FILE *data_file;
    
    long long IDs[ROWS], FOFID[ROWS], ID, fofid;
    double  X[ROWS], Y[ROWS], Z[ROWS], x, y, z;
    int i,j;
    
    data_file = fopen("DL07_milimill_10_xi.csv", "r");
    
    for( i = 0; i < ROWS; i++ ){
            
        fscanf(data_file, "%lli %lli %lf %lf %lf ", &ID, &fofid, &x, &y, &z);
        IDs[i] = ID, FOFID[i] = fofid, X[i] = x, Y[i] = y, Z[i] = z;
    
    }

    // Variables
    int label, DD_array[NBIN];           
    int DD_1h[NBIN], DD_2h[NBIN];                      // 1 and 2 halo terms
    float Dmin = -2, Dmax = 1;                         // Min and max distance in log10 [d] = Mpc
    float Binsize = (Dmax - Dmin)/NBIN;                // Width of bins  
    float r1[3], r2[3], bins[NBIN], dist_array[NBIN];  // 3D positions and bins arrays      
    long long fofid1, fofid2;
    
    // Compute the bins and fill the DD arrays with zeros
    for ( i = 0; i < NBIN; i++ ){
        DD_array[i] = 0, DD_1h[i] = 0, DD_2h[i] = 0;
        bins[i] = Dmin + Binsize * i;
        dist_array[i] = bins[i] + Binsize * .5;       
    }
    
    // DD computed
    for ( i = 0; i < ROWS; i++ ){
        
        r1[0] = X[i], r1[1] = Y[i], r1[2] = Z[i], fofid1 = FOFID[i];
        
        for ( j = i+1; j < ROWS; j++ ){
            
            r2[0] = X[j], r2[1] = Y[j], r2[2] = Z[j], fofid2 = FOFID[j];
            
            label = ToClass(logdist(r1, r2), Dmin, Dmax);
            if (label >= NBIN || label < 0){
                continue;
            }
            DD_array[label]++;
            if (fofid1 != fofid2){
                DD_2h[label]++;
                continue;
            }
            DD_1h[label]++;    
        }
    }
    
    double logxi[NBIN], logxi_1h[NBIN], logxi_2h[NBIN], logsuma[NBIN];
    double rr,xi, xi_1h, xi_2h;
    float r;
    
    // RR and logxi compute
    for ( i = 0, r = Dmin; i < NBIN; i++, r+=Binsize ){
        
        rr = RR(r, ROWS, 0.1);
        xi = 2*DD_array[i]/rr - 1;
        xi_1h = 2*DD_1h[i]/rr - 1;
        xi_2h = 2*DD_2h[i]/rr - 1;
        logxi[i] = log10(xi);
        logxi_1h[i] = log10(xi_1h);
        logxi_2h[i] = log10(xi_2h);
        logsuma[i] = log10(1 + xi_1h + xi_2h);
    }    
    

    // Write the data in a file
    FILE *data_file2;
    data_file2 = fopen("xi_data_10_2terms.txt", "w");
    
    for (i = 0; i < NBIN; i++){
        
        fprintf(data_file2, "%f %lf %lf %lf %lf", dist_array[i], logxi_1h[i], logxi_2h[i], logxi[i], logsuma[i]);
        fprintf(data_file2, "\n");
    }
      
    
    //for ( i = 0; i < NBIN; i++ ){
    //    printf("%f \n ", RR_array[i]);
    //}
    
    fclose(data_file);
    fclose(data_file2);
    
    // Code timing
    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("--- %.3f: loaded ---\n", time_spent);
    
    return 0;
}