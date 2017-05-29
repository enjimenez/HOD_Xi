#include <math.h> 
# define Lbox 500 // Size of millennium simulation

double Test_C(double * A, double * B,int Len_A,int Len_B){
	int i,j;
	double SUMA = 0;
	for(i = 0; i < Len_A; i++)	for(j = 0; j < Len_B; j++)	SUMA = SUMA+A[i]+B[j];
	return SUMA;
}

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

int ToClass (float d,  float Xmin, float Xmax, int NBIN){
    
    int label;
    label = ((d-Xmin)/(Xmax - Xmin)) * NBIN;
    return label;
}

void corr_function_C(int * DD_array, double * X, double * Y, double * Z, int len_data, float Dmin, float Dmax, int NBIN){
    
    int i,j, label;
    float r1[3], r2[3];
    for ( i = 0; i < len_data; i++ ){
        r1[0] = X[i], r1[1] = Y[i], r1[2] = Z[i];
         for ( j = i+1; j < len_data; j++ ){
             r2[0] = X[j], r2[1] = Y[j], r2[2] = Z[j];
             label = ToClass(logdist(r1, r2), Dmin, Dmax, NBIN);
             if (label >= NBIN || label < 0){
                continue;
             }
             DD_array[label]++;
        }
    }
}

void SMF_C(double * left, double * stellarmass, double * SMF, int NBIN, int len_data){
    
    float mass_cut;
    int i,j, cont = 0;
    for (i = 0; i < NBIN; i++){

        mass_cut = left[i];
       
        for (j = 0; j < len_data; j++){
            if (stellarmass[j] > mass_cut){
                cont++;
            }
            else{
                continue;
            }
        double temp = cont;
        SMF[i] = temp / len_data;
        cont = 0;
        }
    }
}

int ilimit (int x, int iNBin){
    if (x<0) return x + iNBin;
    if (x>=iNBin) return x - iNBin;
    return x;
}

void corr_function_C2(int * DD_array, double * X, double * Y, double * Z, int len_data, float Dmin, float Dmax, int NBIN, int iNBin, int ilim){
    
    printf("partio!\n");
    int ilim2, i,j, k,label, tx, ty, tz, x, y, z, first = -99;
    float r1[3], r2[3];
    
    int * ll; 
    int * ix;
    int * iy;
    int * iz;
    printf("partio 2!\n");
    ll = (int*) calloc (len_data, sizeof(int));
    ix = (int*) calloc (len_data, sizeof(int));
    iy = (int*) calloc (len_data, sizeof(int));
    iz = (int*) calloc (len_data, sizeof(int));
    int lfirst [iNBin][iNBin][iNBin];
     printf("partio 3!\n");
    
    for (i=0; i< iNBin; i++) for(j=0; j< iNBin; j++) for(k=0; k< iNBin; k++) lfirst[i][j][k] = -99;
    for ( i= 0; i< len_data; i++){
        
        ix[i] = (int) (X[i]/Lbox*iNBin);
        iy[i] = (int) (Y[i]/Lbox*iNBin);
        iz[i] = (int) (Z[i]/Lbox*iNBin);
        lfirst[ix[i]][iy[i]][iz[i]] = i;
        
        
    }
    for ( i= 0; i< len_data; i++){
        
        ll[i] = lfirst[ix[i]][iy[i]][iz[i]]; 
        lfirst[ix[i]][iy[i]][iz[i]] = i;
    }
     
    printf("llegue al triple for!\n");
    for ( i = 0; i < len_data; i++ ){
        r1[0] = X[i], r1[1] = Y[i], r1[2] = Z[i];
        
         for ( tx = ix[i] - ilim; tx <= ix[i] + ilim; tx++ ){
             
             ilim2 = (int) sqrt((float) ilim*ilim - (ix[i]-tx)*(ix[i]-tx));
             
             for ( ty = iy[i] - ilim2; ty <= iy[i] + ilim2; ty++ ){
                 
                 for ( tz = iz[i] - ilim2; tz <= iz[i] + ilim2; tz++ ){
                
                     x = ilimit(tx, iNBin);
                     y = ilimit(ty, iNBin);
                     z = ilimit(tz, iNBin);
                     
                     j = lfirst[x][y][z];
                     
                     while (j != -99 && first != j){
                         r2[0] = X[j], r2[1] = Y[j], r2[2] = Z[j];
                         label = ToClass(logdist(r1, r2), Dmin, Dmax, NBIN);
                         if (label < NBIN && label >= 0){
                            DD_array[label]++;
                        }
                         j = ll[j];
                         first = lfirst[x][y][z]; 
                    }
           
                 }
             }
        }
    }
}

void corr_function_C3(int * DD_array, double * X, double * Y, double * Z, int len_data, float Dmin, float Dmax, int NBIN, int iNBin, int ilim, int CPU, int NCPU){
    
    int imin = (int) (len_data*(CPU - 1.)/(float) NCPU);
    int imax = (int) (len_data*(CPU)/(float) NCPU);
    
    printf("partio!\n");
    int ilim2, i,j, k,label, tx, ty, tz, x, y, z, first = -99;
    float r1[3], r2[3];
    
    int * ll;
    int * ix;
    int * iy;
    int * iz;
    printf("partio 2!\n");
    ll = (int*) calloc (len_data, sizeof(int));
    ix = (int*) calloc (len_data, sizeof(int));
    iy = (int*) calloc (len_data, sizeof(int));
    iz = (int*) calloc (len_data, sizeof(int));
    int lfirst [iNBin][iNBin][iNBin];
     printf("partio 3!\n");
    
    for (i=0; i< iNBin; i++) for(j=0; j< iNBin; j++) for(k=0; k< iNBin; k++) lfirst[i][j][k] = -99;
    for ( i= 0; i< len_data; i++){
        
        ix[i] = (int) (X[i]/Lbox*iNBin);
        iy[i] = (int) (Y[i]/Lbox*iNBin);
        iz[i] = (int) (Z[i]/Lbox*iNBin);
        lfirst[ix[i]][iy[i]][iz[i]] = i;
        
        
    }
    for ( i= 0; i< len_data; i++){
        
        ll[i] = lfirst[ix[i]][iy[i]][iz[i]]; 
        lfirst[ix[i]][iy[i]][iz[i]] = i;
    }
     
    printf("llegue al triple for!\n");
    for ( i = imin; i < imax; i++ ){
        r1[0] = X[i], r1[1] = Y[i], r1[2] = Z[i];
        
         for ( tx = ix[i] - ilim; tx <= ix[i] + ilim; tx++ ){
             
             ilim2 = (int) sqrt((float) ilim*ilim - (ix[i]-tx)*(ix[i]-tx));
             
             for ( ty = iy[i] - ilim2; ty <= iy[i] + ilim2; ty++ ){
                 
                 for ( tz = iz[i] - ilim2; tz <= iz[i] + ilim2; tz++ ){
                
                     x = ilimit(tx, iNBin);
                     y = ilimit(ty, iNBin);
                     z = ilimit(tz, iNBin);
                     
                     j = lfirst[x][y][z];
                     
                     while (j != -99 && first != j){
                         r2[0] = X[j], r2[1] = Y[j], r2[2] = Z[j];
                         label = ToClass(logdist(r1, r2), Dmin, Dmax, NBIN);
                         if (label < NBIN && label >= 0){
                            DD_array[label]++;
                        }
                         j = ll[j];
                         first = lfirst[x][y][z]; 
                    }
           
                 }
             }
        }
    }
}

void Ngal_C(double * ngal, double * halomass, int * mass_index, float * M, float * N, int len_data){
   
    int i,index;
    float m,n;
    double power;
    for (i = 0; i < len_data; i++ ){
        index = mass_index[i];
        m = M[index];
        n = N[index];
        power = (m*halomass[i]) + n;
        ngal[i] = pow(10.0, power);
    }
}

void Ngal_C2(double * ngal_all, double * ngal_cen, double * ngal_sat, double * halomass, int * mass_index, float * M_all, float * N_all, float * M_cen, float * N_cen, float * M_sat, float * N_sat, int len_data){
   
    int i,index;
    float m_all, m_cen, m_sat, n_all, n_cen, n_sat, hmass;
    double power_all, power_cen, power_sat;
    for (i = 0; i < len_data; i++ ){
        index = mass_index[i];
        hmass = halomass[i];
        m_all = M_all[index], n_all = N_all[index];
        m_cen = M_cen[index], n_cen = N_cen[index];
        m_sat = M_sat[index], n_sat = N_sat[index];
        power_all = (m_all * hmass) + n_all;
        power_cen = (m_cen * hmass) + n_cen;
        power_sat = (m_sat * hmass) + n_sat;
        ngal_all[i] = pow(10, power_all);
        ngal_cen[i] = pow(10, power_cen);
        ngal_sat[i] = pow(10, power_sat);
    }
}
