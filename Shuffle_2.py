import numpy as np
import random

# np.save('../Data/G13_millennium/V.npy', V) It has the array of director vectors
# np.save('../Data/G13_millennium/r.npy', r_vector) # It contains the suffled catalogue

# Todos los datos estan previamente ordenados de acuerdo al fofid y al bin de masa
log_centralmvir = np.log10(np.load('../Data/G13_millennium/centralmvir.npy') + 1e-10) + 10
type_gal = np.load('../Data/G13_millennium/type.npy')
fofid = np.load('../Data/G13_millennium/fofid.npy')
x = np.load('../Data/G13_millennium/x.npy')
y = np.load('../Data/G13_millennium/y.npy')
z = np.load('../Data/G13_millennium/z.npy')
r_vector = np.array([x,y,z]).T
del x,y,z

logm_min = 10.
logm_max = 15.
NBIN = 30
mass_index = ((log_centralmvir - logm_min)/(logm_max - logm_min) * NBIN).astype(int)

# Elimino las filas correspondientes a masas que superen los 10^15 o sean menores a 10^10
mask = (mass_index >= 0) & (mass_index < NBIN)
mass_index = mass_index[mask]
type_gal = type_gal[mask]
fofid = fofid[mask]
r_vector = r_vector[mask]
del mask

N,j = 0,0
V = np.zeros((len(mass_index),3))  # Array con los vectores directores
len_data = len(mass_index)         # Numero de objetos en el catálogo

# ============================== VECTOR DIRECTOR =======================================
# El for calcula el vector director para cada una de las galaxias del catálogo SAM.
# Los vectores directores se calculan restando las posiciones entre las centrales conse-
# cutivas (y que corresponden a un mismo bin de masa). Este vector será el mismo para 
# todas las satelites que se encuentren en el halo analizado. 
# ======================================================================================
for i in range(NBIN): # Recorre todos los bines de masa
    
    
    # Me entrega las posiciones (en el arreglo) de las centrales y que corresponden al bin de masa i
    central_positions = np.where((type_gal == 0) & (mass_index == i))[0] 
    
    # Numero de centrales en el bin de masa i
    num_halos = len(central_positions)
      
    shuffled_pos = np.array(random.sample(central_positions, num_halos))

    # Recorre cada uno de los halos en el bin de masa i
    for j in range(num_halos):
        
        # Posicion (en el arreglo) de una central que se encuentra en un halo de bin de masa i
        pos_i = shuffled_pos[j] 
        
        # Si llegue al ultimo halo entonces usar el primer halo como pos_f
        if j == num_halos - 1: pos_f = shuffled_pos[0] 
         
        # Posicion (en el arreglo) de la siguiente central
        else: pos_f = shuffled_pos[j+1]  
        
        # Vector director entre estas dos centrales
        v_dir = r_vector[pos_f] - r_vector[pos_i]
        
        cond = True
        
        # Recorre todas las galaxias (incluyendo la misma central) del halo
        while cond:
            V[pos_i + N] = v_dir                          # Guarda el vector director para la galaxia (notar que empieza con la central)
            N += 1                                        # Numero de galaxias en el halo
            if pos_i + N == len_data: break               # No se ejecuta solo cuando llegue al final del catálogo
            cond = (fofid[pos_i+N] == fofid[pos_i])       # ¿Sigo en el mismo halo?
            k += 1                                        # Numero total de galaxias en el catálogo (al final es igual a len(mass_index))
            
        N = 0                                             
# ========================= SHUFFLING ======================================== # TODO MAKE IT MORE EFFICIENT!
# Suma los vectores directores a todas las posiciones de las galxias. Asi, cada
# galaxia es desplazada al halo que le sigue en el array. 
# ============================================================================ 

# Correccion para posiciones resultantes x,y o z mayores que 500 o menores que 0 
def isout (r,vr):
    if r + vr < 0:
        return 1
    elif r + vr > 500: 
        return -1
    else:
        return 0

Lbox = 500.   

# Desplazamiento de todas las galaxias
for i in range(len_data):
    
    x,y,z = r_vector[i]
    dx,dy,dz  = V[i]

    r_vector[i][0] = x + dx + isout(x,dx)*Lbox
    r_vector[i][1] = y + dy + isout(y,dy)*Lbox
    r_vector[i][2] = z + dz + isout(z,dz)*Lbox
    