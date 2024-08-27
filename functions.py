
import numpy as np
import mpmath as mp

def matrix_normalizer(t):
    return t / np.amax(t)

def interaction(t): # mpmath
    return mp.log(t[0, 0] / t[0, 1]) / 2

def transfer_matrix(interaction): # mpmath
    j = mp.mpf(interaction)
    t = mp.matrix([[mp.exp(j), mp.exp(-j)],
                   [mp.exp(-j), mp.exp(j)]])
    return matrix_normalizer(t)

def transfer_matrices(lattice_size, interaction, aferro_concentration):
    ferro  = transfer_matrix(interaction)
    aferro = transfer_matrix(-interaction)
    return [ferro for _ in range(int((1 - aferro_concentration) * lattice_size))] + [aferro for _ in range(int(aferro_concentration * lattice_size))]

def mp_multiply(t1, t2): # mpmath
    n = len(t1)
    t = mp.matrix(n)
    for i in range(n):
        for j in range(n):
            t[i, j] = t1[i, j] * t2[i, j]
    return matrix_normalizer(t)

def bond_moving(matrices): # mpmath
    t = matrices[0]
    n = len(matrices)
    for i in range(n - 1):
        t = mp_multiply(t, matrices[i + 1])
    return matrix_normalizer(t)

def decimation(matrices):
    t = matrices[0] * matrices[1] # mpmath (np.dot karşılığı)
    t = matrix_normalizer(t)    
    t = t * matrices[2]
    return matrix_normalizer(t)

def renormalize(matrices):
    np.random.seed(19)
    N = len(matrices)
    renormalized = []
    for k in range(N):
        
        random = []
        for _ in range(27):
            i = np.random.randint(1, N)
            random.append(matrices[i])

        # Renormalize
        bm1 = bond_moving(random[:9])
        bm2 = bond_moving(random[9:18])
        bm3 = bond_moving(random[18:])
        dm = decimation([bm1, bm2, bm3])
        renormalized.append(dm)
    
    return renormalized

def transfer_matrix_counter(matrices): # mpmath
    zero, one = 0.0001, 0.9999
    ferro, aferro, disorder, outofsink = 0, 0, 0, 0
    
    for i in range(len(matrices)):
        left = matrices[i][0, 0]
        right = matrices[i][0, 1]
        
        if left > one and right < zero:
            ferro += 1
        elif left < zero and right > one:
            aferro += 1
        elif left > one and right > one:
            disorder += 1
        else:
            outofsink += 1

    return ferro, aferro, disorder, outofsink

def phase_sink(interaction, aferro_concentration, lattice_size, search):

    limsup = 0.95 * lattice_size
    liminf = 0.05 * lattice_size

    matrices = transfer_matrices(lattice_size, interaction, aferro_concentration)
    phase = "undetermined"
    
    #k = 0
    #while phase == "undetermined":
    #    k += 1
    #    matrices = renormalize(matrices)
    #    ferro, aferro, disorder, outofsink = transfer_matrix_counter(matrices)
        
    #    if disorder > limsup and outofsink < liminf:
    #        phase = "disorder"
    #    if disorder < liminf and outofsink < liminf:
    #        phase = "spin-glass"

    for k in range(1, 30):
        matrices = renormalize(matrices)
        ferro, aferro, disorder, outofsink = transfer_matrix_counter(matrices)

        if disorder > limsup and outofsink < liminf:
            phase = 'disorder'
            if search == 'disorder':
                break
            else: pass
            
        elif ferro > limsup and outofsink < liminf:
            phase = 'ferro'
            if search == 'ferro':
                break
            else: pass

        else: pass

    return k, phase

def critical_point(temparature, aferro_concentration, lattice_size, tolerance, search_direction):
    """
    Calculates the critical point for a given system.
    Parameters:
    - temparature (float): The temperature of the system.
    - aferro_concentration (float): The concentration of aferro particles in the system.
    - lattice_size (int): The size of the lattice.
    - tolerance (float): The tolerance for convergence.
    - search_direction (str): The direction of the search ('vertical' or 'horizontal').
    Returns:
    - critical_point (float): The calculated critical point of the system.
    """
    T = temparature
    p = aferro_concentration
    N = lattice_size
    e = tolerance
    
    # Vertical search
    if search_direction == 'vertical':
        Ti = T
        phase = 'disorder'
    
        print("Initial search:")
        while phase == "disorder":
            Thigh = Ti
            print("Thigh =", Thigh)
            Ti -= 1
            _, phase = phase_sink(1/Ti, p, N, 'disorder')
            
        Tlow = Ti
        print("Tlow =", Tlow, end="\n\n")

        condition = True
        while condition:
            Tmid = (Tlow + Thigh) / 2
            print(f"Tlow = {round(Tlow, 4)}, Tmid = {round(Tmid, 4)}, Thigh = {round(Thigh, 4)} | e = {round(Thigh - Tlow, 4)}")

            _, phase = phase_sink(1/Tmid, p, N, 'disorder')
            if phase == 'disorder':
                Thigh = Tmid
            else:
                print("phase: not disorder!")
                Tlow = Tmid

            condition = abs(Thigh - Tlow) > e

        print(f"Tlow = {round(Tlow, 4)}, Tmid = {round(Tmid, 4)}, Thigh = {round(Thigh, 4)} | e = {round(Thigh - Tlow, 4)}")
        critical_point = (Thigh + Tlow) / 2
    
    # Horizontal search
    if search_direction == 'horizontal':
        pi = p
        phase = "ferro"
        
        print("Initial search:")
        while phase == "ferro":
            plow = pi
            print("plow =", plow)
            pi += 0.1
            _, phase = phase_sink(1/T, pi, N, 'ferro')
        phigh = pi
        print("phigh =", phigh, end="\n\n")

        condition = True
        while condition:
            pmid = (plow + phigh) / 2
            print(f"plow = {round(plow, 4)}, pmid = {round(pmid, 4)}, phigh = {round(phigh, 4)} | e = {round(phigh - plow, 4)}")

            _, phase = phase_sink(1/T, pmid, N, 'ferro')
            if phase == 'ferro':
                plow = pmid
            else:
                print("phase: not ferro!")
                phigh = pmid

            condition = abs(phigh - plow) > e
        
        print(f"plow = {round(plow, 4)}, pmid = {round(pmid, 4)}, phigh = {round(phigh, 4)} | e = {round(phigh - plow, 4)}")
        critical_point = (phigh + plow) / 2

    print()
    print(f"Critical point = {critical_point}")
    return critical_point



