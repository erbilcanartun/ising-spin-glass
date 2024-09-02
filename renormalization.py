
import numpy as np
import mpmath as mp

def matrix_normalizer(t):
    """
    Normalize a matrix by dividing each element by the maximum value in the matrix.
    Parameters:
    t (numpy.ndarray): The input matrix.
    Returns:
    numpy.ndarray: The normalized matrix.
    """

    return t / np.amax(t)

def interaction(t):
    """
        Calculate the interaction between two elements of a matrix.
        Parameters:
        t (mpmath.matrix): A 2x2 matrix representing the transfer matrix.
        Returns:
        mpmath.mpf: The calculated interaction value or coupling constant.
        """
    
    return mp.log(t[0, 0] / t[0, 1]) / 2

def transfer_matrix(interaction):
    """
    Calculates the transfer matrix for the Ising spin glass model.
    Parameters:
    interaction (float): The interaction strength between spins.
    Returns:
    numpy.ndarray: The normalized transfer matrix.
    """

    j = mp.mpf(interaction)
    t = mp.matrix([[mp.exp(j), mp.exp(-j)],
                   [mp.exp(-j), mp.exp(j)]])
    return matrix_normalizer(t)

def transfer_matrices(lattice_size, interaction, aferro_concentration):
    """
    Generates a list of transfer matrices based on the lattice size, interaction strength, and antiferromagnetic bond concentration.
    Parameters:
    lattice_size (int): The size of the lattice.
    interaction (float): The strength of the interaction.
    aferro_concentration (float): The concentration of antiferromagnetic interactions.
    Returns:
    list: A list of transfer matrices.
    """

    ferro  = transfer_matrix(interaction)
    aferro = transfer_matrix(-interaction)
    return [ferro for _ in range(int((1 - aferro_concentration) * lattice_size))] + [aferro for _ in range(int(aferro_concentration * lattice_size))]

def mp_multiply(t1, t2):
    """
    Multiplies two matrices element-wise and returns the result.
    Parameters:
    t1 (mp.matrix): The first matrix.
    t2 (mp.matrix): The second matrix.
    Returns:
    mp.matrix: The element-wise product of t1 and t2.
    """

    n = len(t1)
    t = mp.matrix(n)
    for i in range(n):
        for j in range(n):
            t[i, j] = t1[i, j] * t2[i, j]
    return matrix_normalizer(t)

def bond_moving(matrices):
    """
    Perform bond moving operation on a list of matrices.
    Parameters:
    matrices (list): A list of matrices.
    Returns:
    matrix: The result of the bond moving operation.
    """

    t = matrices[0]
    n = len(matrices)
    for i in range(n - 1):
        t = mp_multiply(t, matrices[i + 1])
    return matrix_normalizer(t)

def decimation(matrices):
    """
    Applies decimation to a list of matrices.
    Parameters:
    matrices (list): A list of matrices.
    Returns:
    matrix: The result of applying decimation to the matrices.
    """

    t = matrices[0] * matrices[1]
    t = matrices[0] * matrices[1]
    t = matrix_normalizer(t)    
    t = t * matrices[2]
    return matrix_normalizer(t)

def renormalize(matrices, method='bd'):
    """
    Renormalizes a list of transfer matrices. This renormalization
    method uses bond-moving operation first, then applies decimation.
    Parameters:
    matrices (list): A list of matrices.
    Returns:
    list: A list of renormalized matrices.
    """
    
    # Note that the seed is fixed for reproducibility
    np.random.seed(19)

    N = len(matrices)
    renormalized = []

    # Bond-moving and decimation operations, respectively
    if method == 'bd':
        # Iterate through the renormalization process for whole transfer matrix set
        for k in range(N):
            
            # Randomly select 27 matrices:
            # b^(d-1) matrices for b=3 and d=3
            random = []
            for _ in range(27):
                i = np.random.randint(1, N)
                random.append(matrices[i])

            # Bond moving operation
            bm1 = bond_moving(random[:9])
            bm2 = bond_moving(random[9:18])
            bm3 = bond_moving(random[18:])
            
            # Decimation operation
            dm = decimation([bm1, bm2, bm3])

            # Append the renormalized matrix to the list
            renormalized.append(dm)

    # Decimation and bond-moving operations, respectively
    if method == 'db':
        # Iterate through the renormalization process for whole transfer matrix set
        for k in range(N):
            
            # Randomly select 27 matrices:
            # b^(d-1) matrices for b=3 and d=3
            random = []
            for _ in range(27):
                i = np.random.randint(1, N)
                random.append(matrices[i])

            # Bond moving operation
            bm1 = bond_moving(random[:9])
            bm2 = bond_moving(random[9:18])
            bm3 = bond_moving(random[18:])
            # Decimation operation
            dm = decimation([bm1, bm2, bm3])
            # Append the renormalized matrix to the list
            renormalized.append(dm)
    
    return renormalized

def transfer_matrix_counter(matrices):
    """
    Counts the number of occurrences of different types of transfer matrices.
    Parameters:
    matrices (list): A list of matrices.
    Returns:
    tuple: A tuple containing the count of ferro, aferro, disorder, and outofsink transfer matrices.
    """

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

def phase_sink(interaction, aferro_concentration, lattice_size, search, method):
    """
    Determines the phase of a system based on the given parameters.
    Parameters:
    interaction (float): The interaction strength between spins.
    aferro_concentration (float): The concentration of antiferromagnetic spin interactions.
    lattice_size (int): The size of the lattice.
    search (str): The phase to search for ('disorder' or 'ferro').
    Returns:
    tuple: A tuple containing the number of iterations performed (k) and the determined phase of the system.
    """

    limsup = 0.95 * lattice_size
    liminf = 0.05 * lattice_size

    matrices = transfer_matrices(lattice_size, interaction, aferro_concentration)
    phase = "undetermined"
    
    # Iterate through the renormalization process
    for k in range(1, 30):
        matrices = renormalize(matrices, method)
        ferro, _, disorder, outofsink = transfer_matrix_counter(matrices)

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

def critical_point(temparature, aferro_concentration, lattice_size, tolerance, search_direction, method):
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
    
        # Initial search to determine the range in which the critical point lies
        print("Initial search:")
        while phase == "disorder":
            Thigh = Ti
            print("Thigh =", Thigh)
            Ti -= 1
            _, phase = phase_sink(1/Ti, p, N, 'disorder', method)
            
        Tlow = Ti
        print("Tlow =", Tlow, end="\n\n")

        # Binary search to find the critical point
        condition = True
        while condition:
            Tmid = (Tlow + Thigh) / 2
            print(f"Tlow = {round(Tlow, 4)}, Tmid = {round(Tmid, 4)}, Thigh = {round(Thigh, 4)} | e = {round(Thigh - Tlow, 4)}")

            _, phase = phase_sink(1/Tmid, p, N, 'disorder', method)
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
        
        # Initial search to determine the range in which the critical point lies
        print("Initial search:")
        while phase == "ferro":
            plow = pi
            print("plow =", plow)
            pi += 0.1
            _, phase = phase_sink(1/T, pi, N, 'ferro', method)
        phigh = pi
        print("phigh =", phigh, end="\n\n")

        # Binary search to find the critical point
        condition = True
        while condition:
            pmid = (plow + phigh) / 2
            print(f"plow = {round(plow, 4)}, pmid = {round(pmid, 4)}, phigh = {round(phigh, 4)} | e = {round(phigh - plow, 4)}")

            _, phase = phase_sink(1/T, pmid, N, 'ferro', method)
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



