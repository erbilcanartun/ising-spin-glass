# Spin-1/2 Ising Spin Glass (3D)

This repository presents a demonstration of global phase diagram calculations using renormalization-group theory for spin$-1/2$ Ising spin glasses in three spatial dimensions. The calculations are performed exactly on a three-dimensional hierarchical model and, in parallel, using the Migdal-Kadanoff approximation on the cubic lattice.

The renormalization process can be implemented in two ways: either applying a bond-moving operation followed by decimation (bd), or performing decimation first, then bond-moving (db). Both methods are included in this project. Therefore, the only distinction between files labeled 'bd' and 'db' lies in the renormalization technique employed. The sole observable consequence of this difference is the scaling of the critical temperatures.

## Contents

- [Theoretical Background](#theoretical-background)
    - [Ising Model](#ising-model)
    - [Spin Glasses](#spin-glasses)
    - [Renormalization Group](#renormalization-group)
- [Methodology](#methodology)

## Theoretical Background

### Ising Model

The Ising model is a simplified mathematical model used to study ferromagnetism. It consists of a lattice of spins, where each spin can take one of two possible values, usually represented as $+1$ (up) or $-1$ (down). The spins interact with their nearest neighbors via an interaction energy that favors either alignment (ferromagnetic interaction) or anti-alignment (antiferromagnetic interaction). 

### Spin Glasses

Spin glasses are magnetic systems characterized by a combination of disorder and frustration. Disorder refers to the random distribution of ferromagnetic and antiferromagnetic interactions between spins. Frustration arises because it is impossible to satisfy all the interactions simultaneously. 

### Renormalization Group

The renormalization group (RG) is a powerful theoretical framework used to study the behavior of systems at different length scales. It involves systematically integrating out degrees of freedom at short length scales to obtain an effective description at longer length scales. 

## Methodology

The code in this repository calculates the phase diagram of the Ising spin glass model in three dimensions ($d=3$) using a renormalization-group approach.

The spin$-1/2$ Ising model is defined by the Hamiltonian:  
    -\beta \mathcal{H} = \sum_{\langle ij\rangle} {E(s_i,s_j)} = \sum_{\langle ij\rangle} {J_{ij} s_i s_j}

and the interactions are represented using a transfer matrix:  
    $T(s_i,s_j) = e^{E(s_i,s_j)}$

As the local renormalization-group transformation, the Migdal-Kadanoff approximate transformation and, equivalently, the exact transformation for the $d=3$ hierarchical lattice is used (Figure). The length rescaling factor of $b=3$ is used to preserve under renormalization group the ferromagnetic-antiferromagnetic symmetry of the system. This local transformation consists in bond moving followed by decimation.

![renormalization](images/renormalization.png "Figure: (a) RG transformation for the
length-rescaling factor of b=3. In this intuitive approximation, bond moving is followed by decimation. (b) Exact RG transformation of the d = 3, b = 3 hierarchical lattice for which the Migdal-Kadanoff RG recursion relations are exact. The construction of a hierarchical lattice proceeds in the opposite direction of its RG solution.")

The quenched randomness is included by keeping, as a distribution, `lattice_size` sets of the nearest-neighbor interaction energies $E(s_i,s_j)$. At the beginning of each renormalization-group trajectory, this distribution is formed from the double-delta distribution characterized by interactions $\pm J$ with probabilities $p$, $(1−p)$. 

- Bond-moving:

    Multiply elements at the same position of $b^{(d-1)}$ randomly chosen transfer matrices. In this case, $b$ (length rescaling factor) is $3$ and $d$ is $3$, so $9$ transfer matrices are multiplied.

- Decimation:

    Perform matrix multiplication of three randomly chosen bond-moved transfer matrices.


After all operations, a distribution of `lattice_size` renormalized transfer matrices is generated. Phases are determined by following trajectories to their asymptotic limit: The asymptotic limit transfer matrices of trajectories starting in the ferromagnetic phase all have $1$ in the corner diagonals and $0$ at all other positions. The asymptotic limit transfer matrices of trajectories starting in the antiferromagnetic phase all have $1$ in the corner antidiagonals and $0$ at all other positions. The asymptotic limit transfer matrices of trajectories starting in the spin-glass phase all have $1$ in the corner diagonals $(s_i=s, s_j=s)$ and $(s_i=-s, s_j=-s)$ or in the corner antidiagonals $(s_i=s, s_j=-s)$ and $(s_i=-s, s_j=-s)$, and $0$ at all other positions. The asymptotic limit transfer matrices of trajectories starting in the disordered phase all have $1$ at all positions.

Phase diagrams are obtained by numerically determining the boundaries, in the unrenormalized system, of these asymptotic flows.

