import numpy as np
import scipy
import scipy.special
import scipy.sparse.linalg as LA
import scipy.sparse as sparse
import copy
import warnings


def get_1d_coords(p, i, j):
    """
    Finds index of site in 1d chain from 2d lattics based on snake
    decomposition, ie 1d coords on 2d lattice look like:
        0, 1, 2, 3,
        7, 6, 5, 4,
        8, 8, 10, 11
    Args:
        p - dictionary that contains the relevant system parameters
        i - row index on 2d lattice
        j - column index on 2d lattice
    Returns:
        reshaped_i - index of site on 1d chain
    """
    if i % 2 == 0:
        reshaped_i = i * p['W'] + j
    else:
        reshaped_i = i * p['W'] + (p['W'] - j - 1)
    return reshaped_i


def get_2d_coords(p, i):
    """
    Finds coordinates of site in 2d lattice from 1d chain based on snake
    decomposition, ie 1d coords on 2d lattice look like:
        0, 1, 2, 3,
        7, 6, 5, 4,
        8, 8, 10, 11
    Args:
        p - dictionary that contains the relevant system parameters
        i - index of site on 1d chain
    Returns:
        reshaped_i - row index on 2d lattice
        reshaped_j - column index on 2d lattice
    """
    reshaped_i = i // p['W']
    if reshaped_i % 2 == 0:
        reshaped_j = i % p['W']
    else:
        reshaped_j = p['W'] - i % p['W'] - 1
    return reshaped_i, reshaped_j


def reshape_projected_states(p, projected_states):
    """
    Reshapes (2xN) projected_states into (2xLxW) projected_states_lattice
    Args:
        p - dictionary that contains the relevant system parameters    Returns:
        projected_states - (2xN) list with the fermion occupation for
            spcies 1,2
    Returns:
        projected_states_lattice - (2xLxW) list with the fermion occupation for
            spcies 1,2
    """
    projected_states_lattice = np.empty((2, p['L'], p['W']))
    for spin_ind in [0, 1]:
        for i, s in enumerate(projected_states[spin_ind]):
            ind = tuple([spin_ind] + list(get_2d_coords(p, i)))
            projected_states_lattice[ind] = s
    return projected_states_lattice


def reshape_state(p, statelist):
    """
    Reshapes 1d statelist into 2d statelattice
    Args:
        p - dictionary that contains the relevant system parameters    Returns:
        statelist - 1d fermion configuration
    Returns:
        statelattice - 2d fermion configuration
    """
    statelattice = np.empty((p['L'], p['W']))
    for i, s in enumerate(statelist):
        statelattice[get_2d_coords(p, i)] = s
    return statelattice


def flatten_state(p, statelattice):
    """
    Reshapes 2d statelattice into 1d statelist
    Args:
        p - dictionary that contains the relevant system parameters    Returns:
        statelattice - 2d fermion configuration
    Returns:
        statelist - 1d fermion configuration
    """
    statelist = np.empty(p['N'])
    for i in range(p['N']):
        s = statelattice[get_2d_coords(p, i)]
        statelist[i] = s
    return statelist


def convert_to_base(num, base):
    """
    Converts integer num into base b format
    ie number to configuration
    Args:
        num - intger to be converted
        base - base to be converted into
    Returns:
        base b representation of num
    """
    convStr = "0123"
    if num < base:
        return str(num)
    else:
        return convert_to_base(num // base, base) + convStr[num % base]


def state_to_int(p, statelist):
    """
    Converts array of fermion-configuration into integer
    Args:
        p - dictionary that contains the relevant system parameters
        statelist - fermion configuration
    Returns:
        out - integer corresponding to state
    """
    # construct unique integer for the fermion configuration defined
    # in statelist
    out = 0
    for ind, val in enumerate(statelist):
        out += val * 4**(p['N'] - ind - 1)
    out = int(out)
    return out


def int_to_state(p, state_number):
    """
    Converts integer to array of fermion-configuration
    Args:
        p - dictionary that contains the relevant system parameters
        state_number - integer corresponding to state
    Returns:
        statelist - fermion configuration
    """

    # convert integer to spin configuration of length L
    statelist = [0] * int(p['N'])
    # get base4 representation of the state number
    state = convert_to_base(state_number, 4)
    index = 1
    for i in range(len(state) - 1, -1, -1):
        statelist[int(p['N'] - index)] = int(state[i])
        index += 1
    return statelist


def count_particles(p, state_number):
    """
    Counts number of particles N1, N2, N12
    Args:
        p - dictionary that contains the relevant system parameters
        state_number - integer corresponding to state
    Returns:
        N1 - number of fermions in 1
        N2 - number of fermions in 2
        N12 - number of doublons (12)
    """
    # convert integer to spin configuration of length L
    state = int_to_state(p, state_number)
    # count 'loose' fermions of type 1
    N1 = state.count(1)
    # count 'loose' fermions of type 2
    N2 = state.count(2)
    # count doublons of type '12'
    N12 = state.count(3)

    # to get total number N1 (N2), add N1 to N12 (N2 to N12)
    N1 += N12
    N2 += N12

    return N1, N2, N12


def state_in_subspace(p, state_number):
    """
    checks if state corresponding to state_number is in symmetry-allowed
    Hilbertspace
    Args:
        p - dictionary that contains the relevant system parameters
        state_number - integer corresponding to state
    Returns:
        accept - boolean if state is (dis)allowed
    """
    # check if state is in relevant subspace, depending on which Rabi-Terms
    # are turned on
    N1, N2, N12 = count_particles(p, state_number)

    # if individually N1 and N2 are not as set out in parameters p, discard
    # state (unless mu1/mu2 non-zero, then we're mixing differnet
    # Number-sectors)
    if (N1 != p['N1']) or (N2 != p['N2']):
        return False
    else:
        return True


def generate_state_table(p):
    """
    generates table of state-integers that are allowed by the symmetries
    of the model
    Args:
        p - dictionary that contains the relevant system parameters
    Returns:
        state_table - list of all state_numbers that belong to the relevant
            Hilbertspace
    """
    # generate list of state_numbers which are allowed by the symmetries
    state_table = []
    for i in range(int(4**p['N'])):
        if p['mu1'] != 0.0 and p['mu2'] != 0.0:
            # if we have non-zero chemical potential, we are considering all
            # (N1,N2) sectors, so need to include all states in state-table
            state_table.append(i)
        elif state_in_subspace(p, i):
            state_table.append(i)
    return state_table


def project_state_into_spinsectors(statelist):
    """
    projects the fermion configuration defined by state into the 2 fermion
    species' subsectors
    Note: 0 = Empty, 1 = Spin-Up, 2 = Spin-Down, 3 = Double (Spin-Up+Spin-Down)
    Args:
        statelist - fermion configuration
    Returns:
        projected_states - (2xN) list with the fermion occupation for
            spcies 1,2
    """
    # rewrite spin configuration into configurations for each species
    # seperately
    L = len(statelist)
    state1 = [0] * L
    state2 = [0] * L

    for i, si in zip(range(L), statelist):
        state1[i] = int((si == 1) or (si == 3))
        state2[i] = int((si == 2) or (si == 3))

    projected_states = np.array([state1, state2])
    return projected_states


def join_spinsectors_into_state(projected_states):
    """
    takes the projection of a state on the 2 spin sectors and joins them into a
    single 2-fermion state configuration
    Args:
        projected_states - (2xN) list with the femrion occupation for
            species 1,2
    Returns:
        statelist - fermion configuration
    """
    # join the binary spin representation for each spin into a single base8
    # representation
    statelist = np.array(projected_states[0]) + 2 * np.array(projected_states[1])
    return statelist


def fermisign(projected_states, i, j, sigma, tau=None):
    """
    Gounts how many fermions are sitting between sites i and j with spin sigma
    in the state defined by stateprojection
    Note: implicit spin-ordering for Fermi sign: all Up before all Down
    Args:
        projected_states - (2xN) list with the fermion occupation for
            species 1-2
        i - inital site
        j - final site
        sigma - initial spin
        tau (optional) - final spin, default same as initial spin
    Returns:
        fermisign - +/-1 for even/odd # of fermions
    """
    ni = 0
    if tau is None or tau == sigma:
        if i < j:
            ni += np.sum(projected_states[sigma - 1, i + 1:j])
        elif i > j:
            ni += 1
            ni += np.sum(projected_states[sigma - 1, j + 1:i])
    else:
        if sigma < tau:
            ni += np.sum(projected_states[sigma - 1, i + 1:])
            ni += np.sum(projected_states[tau - 1, :j])
        else:
            ni += 1
            ni += np.sum(projected_states[tau - 1, j + 1:])
            ni += np.sum(projected_states[sigma - 1, :i])
    return (-1)**ni


def hopping_matrix(p, state_table, bond_type=None):
    """
    Generates the full hopping matrix defined on the allowed
    Hilbertspace subspace
    Args:
        p - dictionary that contains the relevant system parameters
        state_table - list of all state_numbers that belong to the
            relevant Hilbertspace
    Returns:
        hopping - hopping matrix on the relevant Hilbertspace
    """
    dim = len(state_table)
    row = []
    col = []
    data = []
    t = [p['t1'], p['t2']]

    for In in range(dim):
        state = int_to_state(p, state_table[In])
        stateproj = project_state_into_spinsectors(state)

        for sigma in [1, 2]:
            sp_ind = sigma - 1
            for i in range(p['L']):
                for j in range(p['W']):
                    si = get_1d_coords(p, i, j)
                    horizontal = bond_type is None
                    vertical = bond_type is None
                    if bond_type == 1 and j % 2 == 0:
                        horizontal = True
                    elif bond_type == 3 and j % 2 == 1:
                        horizontal = True
                    if bond_type == 2 and i % 2 == 0:
                        vertical = True
                    elif bond_type == 4 and i % 2 == 1:
                        vertical = True

                    # horizontal hopping
                    if horizontal and (j != p['W'] - 1):
                        si_h = get_1d_coords(p, i, j + 1)
                        if stateproj[sp_ind][si] and not stateproj[sp_ind][si_h]:
                            out_state = copy.deepcopy(stateproj)
                            out_state[sp_ind][si] = 0
                            out_state[sp_ind][si_h] = 1
                            new_state = join_spinsectors_into_state(out_state)
                            Out = state_table.index(state_to_int(p, new_state))
                            sgn = fermisign(stateproj, si, si_h, sigma)
                            matrixelement = -1.0 * t[sp_ind] * sgn
                            # store matrix element
                            row.append(Out)
                            col.append(In)
                            data.append(matrixelement)
                            # add the h.c., hopping in the opposite direction
                            row.append(In)
                            col.append(Out)
                            data.append(np.conj(matrixelement))

                    # vertical hopping
                    if vertical and (i != p['L'] - 1):
                        si_v = get_1d_coords(p, i + 1, j)
                        if stateproj[sp_ind][si] and not stateproj[sp_ind][si_v]:
                            out_state = copy.deepcopy(stateproj)
                            out_state[sp_ind][si] = 0
                            out_state[sp_ind][si_v] = 1
                            new_state = join_spinsectors_into_state(out_state)
                            Out = state_table.index(state_to_int(p, new_state))
                            sgn = fermisign(stateproj, si, si_v, sigma)
                            matrixelement = -1.0 * t[sp_ind] * sgn
                            # store matrix element
                            row.append(Out)
                            col.append(In)
                            data.append(matrixelement)
                            # add the h.c., hopping in the opposite direction
                            row.append(In)
                            col.append(Out)
                            data.append(np.conj(matrixelement))

    hopping = sparse.csr_matrix((data, (row, col)),
                                shape=(dim, dim), dtype=complex)
    return hopping


def onsite_interaction_matrix(p, state_table):
    """
    generates the full on-site interaction matrix defined on the
    allowed Hilbertspace subspace
    Args:
        p - dictionary that contains the relevant system parameters
        state_table - list of all state_numbers that belong to the
            relevant Hilbertspace
    Returns:
        interaction - interaction matrix on the relevant Hilbertspace
    """
    dim = len(state_table)
    row = []
    col = []
    data = []

    for In in range(dim):
        state = int_to_state(p, state_table[In])
        stateprojection = project_state_into_spinsectors(state)

        # use bitwise-and (&) to find double occupancies
        state12 = np.array(stateprojection[0]) & np.array(stateprojection[1])

        # construct & store matrix element
        matrixelement = p['U12'] * np.sum(state12)
        if matrixelement != 0.0:
            row.append(In)
            col.append(In)
            data.append(matrixelement)

        del matrixelement

    interaction = sparse.csr_matrix((data, (row, col)),
                                    shape=(dim, dim), dtype=complex)
    return interaction


def chemical_potential_matrix(p, state_table):
    """
    generates the full chemical potential matrix defined on the
    allowed Hilbertspace subspace
    Note: chemical potential term in Hamiltonian is negative,
    i.e. H ~ - mu1*N1 - mu2*N2
    Args:
        p - dictionary that contains the relevant system parameters
        state_table - list of all state_numbers that belong to the relevant
            Hilbertspace
    Returns:
        mu - chemical potential matrix on the relevant Hilbertspace
    """
    dim = len(state_table)
    row = []
    col = []
    data = []

    for In in range(dim):
        state = int_to_state(p, state_table[In])
        stateprojection = project_state_into_spinsectors(state)

        matrixelement = 0.0
        for i in range(2):
            matrixelement += -1.0 * \
                p['mu' + str(i + 1)] * np.sum(stateprojection[i])

        # store matrix element
        if matrixelement != 0.0:
            row.append(In)
            col.append(In)
            data.append(matrixelement)

    mu = sparse.csr_matrix((data, (row, col)),
                           shape=(dim, dim), dtype=complex)
    return mu


def ni_matrix(p, site, species, state_table):
    """
    generates the matrix corresponding to the density operator on site i for
    fermions of type 'species'
    Args:
        p - dictionary that contains the relevant system parameters
        site - site on which density is to be evaluated
        species - spin species in question (1 or 2)
        state_table - list of all state_numbers that belong to the relevant
            Hilbertspace
    Returns:
        ni - n_{i,species} matrix on the relevant Hilbertspace
    """
    dim = len(state_table)
    row = []
    col = []
    data = []

    for In in range(dim):
        state = int_to_state(p, state_table[In])
        stateprojection = project_state_into_spinsectors(state)

        matrixelement = stateprojection[species - 1][site]
        # store matrix element
        if matrixelement != 0.0:
            row.append(In)
            col.append(In)
            data.append(matrixelement)

        del matrixelement

    ni = sparse.csr_matrix((data, (row, col)), shape=(dim, dim), dtype=complex)
    return ni


def CdagC_matrix(p, i, sigma, j, tau, state_table):
    """
    generates the matrix corresponding to the single particles correlations
    Cdag_{i,sigma}C_{j,tau}
    Args:
        p - dictionary that contains the relevant system parameters
        i - inital site
        sigma - intial spin (1 or 2)
        j - final site
        tau - final spin (1 or 2)
        state_table - list of all state_numbers that belong to the relevant
            Hilbertspace
    Returns:
        CdagC - C^{dagger}_{i,sigma}C_{i,tau} matrix on the relevant
            Hilbertspace
    """
    if i == j and sigma == tau:
        return ni_matrix(p, i, sigma, state_table)

    dim = len(state_table)
    row = []
    col = []
    data = []

    for In in range(dim):
        state = int_to_state(p, state_table[In])
        stateprojection = project_state_into_spinsectors(state)

        # only nonzero if (i,sigma) empty and (j,tau) occupied
        if stateprojection[tau - 1][j] and not stateprojection[sigma - 1][i]:
            out_state = copy.deepcopy(stateprojection)

            out_state[sigma - 1][i] = 1
            out_state[tau - 1][j] = 0

            new_state = join_spinsectors_into_state(out_state)
            Out = state_table.index(state_to_int(p, new_state))

            matrixelement = fermisign(stateprojection, i, j, sigma, tau)
            # store matrix element
            row.append(Out)
            col.append(In)
            data.append(matrixelement)

            del out_state, matrixelement

    CdagC = sparse.csr_matrix(
        (data, (row, col)), shape=(dim, dim), dtype=complex)
    return CdagC


def nk_matrix(p, m, species, state_table):
    """
    generates the matrix corresponding to the density operator in momentum
    space for fermions of type 'species'
    Args:
        p - dictionary that contains the relevant system parameters
        m - index of momentum vector km = pi*m/(N+1) (momentum quantisation
            for open boudary conditions), m = 1,...,N
        species - spin species in question (1 or 2)
        state_table - list of all state_numbers that belong to the relevant
            Hilbertspace
    Returns:
        nk - n_{km,species} matrix on the relevant Hilbertspace
    """
    dim = len(state_table)
    row = []
    col = []
    data = []

    nk = sparse.csr_matrix((data, (row, col)), shape=(dim, dim), dtype=complex)
    k = (np.pi * m) / (p['N'] + 1)

    for i in range(1, p['N'] + 1):
        # for j in range(1,p['N']+1):
        #   # Fourier Transform for open boundary conditions goes as ~ sqrt(2/(N+1)) * sin(k*i)
        #     nk += (2./(p['N']+1)) * CdagC_matrix(p, i-1, species, j-1, species, state_table) * np.sin(k*i) * np.sin(k*j)

        # more efficient as it takes only ~ 1/2 the number of summations (make use of the fact that i<j matrices are daggered matrices for i>j!)
        for j in range(i, p['N'] + 1):
            if i == j:
                nk += (2. / (p['N'] + 1)) * CdagC_matrix(p, i - 1, species, j - 1, species, state_table) * np.sin(k * i) * np.sin(k * j)
            else:
                amplitude = CdagC_matrix(
                    p, i - 1, species, j - 1, species, state_table)
                nk += (2. / (p['N'] + 1)) * (amplitude + np.transpose(
                    np.conjugate(amplitude))) * np.sin(k * i) * np.sin(k * j)
    return nk


def make_Hamiltonian(p, state_table):
    """
    Generates full Hamiltonian on the relevant sub-Hilbertspace
    Args:
        p - dictionary that contains the relevant system parameters
        state_table - list of all state_numbers that belong to the relevant
            Hilbertspace
    Returns:
        H - Hamiltonian matrix on the relevant Hilbertspace
    """
    # dim = len(state_table)
    # row = []
    # col = []
    # data = []

    # H = sparse.csr_matrix((data, (row, col)), shape=(dim, dim), dtype=complex)
    H = hopping_matrix(p, state_table)
    H += onsite_interaction_matrix(p, state_table)
    H += chemical_potential_matrix(p, state_table)

    return H


def make_trotter_Hamiltonian(p, state_table, order=1):
    """
    Generates list of Hamiltonians which can be executed sequentially to
    advance one trotter timestep
    Args:
        p - dictionary that contains the relevant system parameters
        state_table - list of all state_numbers that belong to the relevant
            Hilbertspace
    Returns:
        H_list - list of Hamiltonians to be applied sequeationally
            on the relevant Hilbertspace
    """
    H_list = []
    H_list.append(hopping_matrix(p, state_table, bond_type=1))
    H_list.append(hopping_matrix(p, state_table, bond_type=2))
    H_list.append(hopping_matrix(p, state_table, bond_type=3))
    H_list.append(hopping_matrix(p, state_table, bond_type=4))
    H_list.append(onsite_interaction_matrix(p, state_table))
    H_list.append(chemical_potential_matrix(p, state_table))
    return H_list


def calculate_gs(p):
    """
    calculates groundstate of full Hamiltonian on the relevant sub-Hilbertspace
    Args:
        p - dictionary that contains the relevant system parameters for the
            groundstate search
    Returns:
        E0 - groundstate energy
        gs - groundstate vector on the relevant Hilbertspace subspace
        state_table - list of all state_numbers that belong to the relevant
            Hilbertspace
    """
    state_table = generate_state_table(p)
    H = make_Hamiltonian(p, state_table)
    w, v = scipy.sparse.linalg.eigsh(H, k=1, which='SA')

    return w[0], v[:, 0], state_table


def expct_val(Op, psi):
    """
    compute expecation value of operator 'Op' with state 'psi'
    Args:
        Op - operator corresponding to observable to be measured
        psi - state-vector (on sub-Hilbertspace)
    Returns:
        <psi| Op |psi>
    """
    return (psi.conj().T).dot(Op.dot(psi))


def evolve(p, state_table, state, kind="list", correlation_measurement=False,
           trotterised=False):
    """
    evolve 'state' under parameters defined in dictionary 'p'
    Args:
        p - dictionary that contains the relevant system parameters
            for time-evolution
        state - fermion configuration OR state-vector on the relevant
            Hilbertspace
        kind - which kind of state is passed to the function: kind=list
            (default) spin-configuration (productstate) OR kind="ket" arbitrary
            vector in Hilbert-subspace
                OR kind="int" the unique state id in the state_table
        correlation_measurement - boolean (default=False), should
            single-particle correlations (and from this momentum
            distribution n_k) be measured?
    Returns:
        sim - dictionary with the relevant density measurements: N1, N2, N12
        state_table - list of all state_numbers that belong to the relevant
            Hilbertspace
    """
    if p['L'] != 1 and p['W'] != 1 and correlation_measurement:
        warnings.warn('No correlation for lattices, tough break pal')
        correlation_measurement = False
    if kind == "ket":
        psi0 = state
    elif kind == "list":
        # if we parsed a product state, construct ket by identifying the
        # corresponding number of the basis state and putting a 1 into the ket
        psi0 = np.zeros((len(state_table), 1), dtype=complex)
        psi0[state_table.index(state_to_int(p, state))] = 1.
    elif kind == "int":
        psi0 = np.zeros((len(state_table), 1), dtype=complex)
        psi0[state_table.index(state)] = 1.

    time = np.linspace(p['t_initial'], p['t_final'],
                       int(p['t_final'] / p['dt'] + 1))

    # make dictionary with measurement operators
    meas = {}
    for i in range(int(p['N'])):
        for species in range(1, 3):
            meas['N' + str(species) + ' Site ' + str(i + 1)
                 ] = ni_matrix(p, i, species, state_table)

        meas['N12 Site ' + str(i + 1)] = meas['N1 Site ' + str(i + 1)].dot(meas['N2 Site ' + str(i + 1)])

        if correlation_measurement:
            for m in range(1, p['N'] + 1, 1):
                for species in range(1, 3, 1):
                    meas['n_k%i_%i' % (m, species)] = nk_matrix(
                        p, m, species, state_table)

    sim = {}
    sim['Time'] = time
    for key in meas.keys():
        # don't include the correlation measurements here,
        # later by hand for convenience (see below)
        if key[:3] != 'n_k':
            sim['Re(' + key + ')'] = np.zeros(np.shape(time))
    if correlation_measurement:
        sim['k'] = np.pi * \
            np.linspace(start=1, stop=p['N'], num=p['N']) / (p['N'] + 1)
        sim['nk1'] = np.zeros(shape=(len(sim['k']), len(time)))
        sim['nk2'] = np.zeros(shape=(len(sim['k']), len(time)))

    if trotterised:
        H_list = make_trotter_Hamiltonian(p, state_table)
    else:
        H_list = [make_Hamiltonian(p, state_table)]

    # construct time-evolution operators for a single time-step
    U_list = [LA.expm(-1.j * H.tocsc() * p['dt']) for H in H_list]

    # Time Evolution
    for i in range(len(time)):
        # define initial (t=0) state
        if i == 0:
            psi = psi0

        # measurements
        for operator in meas.keys():
            expct = expct_val(meas[operator], psi)  # [0][0]

            # if we do a correlation measurement (nk), put it into the nk1
            # and nk2 matrix by hand
            if operator[:3] == 'n_k':
                kindex = int(operator[3])

                # if np.imag(expct) == 0.0:
                if np.imag(expct) < 1e-12:
                    if operator[-1] == '1':
                        sim['nk1'][kindex - 1, i] = np.real(expct)
                    elif operator[-1] == '2':
                        sim['nk2'][kindex - 1, i] = np.real(expct)
                else:
                    print("Imaginary Measurement %s" % (operator))
            else:
                # if np.imag(expct) == 0.0:
                if np.imag(expct) < 1e-12:
                    sim['Re(' + operator + ')'][i] = np.real(expct)
                else:
                    print("Imaginary Measurement %s" % (operator))

        # apply U to current state psi to get psi(t+dt) = U * psi(t)
        for U in U_list:
            psi = U.dot(psi)

    return sim, state_table


def calculate_average_errors(p, sim_exact, sim):
    sim_errors = {'Time': sim_exact['Time']}
    n_errors = np.empty((3, p['N'], len(sim_errors['Time'])))
    for s in range(p['N']):
        n_errors[0, s] = abs(sim[f'Re(N1 Site {s + 1})'] - sim_exact[f'Re(N1 Site {s + 1})'])
        n_errors[1, s] = abs(sim[f'Re(N2 Site {s + 1})'] - sim_exact[f'Re(N2 Site {s + 1})'])
        n_errors[2, s] = abs(sim[f'Re(N12 Site {s + 1})'] - sim_exact[f'Re(N12 Site {s + 1})'])
    sim_errors['Re(N1)'] = np.mean(n_errors[0], axis=0)
    sim_errors['Re(N2)'] = np.mean(n_errors[1], axis=0)
    sim_errors['Re(N12)'] = np.mean(n_errors[2], axis=0)
    if 'k' in sim_exact and 'k' in sim:
        sim_errors['k'] = sim_exact['k']
        sim_errors['nk1'] = abs(sim['nk1'] -  sim_exact['nk1'])
        sim_errors['nk2'] = abs(sim['nk2'] -  sim_exact['nk2'])
    return sim_errors
