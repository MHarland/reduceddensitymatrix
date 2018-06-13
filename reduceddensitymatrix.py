import numpy as np, itertools as itt


class ReducedDensityMatrix:
    """
    eigenstates is a matrix, with columns as eigenstates
    subsystem_a/b are indices of the subsystems in occupation number representation
    B will be traced out
    state_selection by energy
    """
    def __init__(self, beta, eigenstates, eigenenergies, n_a, n_b, permutations = [], state_selection = None, verbose = False):
        self.beta = beta
        self.eigenstates = np.array(eigenstates)
        self.eigenenergies = np.array(eigenenergies)
        self.n_full_fock = len(eigenenergies)
        self.i_full_fock = range(self.n_full_fock)
        
        self.n_a_sp = n_a
        self.n_a_fock = 2**n_a
        self.i_a_sp = range(n_a)
        self.i_a_fock = range(self.n_a_fock)
        
        self.n_b_sp = n_b
        self.n_b_fock = 2**n_b
        self.i_b_sp = range(n_b)
        self.i_b_fock = range(self.n_b_fock)

        self.permutations = permutations
        
        self.n_ab_sp = self.n_a_sp + self.n_b_sp
        if state_selection is None:
            self.state_selection = self.i_full_fock
        else:
            state_selection = np.array(state_selection)
            self.state_selection = np.argsort(eigenenergies)[state_selection]
        self.n_selection = len(self.state_selection)
        self.i_state_selection = range(self.n_selection)
        self.rho_full = np.zeros([self.n_selection])
        self.c = np.zeros([self.n_selection, self.n_a_fock, self.n_b_fock])
        self.rho_a = np.zeros([self.n_a_fock, self.n_a_fock])
        self.verbose = verbose

    def _calc_rho_full(self):
        if self.verbose: print 'calculating rho_full...'
        expbetae = np.exp(-self.beta * self.eigenenergies[self.state_selection])
        self.rho_full = expbetae / np.sum(expbetae)

    def _get_subspace_string(self, full_binary_string, subspace_indices):
        result = ''
        for a in subspace_indices:
            result += full_binary_string[a]
        return result

    def _calc_c(self):
        if self.verbose: print 'calculating c...'
        for n, a, b, m in itt.product(self.i_state_selection, self.i_a_fock, self.i_b_fock, self.i_full_fock):
            bm = list(np.binary_repr(m, self.n_ab_sp))
            ba = list(np.binary_repr(a, self.n_a_sp))
            bb = list(np.binary_repr(b, self.n_b_sp))
            sign = 1
            for p1, p2 in self.permutations:
                s1, s2 = bm[p1], bm[p2]
                if s1 == '1' and s2 == '1':
                    sign *= -1
                bm[p1] = s2
                bm[p2] = s1
            if ba == bm[:self.n_a_sp] and bb == bm[self.n_a_sp:]:
                self.c[n, a, b] += sign * self.eigenstates[m, self.state_selection[n]]

    def _calc_rho_a(self):
        if self.verbose: print 'calculating rho_a...'
        for a, aa, b, n in itt.product(self.i_a_fock, self.i_a_fock, self.i_b_fock, self.i_state_selection):
            self.rho_a[a, aa] += self.rho_full[n] * self.c[n,a,b] * self.c[n,aa,b].conjugate()

    def calculate(self):
        self._calc_rho_full()
        self._calc_c()
        self._calc_rho_a()

    def __getitem__(self, index):
        return self.rho_a[index]

    def get_entanglement_spectrum(self):
        return np.linalg.eigvalsh(self.rho_a)

    def get_entanglement_entropy(self):
        es = self.get_entanglement_spectrum()
        nozeros = []
        for e in es:
            if not np.allclose(e, 0):
                nozeros.append(e)
        nozeros = np.array(nozeros)
        return -np.sum(nozeros*np.log(nozeros))
