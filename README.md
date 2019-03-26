# reduceddensitymatrix
The class ReducedDensityMatrix depends only on numpy. Thus eigenstates and eigenvectors are arrays. n_a, n_b define subblock sizes, the reduced density matrix will act on the subspace of A and B is traced out. The eigenvectors have a Fockbasis (order of single particle states) that has to be permuted pairwisely into A and B subblocks. For performance only some eigenstates can be taken into account (state_selection).
mutual_information.py and plaquette.py require EasyED, mutual_information_plt.py requires mpltotex
