import numpy as np, h5py
from EasyED.ensembles import GrandcanonicalEnsemble
from EasyED.hamiltonians import Hubbard
from reduceddensitymatrix import ReducedDensityMatrix as RDM


def get_mutual_information(rhoa, rhob, rhoab):
    ra = rhoa.get_entanglement_entropy()
    rb = rhob.get_entanglement_entropy()
    rab = rhoab.get_entanglement_entropy()
    return ra+rb-rab

def get_u_mu24(tnnn):
    with h5py.File('intersectionstprime.h5', 'r') as data:
        i_tnnn = np.where(data['tnnns'][:] == tnnn)[0]
        phaseboundary = data['intersections'][i_tnnn, 1]
    return phaseboundary

beta = 30
t = -1
tp = .3
hop = [[0,t,t,tp],[t,0,tp,t],[t,tp,0,t],[tp,t,t,0]]
umu = get_u_mu24(0.3)
results = []

for u, mu in zip(umu[0,:,0], umu[0,:,1]):
    print u, mu
    hamiltonian = Hubbard(hop, u) # basis: [['up', 'dn'], [0, 1, 2, 3]]
    plaq = GrandcanonicalEnsemble(hamiltonian, beta, mu)
    plaq.calcEigensystem()
    energies = plaq.hamiltonian.eigenEnergies
    eigenstates = plaq.hamiltonian.eigenStates.toarray()
    rho_0 = RDM(beta, eigenstates, energies, 2, 6, [(3,4),(2,3),(1,2)], range(9))
    rho_0.calculate()
    rho_01 = RDM(beta, eigenstates, energies, 4, 4, [(3,4),(2,3),(1,2),(4,5),(3,4)], range(9))
    rho_01.calculate()
    rho_03 = RDM(beta, eigenstates, energies, 4, 4, [(3,4),(2,3),(1,2),(3,4),(2,3),(6,7),(5,6),(4,5),
                                                 (3,4)], range(9))
    rho_03.calculate()
    i01 = get_mutual_information(rho_0, rho_0, rho_01)
    i03 = get_mutual_information(rho_0, rho_0, rho_03)
    results.append([beta, t, tp, u, mu, i01, i03])
    np.save('mutual_information.npy', np.array(results))
print 'mutual_information.npy ready'
