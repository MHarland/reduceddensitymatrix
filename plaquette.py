import numpy as np
from EasyED.ensembles import GrandcanonicalEnsemble
from EasyED.hamiltonians import Hubbard, Hubbard2
from EasyED.util import report
from reduceddensitymatrix import ReducedDensityMatrix as RDM


beta = 10
t = -1
tp = 0
u = 3
mu = u*.5

hop = [[0,t,t,tp],[t,0,tp,t],[t,tp,0,t],[tp,t,t,0]]
hamiltonian = Hubbard2(hop, u)
plaq = GrandcanonicalEnsemble(hamiltonian, beta, mu)
plaq.calcEigensystem()
e0 = plaq.hamiltonian.getGroundStateEnergy()
energies = plaq.hamiltonian.eigenEnergies
eigenstates = plaq.hamiltonian.eigenStates.toarray()
rho_dim = RDM(beta, eigenstates, energies, 4, 4, [], [0])
rho_dim.calculate()
print rho_dim.get_entanglement_spectrum()
print rho_dim.get_entanglement_entropy()

hamiltonian = Hubbard(hop, u)
print hamiltonian.singleParticleBasis
plaq = GrandcanonicalEnsemble(hamiltonian, beta, mu)
plaq.calcEigensystem()
e0 = plaq.hamiltonian.getGroundStateEnergy()
energies = plaq.hamiltonian.eigenEnergies
eigenstates = plaq.hamiltonian.eigenStates.toarray()
rho_0 = RDM(beta, eigenstates, energies, 4, 4, [(3,4),(2,3),(1,2),(4,5),(3,4)], [0])
rho_0.calculate()
print rho_0.get_entanglement_spectrum()
print rho_0.get_entanglement_entropy()

hamiltonian = Hubbard(hop, u)
plaq = GrandcanonicalEnsemble(hamiltonian, beta, mu)
plaq.calcEigensystem()
e0 = plaq.hamiltonian.getGroundStateEnergy()
energies = plaq.hamiltonian.eigenEnergies
eigenstates = plaq.hamiltonian.eigenStates.toarray()
rho_0 = RDM(beta, eigenstates, energies, 4, 4, [(3,4),(2,3),(1,2),(2,3),(5,6),(4,5),(3,4)], [0])
rho_0.calculate()
print rho_0.get_entanglement_spectrum()
print rho_0.get_entanglement_entropy()


"""
hop = [[0,t,tp,t],[t,0,t,tp],[tp,t,0,t],[t,tp,t,0]]
hamiltonian = Hubbard(hop, u)
plaq = GrandcanonicalEnsemble(hamiltonian, beta, mu)
plaq.calcEigensystem()
e0 = plaq.hamiltonian.getGroundStateEnergy()
energies = plaq.hamiltonian.eigenEnergies
eigenstates = plaq.hamiltonian.eigenStates.toarray()
rho_0 = RDM(eigenstates, energies, [0,1,2,3], [4,5,6,7], beta)
rho_0.calculate()
print rho_0.get_entanglement_spectrum()
print rho_0.get_entanglement_entropy()

hamiltonian = Hubbard(hop, u)
plaq = GrandcanonicalEnsemble(hamiltonian, beta, mu)
plaq.calcEigensystem()
e0 = plaq.hamiltonian.getGroundStateEnergy()
energies = plaq.hamiltonian.eigenEnergies
eigenstates = plaq.hamiltonian.eigenStates.toarray()
rho_1 = RDM(eigenstates, energies, [0,1,2,3], [4,5,6,7], beta)
rho_1.calculate()
print rho_1.get_entanglement_spectrum()
print rho_1.get_entanglement_entropy()
"""
