

from pysb import *
from pysb.util import alias_model_components

from earm import shared
from earm import lopez_modules
from earm import albeck_modules

"""_________________________________________________________________________________________________________________________________________"""

Model()
# Bim, states:
Monomer('Bim', )
#Puma, states: Mitochondrial and Cytoplasmic
Monomer('Puma', ['bf', 'state'], {'state': ['M', 'C']})
# p53, states: Nuclear and Cytoplasmic
Monomer('p53', ['bf', 'state'], {'state': ['N', 'C']})
# Bax, states: Cytoplasmic, Mitochondrial, Active
# sites 's1' and 's2' are used for pore formation
Monomer('Bax', ['bf', 's1', 's2', 'state'], {'state':['C', 'M', 'A']})
# Bak, states: inactive and Mitochondrial, Active (and mitochondrial)
# Bak is constitutively inserted into the outer mitochondrial membrane
# sites 's1' and 's2' are used for pore formation
Monomer('Bak', ['bf', 's1', 's2', 'state'], {'state':['M', 'A']})
# Bid, states: Untruncated, Truncated, truncated and Mitochondrial
Monomer('Bid', ['bf', 'state'], {'state':['U', 'T', 'M']})
# **Anti-Apoptotics**
Monomer('Bcl2', ['bf'])
Monomer('BclxL', ['bf', 'state'], {'state': ['C', 'M']})
Monomer('Mcl1', ['bf', 'state'], {'state': ['M', 'C']})


def declare_initial_conditions():
    """Declare initial conditions for Bcl-2 family proteins, Cyto c, and Smac.
    """
    Parameter('Bid_0'   , 4.0e4) # Bid
    Parameter('BclxL_0' , 2.0e4) # cytosolic BclxL
    Parameter('Mcl1_0'  , 2.0e4) # Mitochondrial Mcl1
    Parameter('Bcl2_0'  , 2.0e4) # Mitochondrial Bcl2
    Parameter('Bax_0'   , 0.8e5) # Bax
    Parameter('Bak_0'   , 0.2e5) # Bak

    alias_model_components()

    Initial(Bid(bf=None, state='U'), Bid_0)
    Initial(Bax(bf=None, s1=None, s2=None, state='C'), Bax_0)
    Initial(Bak(bf=None, s1=None, s2=None, state='M'), Bak_0)
    Initial(Bcl2(bf=None), Bcl2_0)
    Initial(BclxL (bf=None, state='C'), BclxL_0)
    Initial(Mcl1(bf=None, state='M'), Mcl1_0)


def effectors_bind_anti_apoptotics():
    """Binding of Bax and Bak to Bcl2, BclxL, and Mcl1.

    Affinities of Bak for Bcl-xL and Mcl-1 are taken from Willis et al.

    Preferential affinity of Bax for Bcl-2 and Bcl-xL were taken from Zhai et
    al.  Bax:Bcl2 and Bax:Bcl-xL affinities were given order of magnitude
    estimates of 10nM.

    See comments on units for :py:func:`tBid_binds_all_anti_apoptotics`.

    Willis, S. N., Chen, L., Dewson, G., Wei, A., Naik, E., Fletcher, J. I.,
    Adams, J. M., et al. (2005). Proapoptotic Bak is sequestered by Mcl-1 and
    Bcl-xL, but not Bcl-2, until displaced by BH3-only proteins. Genes &
    Development, 19(11), 1294-1305. `doi:10.1101/gad.1304105`

    Zhai, D., Jin, C., Huang, Z., Satterthwait, A. C., & Reed, J. C. (2008).
    Differential regulation of Bax and Bak by anti-apoptotic Bcl-2 family
    proteins Bcl-B and Mcl-1. The Journal of biological chemistry, 283(15),
    9580-9586.  `doi:10.1074/jbc.M708426200`
    """

    bind_table([[                            Bcl2,  BclxL(state='M'),         Mcl1],
                [Bax(active_monomer), 10e-9*N_A*V,       10e-9*N_A*V,         None],
                [Bak(active_monomer),        None,       50e-9*N_A*V,  10e-9*N_A*V]],
               kf=1e6/(N_A*V))




