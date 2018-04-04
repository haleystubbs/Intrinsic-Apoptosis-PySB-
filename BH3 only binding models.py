
"""Hit and Run Models of Bax and Bak Activation by PUMA, BIM, tBid, and p53"""

# PUMA, tBid, and BIM directly interact with Bax and Bak to catalyze oligomerizatoin/pore formation
# numerous reports have suggested that these three proteins as well as p53 transiently and dynamically bind
# Bax and Bak that leads to a stepwise series of conformational changes that lead to activation so that
# Bax and Bak can oligomerize and form pores for MOMP
# Here, several different mechanisms for each of the proteins are written in order to investigate the nature
# of these direct interactions

# References:
# Kim, H. et al. Stepwise Activation of BAX and BAK by tBID, BIM, and PUMA Initiates Mitochondrial Apoptosis. Mol. Cell 36, 487–499 (2009).


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



