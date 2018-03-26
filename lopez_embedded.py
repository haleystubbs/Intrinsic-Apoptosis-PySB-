"""
Model M1a: Extrinsic apoptosis model with expanded, "embedded together"
model of MOMP.
"""

from pysb import *
from pysb.util import alias_model_components

from earm import shared
from earm import lopez_modules
from earm import albeck_modules

"""_________________________________________________________________________________________________________________________________________"""

Model()

def declare_p53_upregulated_monomers()
    # declare some monomers that p53 upregulates transcription for that are not already included in EARM
    Monomer('Puma', ['bf', 'state'], {'state': ['M', 'C']})     #binds all prosurvival Bcl-2 members, along with Bim, and tBid
    Monomer('TRAILR2')
    Monomer('PIDD')
    Monomer('p53', ['bf', 'state'], {'state': ['N', 'C']})      # has a single binding site and nuclear and cytoplasmic states

    alias_model_components()

    Annotation(Puma, 'http://www.uniprot.org/uniprot/Q9BXH1')
    Annotation(PIDD, 'http://www.uniprot.org/uniprot/Q9HB75')
    Annotation(p53AIP1, 'http://www.uniprot.org/uniprot/A7LAQ2')

def p53_upregulated_parameters()
    # add parameters/rates for all the proteins that p53 transcriptionally upregulates
    Parameter('upregulate_Puma_ks')
    Parameter('upregulate_Bax_ks')
    Parameter('upregulate_Noxa_ks')
    Parameter('upregulate_Bak_ks')
    Parameter('upregulate_Bid_ks')
    Parameter('upregulate_Apaf_ks')
    alias_model_components()

def Puma_binds_Bcl2_members()
    #DNA damage leads to Bim and Bmf phosphorylation, these two Bcl-2 family members are released from motor proteins that usually sequester them
    Equilibrate(Puma(bf=None, state='C'), Puma(bf=None, state='M'))
    #BclxL binds Puma in the cytosol and the mitochondria, blobking Puma'a ability to activate Bax and Bak
    Rule('Puma_bind_BclxL_C', Puma(bf=None, state-'C') + BclxL(bf=None, state='C') | Puma(bf=1, state='C') % BclxL(bf=1, state='C'))
    Rule('Puma_bind_BclxL_M', Puma(bf=None, state-'M') + BclxL(bf=None, state='M') | Puma(bf=1, state='M') % BclxL(bf=1, state='M'))
    #



def p53_upregulates_expression()
    #p53 acts as a transcfription factor to upregulate expression of the following proteins
    #create observable for the concentration of P53 and multiply the upregulated synthesis rates  by the conc of P53
    Observable('Obs_p53', p53(bf=None, state='C'))
    #equilirate p53 between the nucleus and the cytoplasm
    Equilibrate(p53(bf=None, state='N'), p53(bf=None, state='C'))
    Rule('p53_upregulate_Puma', None >> Puma(bf=None, state='C'), upregulate_Puma_ks*Obs_p53)
    Rule('p53_upregulate_Bax', None >> Bax(bf=None, s1=None, s2=None, state='C'), upregulate_Bax_ks*Obs_p53)
    Rule('p53_upregulate_Noxa', None >> Noxa(bf=None, state='C'), upregulate_Noxa_ks*Obs_p53)
    Rule('p53_upregulate_Bak', None >> Bak(bf=None, s1=None, s2=None, state='M'), upregulate_Bak_ks*Obs_p53)
    Rule('p53_upregulate_Bid', None >> Bid(bf=None, state='U'), upregulate_Bid_ks*Obs_p53)
    Rule('p53_upregulate_Apaf', None >> Apaf(bf=None, state='I'), upregulate_Apaf_ks*Obs_p53)
    #BclxL binds p53 in the cytosol and keeps p53 from upregulated expression of the above proteins
    Rule('BclxL_bind_p53', BclxL(bf=None, state='C') + p53(bf=None, state='C') | BclxL(bf=1, state='C') % p53(bf=1, state='C'))

"""_________________________________________________________________________________________________________________________________________"""

# Declare monomers
albeck_modules.ligand_to_c8_monomers()
lopez_modules.momp_monomers()
albeck_modules.apaf1_to_parp_monomers()

# Generate the upstream and downstream sections
albeck_modules.rec_to_bid()
albeck_modules.pore_to_parp()

# The specific MOMP model to use
lopez_modules.embedded()

# Declare observables
shared.observables()

print(model.monomers)