
"""Monomer Abbreviations"""

# Apaf1 = Apoptotic protease activating factor 1, conc here is representative of activated Apaf1 in apoptosomes
# C9 = caspase 9, 2 binding sites for Apaf1 and XIAP, can be pro or cleaved by Apaf1
# XIAP = X linked inhibitor of apoptosis protein
# C3 = caspase 3, 2 binding sites for C9 and XIAP, can be pro or cleaved by C9



# import the pysb module and all its methods and functions
from pysb import *
import matplotlib.pyplot as plt
from pysb.bng import *
from pysb.simulator import ScipyOdeSimulator
import pylab as pl

# Instantiate a model
Model()

"""Declare Monomers"""
Monomer('Apaf1', ['bC9'])
# create monomer Apaf1 with single binding site for CC
Monomer('C9', ['bApaf1', 'bXIAP', 'bC3', 'm'], {'m': {['p', 'c']}})
# create monomer C9 with a binding and cleavage site, the cleavage site can be pro (before cleavage, or cleaved (after cleavage)
Monomer('XIAP', ['bC9', 'bC3'])
# create monomer XIAP with binding site a to bind to C9 and binding site b to bind C3
Monomer('C3', ['bC9','bXIAP','m'], {'m':['p','c']})
# create monomer C3 with binding site a for C9, b binding site for XIAP, and c cleavage site with states pro or cleaved


"""Parameters"""
# Parameters of the forward, reverse and catalysis rate constants listed in order of appearance in rules below
# kc = catalytic rate
# kf = forward reaction rate
# kr = reverse reaction rate
# ks = rate of synthesis
# kd = rate of degradation
Parameter('Apaf1_bind_C9p_kf', 2e-03)                     # nM^-1 s^-1
Parameter('Apaf1_bind_C9p_kr', 0.1)                       # s^-1
Parameter('Apaf1_bind_C9c_kf', 2e-03)                     # nM^-1 s^-1
Parameter('Apaf1_bind_C9c_kr', 0.1)                       # s^-1
Parameter('Apaf1_C9c_bind_XIAP_kf', 1e-03)                # nM^-1 s^-1
Parameter('Apaf1_C9c_bind_XIAP_kr', 1e-03)                # s^-1
Parameter('Apaf1_C9p_bind_XIAP_kf', 1e-03)                # nM^-1 s^-1
Parameter('Apaf1_C9p_bind_XIAP_kr', 1e-03)                # s^-1
Parameter('C9p_XIAP_bind_Apaf1_kf', 2e-03)                # nM^-1 s^-1
Parameter('C9p_XIAP_bind_Apaf1_kr', 0.1)                  # s^-1
Parameter('C9c_XIAP_bind_Apaf1_kf', 2e-03)                # nM^-1 s^-1
Parameter('C9c_XIAP_bind_Apaf1_kr', 0.1)                  # s^-1
Parameter('C9p_bind_XIAP_kf', 1e-03)                      # nM^-1 s^-1
Parameter('C9p_bind_XIAP_kr', 1e-03)                      # s^-1
Parameter('C9c_bind_XIAP_kf', 1e-03)                      # nM^-1 s^-1
Parameter('C9c_bind_XIAP_kr', 1e-03)                      # s^-1
Parameter('C3c_bind_XIAP_kf', 3e-03)                      # nM^-1 s^-1
Parameter('C3c_bind_XIAP_kr', 1e-03)                      # s^-1
Parameter('C9p_cleave_C3p_kc', 5e-06)                     # nM^-1 s^-1
Parameter('C9c_cleave_C3p_kc', 5e-05)                     # nM^-1 s^-1
Parameter('C9p_Apaf1_cleave_C3p_kc', 3.5e-04)             # nM^-1 s^-1
Parameter('C9c_Apaf1_cleave_C3p_kc', 3.5e-03)             # nM^-1 s^-1
Parameter('C3c_cleave_C9p_kc', 2e-04)                     # nM^-1 s^-1
Parameter('C3c_cleave_C9p_Apaf1_kc', 2e-04)               # nM^-1 s^-1
Parameter('degrade_Apaf1_kd', 1e-03)                      # s^-1
Parameter('degrade_C9p_kd', 1e-03)                        # s^-1
Parameter('synthesize_XIAP_ks', 0.04)                     # nM s^-1
Parameter('degrade_C3p_kd', 1e-03)                        # s^-1
Parameter('synthesize_Apaf1_ks', 0.02)                    # nM s^-1
Parameter('synthesize_C9p_ks', 0.02)                      # nM s^-1
Parameter('degrade_C9c_kd', 1e-03)                        # s^-1
Parameter('synthesize_C3p_ks', 0.2)                       # nM s^-1
Parameter('degrade_C3c_kd', 1e-03)                        # s^-1
Parameter('degrade_XIAP_kd', 1e-03)                       # s^-1





"""Rules for activated Apaf1 equilibrium binding"""
Rule('Apaf1_bind_C9p', Apaf1(bC9=None) + C9(bApaf1=None, m='p') | Apaf1(bC9=1) % C9(bApaf1=1, m='p'), Apaf1_bind_C9p_kf, Apaf1_bind_C9p_kr)
# Activated Apaf1 reversibly binds procaspase 9
Rule('Apaf1_bind_C9c', Apaf1(bC9=None) + C9(bApaf1=None, m='c') | Apaf1(bC9=1) % C9(bApaf1=1, m='p'), Apaf1_bind_C9c_kf, Apaf1_bind_C9c_kr)
# Activated Apaf1 reversibly binds caspase 9
Rule('Apaf1_C9c_bind_XIAP', Apaf1(bC9=1) % C9(bApaf1=1, bXIAP=None, m='c') + XIAP(bC9=None) | Apaf1(bC9=1) % C9(bApaf1=1, bXIAP=1, m='c') % XIAP(bC9=1), Apaf1_C9c_bind_XIAP_kf, Apaf1_C9c_bind_XIAP_kr)
# the Apaf1 caspase 9 complex binds XIAP
Rule('Apaf1_C9p_bind_XIAP', Apaf1(bC9=1) % C9(bApaf1=1, bXIAP=None, m='p') + XIAP(bC9=None) | Apaf1(bC9=1) % C9(bApaf1=1, bXIAP=1, m='p') % XIAP(bC9=1), Apaf1_C9p_bind_XIAP_kf, Apaf1_C9p_bind_XIAP_kr)
# the Apaf1 procaspase 9 complex binds XIAP
Rule('C9p_XIAP_bind_Apaf1', Apaf1(bC9=None) + C9(bApaf1=None, bXIAP=1, m='p') % XIAP(bC9=1) | Apaf1(bC9=1) % C9(bApaf1=1, bXIAP=1, m='p') % XIAP(bC9=1), C9p_XIAP_bind_Apaf1_kf, C9p_XIAP_bind_Apaf1_kr)
# the procaspase 9 XIAP complex binds Afap1
Rule('C9c_XIAP_bind_Apaf1', Apaf1(bC9=None) + C9(bApaf1=None, bXIAP=1, m='c') % XIAP(bC9=1) | Apaf1(bC9=1) % C9(bApaf1=1, bXIAP=1, m='c') % XIAP(bC9=1), C9c_XIAP_bind_Apaf1_kf, C9c_XIAP_bind_Apaf1_kr)
# the caspase 9 XIAP complex binds Apaf1


"""Rules for C9 & C3 equilibrium binding"""
Rule('C9p_bind_XIAP', C9(bXIAP=None, m='p') + XIAP(bC9=None) | C9(bXIAP=1, m='p') % XIAP(bC9=1), C9p_bind_XIAP_kf, C9p_bind_XIAP_kr)
# Procaspase 9 binds XIAP
Rule('C9c_bind_XIAP', C9(bXIAP=None, m='c') + XIAP(bC9=None) | C9(bXIAP=1, m='c') % XIAP(bC9=1), C9c_bind_XIAP_kf, C9c_bind_XIAP_kr)
# Caspase 9 binds XIAP
Rule('C3c_bind_XIAP', C3(bXIAP=None, m='c') + XIAP(bC3=None) | C3(bXIAP=1, m='c') % XIAP(bC3=1), C3c_bind_XIAP_kf, C3c_bind_XIAP_kr)
# Caspase 3 binds XIAP


"""Rules for Catalysis of Procaspase 3 Cleavage"""
Rule('C9p_cleave_C3p', C3(m='p') + C9(m='p') >> C3(m='c') + C9(m='p'), C9p_cleave_C3p_kc)
# procaspase 9 catalyses cleavage of procaspase 3 to caspase 3
Rule('C9c_cleave_C3p', C3(m='p') + C9(m='c') >> C3(m='c') + C9(m='c'), C9c_cleave_C3p_kc)
# Caspase 9 catalyses cleavage of procaspase 3 to caspase 3
Rule('C9p_Apaf1_cleave_C3p', C3(m='p') + Apaf1(bC9=1) % C9(bApaf1=1, m='p') >> C3(m='c') + Apaf1(bC9=1) % C9(bApaf1=1, m='p'), C9p_Apaf1_cleave_C3p_kc)
# Complex of activated Apaf1 and procaspase 9 catalyze cleavage of procaspase 3 to caspase 3
Rule('C9c_Apaf1_cleave_C3p', C3(m='p') + Apaf1(bC9=1) % C9(bApaf1=1, m='c') >> C3(m='c') + Apaf1(bC9=1) % C9(bApaf1=1, m='c'), C9c_Apaf1_cleave_C3p_kc)
# Complex of activated Apaf1 and caspase 9 catalyze cleavage of procaspase 3 to caspase 3


"""Rules for Catalysis of Procaspase 9 Cleavage"""
Rule('C3c_cleave_C9p', C9(m='p') + C3(m='c') >> C9(m='c') + C3(m='c'), C3c_cleave_C9p_kc)
# Caspase 3 catalyzes the cleavage of procaspase 9 to caspase 9
Rule('C3c_cleave_C9p_Apaf1', C3(m='c') + C9(bApaf1=1, m='p') % Apaf1(bC9=1) >> C3(m='c') + C9(bApaf1=1, m='c') % Apaf1(bC9=1), C3c_cleave_C9p_Apaf1_kc)
# Caspase 3 catalyzes the cleavage of procaspase 9 bound to Apaf1 to caspase 9 bound to Apaf1


"""Rules for the Synthesis and Degradation of all single protein species"""
Rule('synthesize_Apaf1', None >> Apaf1(), synthesize_Apaf1_ks)
# Make activated Apaf1
Rule('synthesize_C9p', None >> C9(m='p'), synthesize_C9p_ks)
# Make procaspase 9
Rule('synthesize_XIAP', None >> XIAP(), synthesize_XIAP_ks)
# Make XIAP
Rule('synthesize_C3p', None >> C3 (m='p'), synthesize_C3p_ks)
# Make procaspase 3
Rule('degrade_Apaf1', Apaf19() >> None, degrade_Apaf1_kd)
# Degrade Apaf1
Rule('degrade_C9p', C9(m='p') >> None, degrade_C9p_kd)
# Degrade procaspase 9
Rule('degrade_C9c', C9(m='c') >> None, degrade_C9c_kd)
# Degrade caspase 9
Rule('degrade_C3p', C3(m='p') >> None, degrade_C3p_kd)
# Degrade procaspase 3
Rule('degrade_C3c', C3(m='c') >> None, degrade_C3c_kd)
# Degrade caspase 3
Rule('degrade_XIAP', XIAP() >> None, degrade_XIAP_kd)
# Degrade XIAP


"""Rules for Degradation of Complexes"""


