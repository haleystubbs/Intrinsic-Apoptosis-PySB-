
"""Monomer Abbreviations"""

# Apaf1 = Apoptotic protease activating factor 1, conc here is representative of activated Apaf1 in apoptosomes
# C9 = caspase 9, 2 binding sites for Apaf1 and XIAP, can be pro or cleaved by Apaf1
# XIAP = X linked inhibitor of apoptosis protein
# C3 = caspase 3, 2 binding sites for C9 and XIAP, can be pro or cleaved by C9

"""Definitions listed in the program"""
# momp_monomers : declares the monomers Bax, Bak, Bid, Bcl2, BclxL, Mc11, Bad, Noxa, CytoC, and Smac
# downstream_monomers : declares the monomers Apaf1, C9, C3, and XIAP
# declare_initial_conditions : declares the parameters for the initial conditions and sets the initial conditions to that parameter
# declare_reaction_rates : declare the forward and reverse binding rates, catalytic rates, and synthesis and degradation rates
# translocate_tBid_Bax_BclxL : tBid, Bax, and BclXL translocate to the mitochondrial membrane
# tBid_activates_Bax_and_Bak : bid activates bax translocating it to the membrane; tBid activates Bax and Bak
# tBid_binds_all_anti_apoptotics : tBid binds and inhibits Bcl2, Mcl1, and Bcl-XL.
# sensitizers_bind_anti_apoptotics : Binding of Bad and Noxa to Bcl2, Mcl1, and Bcl-XL.
# effectors_bind_anti_apoptotics : Binding of Bax and Bak to Bcl2, BclxL, and Mcl1.
# lopez_pore_formation(do_pore_transport=True) : Pore formation and transport process used by all modules
# C3_and_C9_cleavage_by_apoptosome_with_XIAP_inhibition : Apaf1 binds C9, and the cleavage of C3 and C9 are catalyzed and XIAP blocks activity



# import the pysb module and all its methods and functions
from pysb import *
import matplotlib.pyplot as plt
from pysb.bng import *
from pysb.bng import generate_equations
from pysb.bng import generate_network
from pysb.simulator import ScipyOdeSimulator
import pylab as pl
import pandas as pd
from pysb.util import alias_model_components

# _______________________________________________________________________________________________________________________________________________________________________________________________#

# Instantiate a model
Model()

"""Declare Caspases and related Monomers"""
def downstream_monomers():
    Monomer('Apaf1', ['bC9'])
    # create monomer Apaf1 with single binding site for CC
    Monomer('C9', ['bApaf1', 'bXIAP', 'm'], {'m': ['p', 'c']})
    # create monomer C9 with a binding and cleavage site, the cleavage site can be pro (before cleavage, or cleaved (after cleavage)
    Monomer('XIAP', ['bC9C3'])
    # create monomer XIAP with a single binding site for caspase 3, caspase 9, procaspase 3, and procaspase 9, wihch bind competitively, inhiitis caspase 3 and 9
    Monomer('C3', ['bXIAP','m'], {'m':['p','c']})
    # create monomer C3 with binding site a for C9, b binding site for XIAP, and c cleavage site with states pro or cleaved

    alias_model_components()

def momp_monomers():
    """Activator Monomers"""
    Monomer('Bid', ['bf', 'state'], {'state':['U', 'T', 'M']})
    # create monomer Bid of the Bcl2 family with states U for untruncated, T for truncated, and M for truncated and in the mitochondria, bf is a binding site for protein protein interactions

    """Effector Monomers"""
    Monomer('Bax', ['bf', 's1', 's2', 'state'], {'state':['C', 'M', 'A']})
    # Bax, states: Cytoplasmic, Mitochondrial, Active
    # sites 's1' and 's2' are used for pore formation, bf is a binding site for protein protein interactions
    Monomer('Bak', ['bf', 's1', 's2', 'state'], {'state':['M', 'A']})
    # Bak, states: inactive and Mitochondrial, Active (and mitochondrial)
    # sites 's1' and 's2' are used for pore formation, bf is a binding site for protein protein interactions

    """Anti-Apoptotic Monomers"""
    Monomer('Bcl2', ['bf','state'], {'state':['C','M']})
    Monomer('BclxL', ['bf', 'state'], {'state':['C', 'M']})
    Monomer('Mcl1', ['bf', 'state'], {'state':['C', 'M']})

    """Sensitizer Monomers"""
    # C and M stand for cytosolic and mitochondrial states, bf is a binding site for protein protein interactions
    Monomer('Bad', ['bf', 'state'], {'state':['C', 'M']})
    Monomer('Noxa', ['bf', 'state'], {'state': ['C', 'M']})

    """Cytochrome C and Smac"""
    # M is for mitochondrial, C for cytosolic, and A is for active, bf is a binding site for protein protein interactions
    Monomer('CytoC', ['bf', 'state'], {'state': ['M', 'C', 'A']})
    Monomer('Smac', ['bf', 'state'], {'state': ['M', 'C', 'A']})

    alias_model_components()

# ______________________________________________________________________________________________________________________________________________________________________________________________#


"""Set Initial Conditions"""
def declare_initial_conditions():
    Parameter('Apaf1_0', 20)          # Apaf1
    Parameter('C9p_0', 20)            # Procaspase 9
    Parameter('XIAP_0', 40)           # XIAP
    Parameter('C3p_0', 200)           # Procaspase 3
    Parameter('Bid_0', 4.0e4)         # Bid
    Parameter('BclxL_0', 2.0e4)       # cytosolic BclxL
    Parameter('Mcl1_0', 2.0e4)        # Mitochondrial Mcl1
    Parameter('Bcl2_0', 2.0e4)        # Mitochondrial Bcl2
    Parameter('Bad_0', 1.0e3)         # Bad
    Parameter('Noxa_0', 1.0e3)        # Noxa
    Parameter('CytoC_0', 5.0e5)       # cytochrome c
    Parameter('Smac_0', 1.0e5)        # Smac
    Parameter('Bax_0', 0.8e5)         # Bax
    Parameter('Bak_0', 0.2e5)         # Bak

    alias_model_components()

    Initial(Apaf1(bC9=None), Apaf1_0)
    Initial(C9(bXIAP=None, bApaf1=None, m='p'), C9p_0)
    Initial(C3(bXIAP=None, m='p'), C3p_0)
    Initial(XIAP(bC9C3=None), XIAP_0)
    Apaf1_num = [5, 6.5, 10, 20]
    #these are the different concentrations of Apaf1
    Initial(Bid(bf=None, state='U'), Bid_0)
    Initial(Bad(bf=None, state='M'), Bad_0)
    Initial(Bax(bf=None, s1=None, s2=None, state='C'), Bax_0)
    Initial(Bak(bf=None, s1=None, s2=None, state='M'), Bak_0)
    Initial(Bcl2(bf=None, state='M'), Bcl2_0)
    Initial(BclxL(bf=None, state='M'), BclxL_0)
    Initial(Mcl1(bf=None, state='M'), Mcl1_0)
    Initial(Noxa(bf=None, state='M'), Noxa_0)
    Initial(CytoC(bf=None, state='M'), CytoC_0)
    Initial(Smac(bf=None, state='M'), Smac_0)


def declare_reaction_rates
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
    Parameter('degrade_C9p_XIAP_kd', 1e-03)                   # s^-1
    Parameter('degrade_Apaf1_C9p_XIAP_kd', 1e-03)             # s^-1
    Parameter('degrade_Apaf1_C9p_kd', 1e-03)                  # s^-1
    Parameter('degrade_C3c_XIAP_kd', 1e-03)                   # s^-1
    Parameter('degrade_C9c_XIAP_kd', 1e-03)                   # s^-1
    Parameter('degrade_Apaf1_C9c_kd', 1e-03)                  # s^-1
    Parameter('degrade_Apaf1_C9c_XIAP_kd', 1e-03)             # s^-1

    alias_model_components()

# ______________________________________________________________________________________________________________________________________________________#

"""MOMP modules from EARM"""

def translocate_tBid_Bax_BclxL():
    """tBid, Bax and BclXL translocate to the mitochondrial membrane."""
    equilibrate(Bid(bf=None, state='T'), Bid(bf=None, state='M'), [1e-1, 1e-3])

def tBid_activates_Bax_and_Bak():
    """tBid activates Bax and Bak."""
    # bid activates bax translocating it to the membrane
    catalyze(Bid(state='M'), Bax(state='C'), Bax(state='M'), activation_rates)
    catalyze(Bid(state='M'), Bax(state='M'), Bax(state='A'), activation_rates)
    catalyze(Bid(state='M'), Bak(state='M'), Bak(state='A'), activation_rates)

def tBid_binds_all_anti_apoptotics():
    # tBid binds and inhibits Bcl2, Mcl1, and Bcl-XL.
    bind_table([[Bcl2(state='M'), BclxL(state='M'), Mcl1(state='M')],
                [Bid(state='M'), 66e-9 * N_A * V, 12e-9 * N_A * V, 10e-9 * N_A * V]],
               kf=1e6 / (N_A * V))

def sensitizers_bind_anti_apoptotics():
    # Binding of Bad and Noxa to Bcl2, Mcl1, and Bcl-XL.
    bind_table([[                        Bcl2,  BclxL(state='M'),  Mcl1(state='M')],
                [Bad(state='M'),  11e-9*N_A*V,       10e-9*N_A*V,             None],
                [Noxa(state='M'),        None,              None,      19e-9*N_A*V]],
               kf=1e6/(N_A*V))

def effectors_bind_anti_apoptotics():
    # Binding of Bax and Bak to Bcl2, BclxL, and Mcl1.
    bind_table([[                            Bcl2,  BclxL(state='M'),Mcl1(state='M')],
                [Bax(active_monomer), 10e-9*N_A*V,       10e-9*N_A*V,         None],
                [Bak(active_monomer),        None,       50e-9*N_A*V,  10e-9*N_A*V]],
               kf=1e6/(N_A*V))

def lopez_pore_formation(do_pore_transport=True):
    # Pore formation and transport process used by all modules.
    alias_model_components()

    # Rates
    pore_max_size = 4
    pore_rates = [[2.040816e-04,  # 1.0e-6/v**2
                   1e-3]] * (pore_max_size - 1)
    pore_transport_rates = [[2.857143e-5, 1e-3, 10]] # 2e-6 / v?

    # Pore formation by effectors
    assemble_pore_sequential(Bax(bf=None, state='A'), pore_max_size, pore_rates)
    assemble_pore_sequential(Bak(bf=None, state='A'), pore_max_size, pore_rates)

    # CytoC, Smac release
    if do_pore_transport:
        pore_transport(Bax(bf=None, state='A'), pore_max_size, CytoC(state='M'),
                       CytoC(state='C'), pore_transport_rates)
        pore_transport(Bax(bf=None, state='A'), pore_max_size, Smac(state='M'),
                       Smac(state='C'), pore_transport_rates)
        pore_transport(Bak(bf=None, state='A'), pore_max_size, CytoC(state='M'),
                       CytoC(state='C'), pore_transport_rates)
        pore_transport(Bak(bf=None, state='A'), pore_max_size, Smac(state='M'),
                       Smac(state='C'), pore_transport_rates)

# ________________________________________________________________________________________________________________________________________________________________________________________#
def C3_and_C9_cleavage_by_apoptosome_with_XIAP_inhibition():
    """Rules for activated Apaf1 equilibrium binding"""
    #Reaction 1, 8, 12, 10, 13, 14 from paper Bistability in Caspase Activation, 2006, PLoS Computational Biology
    Rule('Apaf1_bind_C9p', Apaf1(bC9=None) + C9(bApaf1=None, bXIAP=None, m='p') | Apaf1(bC9=1) % C9(bApaf1=1, bXIAP=None, m='p'), Apaf1_bind_C9p_kf, Apaf1_bind_C9p_kr)
    # Activated Apaf1 reversibly binds procaspase 9
    Rule('Apaf1_bind_C9c', Apaf1(bC9=None) + C9(bApaf1=None, bXIAP=None, m='c') | Apaf1(bC9=1) % C9(bApaf1=1, bXIAP=None, m='c'), Apaf1_bind_C9c_kf, Apaf1_bind_C9c_kr)
    # Activated Apaf1 reversibly binds caspase 9
    Rule('Apaf1_C9c_bind_XIAP', Apaf1(bC9=1) % C9(bApaf1=1, bXIAP=None, m='c') + XIAP(bC9C3=None) | Apaf1(bC9=1) % C9(bApaf1=1, bXIAP=2, m='c') % XIAP(bC9C3=2), Apaf1_C9c_bind_XIAP_kf, Apaf1_C9c_bind_XIAP_kr)
    # the Apaf1 caspase 9 complex binds XIAP
    Rule('Apaf1_C9p_bind_XIAP', Apaf1(bC9=1) % C9(bApaf1=1, bXIAP=None, m='p') + XIAP(bC9C3=None) | Apaf1(bC9=1) % C9(bApaf1=1, bXIAP=2, m='p') % XIAP(bC9C3=2), Apaf1_C9p_bind_XIAP_kf, Apaf1_C9p_bind_XIAP_kr)
    # the Apaf1 procaspase 9 complex binds XIAP
    Rule('C9p_XIAP_bind_Apaf1', Apaf1(bC9=None) + C9(bApaf1=None, bXIAP=1, m='p') % XIAP(bC9C3=1) | Apaf1(bC9=2) % C9(bApaf1=2, bXIAP=1, m='p') % XIAP(bC9C3=1), C9p_XIAP_bind_Apaf1_kf, C9p_XIAP_bind_Apaf1_kr)
    # the procaspase 9 XIAP complex binds Afap1
    Rule('C9c_XIAP_bind_Apaf1', Apaf1(bC9=None) + C9(bApaf1=None, bXIAP=1, m='c') % XIAP(bC9C3=1) | Apaf1(bC9=2) % C9(bApaf1=2, bXIAP=1, m='c') % XIAP(bC9C3=1), C9c_XIAP_bind_Apaf1_kf, C9c_XIAP_bind_Apaf1_kr)
    # the caspase 9 XIAP complex binds Apaf1

    """Rules for C9 & C3 equilibrium binding"""
    # Reactions 9, 11, and 15 from paper Bistability in Caspase Activation, 2006, PLoS Computational Biology
    Rule('C9p_bind_XIAP', C9(bXIAP=None, bApaf1=None, m='p') + XIAP(bC9C3=None) | C9(bXIAP=1, bApaf1=None, m='p') % XIAP(bC9C3=1), C9p_bind_XIAP_kf, C9p_bind_XIAP_kr)
    # Procaspase 9 binds XIAP
    Rule('C9c_bind_XIAP', C9(bXIAP=None, bApaf1=None, m='c') + XIAP(bC9C3=None) | C9(bXIAP=1, bApaf1=None, m='c') % XIAP(bC9C3=1), C9c_bind_XIAP_kf, C9c_bind_XIAP_kr)
    # Caspase 9 binds XIAP
    Rule('C3c_bind_XIAP', C3(bXIAP=None, m='c') + XIAP(bC9C3=None) | C3(bXIAP=1, m='c') % XIAP(bC9C3=1), C3c_bind_XIAP_kf, C3c_bind_XIAP_kr)
    # Caspase 3 binds XIAP

    """Rules for Catalysis of Procaspase 3 Cleavage"""
    # Reactions 2, 3, 6, and 7 from paper Bistability in Caspase Activation, 2006, PLoS Computational Biology
    Rule('C9p_cleave_C3p', C3(bXIAP=None, m='p') + C9(bApaf1=None, bXIAP=None, m='p') >> C3(bXIAP=None, m='c') + C9(bApaf1=None, bXIAP=None, m='p'), C9p_cleave_C3p_kc)
    # procaspase 9 catalyses cleavage of procaspase 3 to caspase 3
    Rule('C9c_cleave_C3p', C3(bXIAP=None, m='p') + C9(bApaf1=None, bXIAP=None, m='c') >> C3(bXIAP=None, m='c') + C9(bApaf1=None, bXIAP=None, m='c'), C9c_cleave_C3p_kc)
    # Caspase 9 catalyses cleavage of procaspase 3 to caspase 3C9p_XIAP_bind_Apaf1_kr
    Rule('C9p_Apaf1_cleave_C3p', C3(bXIAP=None, m='p') + Apaf1(bC9=1) % C9(bApaf1=1, bXIAP=None, m='p') >> C3(bXIAP=None, m='c') + Apaf1(bC9=1) % C9(bApaf1=1, bXIAP=None, m='p'), C9p_Apaf1_cleave_C3p_kc)
    # Complex of activated Apaf1 and procaspase 9 catalyze cleavage of procaspase 3 to caspase 3
    Rule('C9c_Apaf1_cleave_C3p', C3(bXIAP=None, m='p') + Apaf1(bC9=1) % C9(bApaf1=1, bXIAP=None, m='c') >> C3(bXIAP=None, m='c') + Apaf1(bC9=1) % C9(bApaf1=1, bXIAP=None, m='c'), C9c_Apaf1_cleave_C3p_kc)
    # Complex of activated Apaf1 and caspase 9 catalyze cleavage of procaspase 3 to caspase 3

    """Rules for Catalysis of Procaspase 9 Cleavage"""
    # Reactions 4 and 5 from paper Bistability in Caspase Activation, 2006, PLoS Computational Biology
    Rule('C3c_cleave_C9p', C9(bApaf1=None, bXIAP=None, m='p') + C3(bXIAP=None, m='c') >> C9(bApaf1=None, bXIAP=None, m='c') + C3(bXIAP=None, m='c'), C3c_cleave_C9p_kc)
    # Caspase 3 catalyzes the cleavage of procaspase 9 to caspase 9
    Rule('C3c_cleave_C9p_Apaf1', C3(bXIAP=None, m='c') + C9(bApaf1=1, bXIAP=None, m='p') % Apaf1(bC9=1) >> C3(bXIAP=None, m='c') + C9(bApaf1=1, bXIAP=None, m='c') % Apaf1(bC9=1), C3c_cleave_C9p_Apaf1_kc)
    # Caspase 3 catalyzes the cleavage of procaspase 9 bound to Apaf1 to caspase 9 bound to Apaf1


def synthesize_degrade_downstream_monomers():
    #Rules for the Synthesis and Degradation of all single protein species
    # Reactions 16, 17, 18, 22, 23, and 26 from paper Bistability in Caspase Activation, 2006, PLoS Computational Biology
    Rule('synthesize_Apaf1', None >> Apaf1(bC9=None), synthesize_Apaf1_ks)
    # Make activated Apaf1
    Rule('synthesize_C9p', None >> C9(bApaf1=None, bXIAP=None, m='p'), synthesize_C9p_ks)
    # Make procaspase 9
    Rule('synthesize_XIAP', None >> XIAP(bC9C3=None), synthesize_XIAP_ks)
    # Make XIAP
    Rule('synthesize_C3p', None >> C3(bXIAP=None, m='p'), synthesize_C3p_ks)
    # Make procaspase 3
    Rule('degrade_Apaf1', Apaf1(bC9=None) >> None, degrade_Apaf1_kd)
    # Degrade Apaf1
    Rule('degrade_C9p', C9(bApaf1=None, bXIAP=None, m='p') >> None, degrade_C9p_kd)
    # Degrade procaspase 9
    Rule('degrade_C9c', C9(bApaf1=None, bXIAP=None, m='c') >> None, degrade_C9c_kd)
    # Degrade caspase 9
    Rule('degrade_C3p', C3(bXIAP=None, m='p') >> None, degrade_C3p_kd)
    # Degrade procaspase 3
    Rule('degrade_C3c', C3(bXIAP=None, m='c') >> None, degrade_C3c_kd)
    # Degrade caspase 3
    Rule('degrade_XIAP', XIAP(bC9C3=None) >> None, degrade_XIAP_kd)
    # Degrade XIAP

def degrade_downstream_complexes():
    """Rules for Degradation of Complexes"""
    # Reactions 19, 20, 21, 24, 25, 27, and 28 from paper Bistability in Caspase Activation, 2006, PLoS Computational Biology
    Rule('degrade_C9p_XIAP', C9(bXIAP=1, bApaf1=None, m='p') % XIAP(bC9C3=1) >> None, degrade_C9p_XIAP_kd)
    #Degrade the procaspase 9 XIAP complex
    Rule('degrade_Apaf1_C9p_XIAP', C9(bXIAP=1, bApaf1=2, m='p') % XIAP(bC9C3=1) % Apaf1(bC9=2) >> None, degrade_Apaf1_C9p_XIAP_kd)
    # Degrade the activated Apaf1, procaspase 9, XIAP complex
    Rule('degrade_Apaf1_C9p', Apaf1(bC9=1) % C9(bApaf1=1, bXIAP=None, m='p') >> None, degrade_Apaf1_C9p_kd)
    # Degrade the Apaf1 procaspase 9 complex
    Rule('degrade_C3c_XIAP', C3(bXIAP=1, m='c') % XIAP(bC9C3=1) >> None, degrade_C3c_XIAP_kd)
    # Degrade the caspase 3 XIAP complex
    Rule('degrade_C9c_XIAP', C9(bApaf1=None, bXIAP=1, m='c') % XIAP(bC9C3=1) >> None, degrade_C9c_XIAP_kd)
    # Degrade the caspase 9 XIAP complex
    Rule('degrade_Apaf1_C9c', Apaf1(bC9=1) % C9(bXIAP=None, bApaf1=1, m='c') >> None, degrade_Apaf1_C9c_kd)
    # Degrade the activated Apaf1 caspase 9 complex
    Rule('degrade_Apaf1_C9c_XIAP', Apaf1(bC9=1) % C9(bXIAP=2, bApaf1=1, m='c') % XIAP(bC9C3=2) >> None, degrade_Apaf1_C9c_XIAP_kd)
    # Degrade the activated Apaf1 caspase 9 XIAP complex


# ______________________________________________________________________________________________________________________________________________________________________________________________#


"""Set the observables and Setup a Simulation"""
Observable('ObsC3c', C3(bXIAP=None, m='c'))
Observable('Apaf1_obs', Apaf1(bC9=None))




generate_network(model)
generate_equations(model)
print(model.species)
for i,sp in enumerate(list(model.species)):
    print i,":",sp

# for i,ode in enumerate(list(model.odes)):
#     print i, ":", ode

#Dictionary to substitute in species names
species_dict = {0 : 'Apaf1',
                1 : 'C9p',
                2 : 'XIAP',
                3 : 'C3p',
                4 : 'Apaf1-C9p',
                5 : 'C9p-XIAP',
                6 : 'C3c',
                7 : 'Apaf1-C9p-XIAP',
                8 : 'C3c-XIAP',
                9 : 'C9c',
                10 : 'Apaf1-C9c',
                11 : 'Apaf1-C9c-XIAP',
                12 : 'C9c-XIAP'
 }

for ode in list(model.odes):
    for i in range(len(model.species)):
        ode = re.sub(r'\b__s%d\b'%i, species_dict[i], str(ode))
    print ode



"""Set the observables and Setup a Simulation"""
# Observable('ObsC3c', C3(bXIAP=None, m='c'))
# Observable('Apaf1_obs', Apaf1(bC9=None))


faddnum = []
t = pl.linspace(0,4800,4801)
simres = ScipyOdeSimulator(model, tspan=t)
simres_result = simres.run(initials = {Apaf1(bC9=None): Apaf1_num})
simres_df = simres_result.dataframe
# yout = simres_result.all
print simres_df['ObsC3c']



# """Plot the Observables"""
pl.ion()
pl.figure()
for n in range(0,4):
    plt.plot(t/60, simres_df.loc[n]['ObsC3c'].iloc[:], lw=1)
#pl.plot(t/60, simres_df['ObsC3c'], label='Activated Caspase 3')
# #pl.plot(t/60, simres.observables['ObsC3c'], label='Activated Caspase 3')
# pl.plot(t/60, yout['__s6'], label='Activated Caspase 3 for real')
# pl.plot(t/60, simres.observables['Apaf1_obs'], label='Apaf1')
# pl.plot(t/60, yout['__s0'], label='Apaf1sp')
# pl.legend(loc='best')
# pl.xlabel('Time (s)')
# pl.ylabel('Concentration in nM')
# pl.show()

