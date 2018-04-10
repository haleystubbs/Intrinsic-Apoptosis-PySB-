

from pysb import *
from pysb.util import alias_model_components
from sympy import sin
from sympy import Piecewise
from pysb.macros import catalyze
import matplotlib.pyplot as plt
from pysb.bng import *
from pysb.bng import generate_equations
from pysb.bng import generate_network
from pysb.simulator import ScipyOdeSimulator
import pylab as pl
import pandas as pd

from earm import shared
from earm import lopez_modules
from earm import albeck_modules
"""_________________________________________________________________________________________________________________________________________"""

Model()

Monomer('mdm2', ['bf'])
Monomer('p53', ['bf', 'state'], {'state': ['I', 'A']})
Monomer('Puma', ['bf', 'state'], {'state': ['I', 'A']})
Monomer('BIM', ['bf', 'state'], {'state': ['I', 'A']})
# Bax, states: Cytoplasmic, Mitochondrial, Active
# sites 's1' and 's2' are used for pore formation
Monomer('t')
Monomer('time_delay_mon')
Monomer('signal1')
Monomer('signal2')
Monomer('signal3')
Monomer('signal4')
Monomer('signal5')
Monomer('signal6')
Monomer('signal7')
Monomer('signal8')
Monomer('signal9')

albeck_modules.ligand_to_c8_monomers()
lopez_modules.momp_monomers()
albeck_modules.apaf1_to_parp_monomers()

Parameter('time_delay_k', 1000)
Parameter('activate_p53_k', 5000)             # per hour
Parameter('transcribe_Puma_k', 900)
Parameter('transcribe_BIM_k', 900)
Parameter('transcribe_Bid_k', 900)
Parameter('transcribe_Bax_k', 900)
Parameter('transcribe_Bak_k', 900)
Parameter('transcribe_mdm2_k', 900)
# Parameter('mdm2_bind_p53I_kf', )
# Parameter('mdm2_bind_p53I_kr', )
# Parameter('mdm2_bind_p53A_kf', )
# Parameter('mdm2_bind_p53A_kr', )
# Parameter('signal1_exp', 2)
# Parameter('synthesize_mdm2_independently_k', synthesize_mdm2_independently_k_exp)
Parameter('degrade_mdm2_independently_k', 1000)
# Parameter('synthesize_p53I_k', synthesize_p53I_k_exp)
Parameter('degrade_p53I_independently_k', 2000)
# Parameter('mdm2_degrade_p53I_k', mdm2_degrade_p53I_k_exp)
# Parameter('mdm2_degrade_p53A_k', mdm2_degrade_p53A_k_exp)
# Parameter('signal_degrade_mdm2', signal_degrade_mdm2_k_exp)
activation_rates = [        1e-7, 1e-3, 1]
Parameter('n', 4)                            # Hill coefficient
Parameter('T', 900)                           # the signal concentration at half maximal p53 activation = 0.75(maximum concentration of signal
Parameter('frequency1', 1)                   # frequency of signal sine wave in units of per hour
Parameter('frequency2', 3)                   # frequency of signal sine wave in units of per hour
Parameter('frequency3', 6)                   # frequency of signal sine wave in units of per hour
Parameter('max_signal', 1200)                 # maximum amount of signal throughout the simulation, also the amplitude for the oscillatory signal nd height of pulsed signals
Parameter('time_rate', 0.2)
Expression('simulated_conc_unit', (3/4)*max_signal)
Parameter('Mmax', 3000)                    # the maximum amount of mdm2
Parameter('Pmax', 3000)                    # the maximum amount of p53
Parameter('time_delay_initial', 0.7)
Rule('time', None >> t(), time_rate)         #time points every 0.2 of an hour
Parameter('signal3_rate', 0)
Parameter('mdm2_initial', 500)
Parameter('Puma_initial', 1000)
Parameter('p53I_initial', 800)
Parameter('p53A_initial', 0)
# Initial(p53(bf=None, state='I'), p53I_initial)
Initial(mdm2(bf=None), mdm2_initial)
Initial(signal3(), max_signal)
Initial(signal7(), max_signal)
Initial(signal8(), max_signal)
Initial(signal9(), max_signal)
Initial(time_delay_mon(), time_delay_initial)
Initial(Puma(bf=None, state='I'), Puma_initial)
Initial(p53(bf=None, state='A'), p53A_initial)
Observable('Obs_p53', p53(bf=None, state='A'))
Observable('Obs_time', t())
# Observable('Obs_signal1', signal1)
# Observable('Obs_signal2', signal2())
Observable('Obs_signal3', signal3())
Observable('Obs_signal4', signal4())
Observable('Obs_signal5', signal5())
Observable('Obs_signal6', signal6())
Observable('Obs_signal7', signal7())
Observable('Obs_signal8', signal8())
Observable('Obs_signal9', signal9())
Observable('time_delay', time_delay_mon())
Observable('Obs_p53I', p53(state='I'))
Observable('Obs_p53A', p53(state='A'))
Observable('Obs_p53I_unbound', p53(bf=None, state='I'))
Observable('Obs_p53A_unbound', p53(bf=None, state='A'))
Observable('Obs_mdm2', mdm2())
Observable('Obs_mdm2_unbound', mdm2(bf=None))
Observable('Obs_p53_mdm2', p53(bf=1) % mdm2(bf=1))
Observable('Obs_Puma', Puma(state='A'))
Observable('Obs_time_delay', time_delay_mon())
Observable('Obs_CytC_release', CytoC(state='C') and CytoC(state='A'))


"""__________________________________________________________________________________________________________________________________________________________"""

Expression('signal_factor', Obs_signal4 ** n / (Obs_signal4 ** n + T ** n))

"""Expressions for rates as a function of maximum protein concentration"""


Expression('synthesize_mdm2_independently_k_exp', 0.27*Mmax)
Expression('synthesize_p53I_k_exp', 0.87*Pmax)
Expression('mdm2_degrade_p53I_k_exp', (3.76/Mmax)*Obs_mdm2_unbound)
Expression('mdm2_degrade_p53A_k_exp', 1.05/Mmax*Obs_mdm2_unbound)
Expression('signal_degrade_mdm2_k_exp', 0.67/max_signal*Obs_signal4)
Rule('time_delay_rule', None >> time_delay_mon(), time_delay_k)
Expression('activate_p53_k_exp', activate_p53_k*signal_factor)
Expression('transcribe_BIM_k_exp', transcribe_BIM_k*Obs_p53A*time_delay)
Expression('transcribe_Puma_k_exp', transcribe_Puma_k*Obs_p53A*time_delay)
Expression('transcribe_Bid_k_exp', transcribe_Bid_k*Obs_p53A*time_delay)
Expression('transcribe_Bax_k_exp', transcribe_Bax_k*Obs_p53A*time_delay)
Expression('transcribe_Bak_k_exp', transcribe_Bak_k*Obs_p53A*time_delay)
Expression('transcribe_mdm2_k_exp', transcribe_mdm2_k*Obs_p53A*time_delay)



"""______________________________________________________________________________________________________________________________________________________________-_"""

"""Signal Component for ODEs and Signal Input Functions"""

# Expression('signal1_exp', Piecewise((max_signal, (1 <= Obs_time) & (Obs_time <= 2)),
#                                 (max_signal, (6 <= Obs_time) & (Obs_time <= 7)),
#                                 (max_signal, (11 <= Obs_time)& (Obs_time <= 12)),
#                                 (max_signal, (16 <= Obs_time)& (Obs_time <= 17)),
#                                 (max_signal, (21 <= Obs_time)& (Obs_time <= 22)),
#                                 (0, True)))
# Expression('signal2_exp', Piecewise((max_signal, (1 <= Obs_time) & (Obs_time <= 4)),
#                                 (max_signal, (8 <= Obs_time) & (Obs_time <= 11)),
#                                 (max_signal, (15 <= Obs_time)& (Obs_time <= 19)),
#                                 (0, True)))
Expression('signal3_exp', max_signal)
Expression('signal4_exp', (max_signal*frequency1/10*3.14159)*sin(frequency1*3.14159*Obs_time))
Expression('signal5_exp', (max_signal*frequency2/10*3.14159)*sin(frequency2*3.14159*Obs_time))
Expression('signal6_exp', (max_signal*frequency3/10*3.14159)*sin(frequency3*3.14159*Obs_time))
Expression('signal7_exp', (max_signal*1.3863) * 2.71828**(-1.3863*Obs_time))
Expression('signal8_exp', (max_signal*0.693147) * 2.71828**(-0.693147*Obs_time))
Expression('signal9_exp', (max_signal*0.11552) * 2.71828**(-0.11552*Obs_time))


# Rule('generate_signal1', None >> signal1(), signal1_exp)
# Rule('generate_signal2', None >> signal2(), signal2_exp)
Rule('generate_signal3', None >> signal3(), signal3_rate)
Rule('generate_signal4', None >> signal4(), signal4_exp)
Rule('generate_signal5', None >> signal5(), signal5_exp)
Rule('generate_signal6', None >> signal6(), signal6_exp)
Rule('generate_signal7', signal7() >> None, signal7_exp)
Rule('generate_signal8', signal8() >> None, signal8_exp)
Rule('generate_signal9', signal9() >> None, signal9_exp)


"""_____________________________________________________________________________________________________________________________________________________________"""



"""Rules for p53 mdm2 loop"""



Rule('synthesize_p53I', None >> p53(bf=None, state='I'), synthesize_p53I_k_exp)
# synthesize mdm2 independently of p53 transcriptional activity
Rule('synthesize_mdm2_independently', None >> mdm2(bf=None), synthesize_mdm2_independently_k_exp)
# degrade mdm2 independently of signaling
Rule('degrade_mdm2_independently', mdm2(bf=None) >> None, degrade_mdm2_independently_k)
# degrade p53I independently of mdm2
Rule('degrade_p53I_independently', p53(bf=None, state='I')>> None, degrade_p53I_independently_k)
# p53 is activated by the signal
Rule('activate_p53', p53(bf=None, state='I') >> p53(bf=None, state='A'), activate_p53_k_exp)
# mdm2 binds p53I and degrades it
Rule('mdm2_degrade_p53I', p53(bf=None, state='I') >> None, mdm2_degrade_p53I_k_exp)
# mdm2 binds and degrades p53A
Rule('mdm2_degrade_p53A', p53(bf=None, state='A') >> None, mdm2_degrade_p53A_k_exp)
# mdm2 binds inactive p53 and prevents it from being activated
Rule('signal_degrade_mdm2', mdm2(bf=None) >> None, signal_degrade_mdm2_k_exp)



"""Transcribe proteins"""

Rule('transcribe_mdm2', None >> mdm2(bf=None), transcribe_mdm2_k_exp)
# p53 activates PUMA by transcribing it
Rule('transcribe_Puma', Puma(bf=None, state='I') >> Puma(bf=None, state='A'), transcribe_Puma_k_exp)
# p53 activates BIM by transcribing it
Rule('transcribe_BIM', BIM(bf=None, state='I') >> BIM(bf=None, state='A'), transcribe_BIM_k_exp)
# p53 activates Bid by transcribing it
Rule('transcribe_Bid', None >> Bid(bf=None, state='U'), transcribe_Bid_k_exp)
# p53 is activated and starts transcribing Bak, and Bax
Rule('transcribe_Bax', None >> Bax(bf=None, s1=None, s2=None, state='C'), transcribe_Bax_k_exp)
Rule('transcribe_Bak', None >> Bak(bf=None, s1=None, s2=None, state='M'), transcribe_Bak_k_exp)



"""BIM and PUMA activate Bax and Bak to initiate MOMP"""
# BIM and Puma activate Bax so that it inserts into the mitochondrial membrane
catalyze(BIM(state='A'),  'bf', Bax(state='C'), 'bf', Bax(state='M'), activation_rates)
catalyze(Puma(state='A'), 'bf', Bax(state='C'), 'bf', Bax(state='M'), activation_rates)
# BIM and Puma activate Bax in the mitochondrial membrane so it will be ready to form pores
catalyze(BIM(state='A'),  'bf', Bax(state='M'), 'bf', Bax(state='A'), activation_rates)
catalyze(Puma(state='A'), 'bf', Bax(state='M'), 'bf', Bax(state='A'), activation_rates)
# BIM and Puma activate Bak so that it will be ready to form pores
catalyze(BIM(state='A'),  'bf', Bak(state='M'), 'bf', Bak(state='A'), activation_rates)
catalyze(Puma(state='A'), 'bf', Bak(state='M'), 'bf', Bak(state='A'), activation_rates)


"""Call MOMP modules from Earm to connect to downstream effectors"""



# Generate the upstream and downstream sections
albeck_modules.rec_to_bid()
albeck_modules.pore_to_parp()

# The specific MOMP model to use
lopez_modules.embedded()

shared.observables()


"""_________________________________________________________________________________________________________________________________________________________________"""

# for i, x in enumerate(signal_expression):       # this will assign an index number to each of the signals in the signal_expression array
#     print  '%s is %d' % (x, i)
#
# for x in signal_expression:                    # loop through the different signals in the signal_expression array and save all the data to a separate array
#     Rule('generate_signal', None >> signal(), x)
#     Observable('signal%d' % (i), signal())
#     Expression('signal_factor', Obs_signal ** n / (Obs_signal ** n + T ** n))
#     signal_%d % (i) = [x]
#     # p53_mdm2_loop()
#     alias_model_components()






"""Set up a Simulation"""

faddnum = []
simulation_time = pl.linspace(0, 24, 240)
simres = ScipyOdeSimulator(model, tspan=simulation_time).run()
yout = simres.all


"""Plot the Observables"""
pl.ion()
pl.figure()
pl.plot(simulation_time, yout['Obs_signal4'], label='Signal 4')
# pl.plot(simulation_time, yout['Obs_signal4'], label='Signal 4')
# pl.plot(simulation_time, yout['Obs_signal5'], label='Signal 5')
# pl.plot(simulation_time, yout['Obs_signal6'], label='Signal 6')
# pl.plot(simulation_time, yout['Obs_signal7'], label='Signal 7')
# pl.plot(simulation_time, yout['Obs_signal8'], label='Signal 8')
# pl.plot(simulation_time, yout['Obs_signal9'], label='Signal 9')
# pl.plot(simulation_time, yout['Obs_p53I'], label='p53 inactive')
pl.plot(simulation_time, yout['Obs_p53A'], label='p53 active')
# pl.plot(simulation_time, yout['Obs_mdm2_unbound'], label='free mdm2')
# pl.plot(simulation_time, yout['Obs_p53_mdm2'], label='p53 bound mdm2')
# pl.plot(simulation_time, yout['Obs_Puma'], label='Puma')
pl.plot(simulation_time, yout['Obs_CytC_release'], label='CytoC Release')
# pl.plot(simulation_time, yout['Obs_time_delay'], label='time delay')
pl.legend(loc='best')
pl.xlabel('Time (h)')
pl.ylabel('Concentration in uM')
pl.show(block=True)


