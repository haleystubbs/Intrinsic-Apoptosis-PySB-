

from pysb import *
from pysb.util import alias_model_components
from sympy import sin
from sympy import Piecewise
from pysb.macros import catalyze
from sympy.abc import x
"""_________________________________________________________________________________________________________________________________________"""

Model()

Monomer('mdm2', ['bf', 'state'], {'state':['I', 'A']})
Monomer('p53', ['bf', 'state'], {'state': ['I', 'A']})
Monomer('PUMA', ['bf', 'state'], {'state': ['I', 'A']})
Monomer('BIM', ['bf', 'state'], {'state': ['I', 'A']})
Monomer('Bid', ['bf', 'state'], {'state':['U', 'T', 'M']})      # untruncated, truncated, and truncated and mitochondrial
# Bax, states: Cytoplasmic, Mitochondrial, Active
# sites 's1' and 's2' are used for pore formation
Monomer('Bax', ['bf', 's1', 's2', 'state'], {'state': ['C', 'M', 'A']})
# Bak, states: inactive and Mitochondrial, Active (and mitochondrial)
# sites 's1' and 's2' are used for pore formation
Monomer('Bak', ['bf', 's1', 's2', 'state'], {'state': ['M', 'A']})
# create a monomer t for time that will be produced at a rate of 0.2 per hour
Monomer('t')
Monomer('signal')


# Parameter('activate_p53_k', 10)             # per hour
# Parameter('recycle_p53_k', )                # recycle rate equal to __ rate, the limiting step between degradation and synthesis
# Parameter('transcribe_Puma_k', )
# Parameter('transcribe_BIM_k', )
# Parameter('transcribe_Bid_k', )
# Parameter('transcribe_Bax_k', )
# Parameter('transcribe_Bak_k', )
# Parameter('mdm2_bind_p53I_kf', )
# Parameter('mdm2_bind_p53I_kr', )
# Parameter('mdm2_bind_p53A_kf', )
# Parameter('mdm2_bind_p53A_kr', )
Parameter('n', 4)                           # Hill coefficient
Parameter('T', 75)                          # the signal concentration at half maximal p53 activation = 0.75(maximum concentration of signal
Parameter('frequency1', 1)                   # frequency of signal sine wave in units of per hour
Parameter('frequency2', 3)                   # frequency of signal sine wave in units of per hour
Parameter('frequency3', 6)                   # frequency of signal sine wave in units of per hour
Parameter('max_signal', 1.2)                   # maximum amount of signal throughout the simulation, also the amplitude for the oscillatory signal nd height of pulsed signals
Parameter('time_delay_expression', 0.7)
Parameter('time_rate', 0.2)
Expression('simulated_conc_unit', (3/4)*max_signal)

Rule('time', None >> t(), time_rate)         #time points every 0.2 of an hour

Observable('Obs_p53', p53(bf=None, state='A'))
Observable('Obs_signal', signal)
Observable('Obs_time', t)

"""________________________________________________________________________________________________________________________"""

"""Signal Component for ODEs and Signal Input Functions"""

signal_expression = [Expression('signal1', Piecewise((max_signal, (1 <= Obs_time) & (Obs_time <= 2)),
                                                     (max_signal, (6 <= Obs_time) & (Obs_time <= 7)),
                                                     (max_signal, (11 <= Obs_time)& (Obs_time <= 12)),
                                                     (max_signal, (16 <= Obs_time)& (Obs_time <= 17)),
                                                     (max_signal, (21 <= Obs_time)& (Obs_time <= 22)),
                                                     (0, True))),
                     Expression('signal2', Piecewise((max_signal, (1 <= Obs_time) & (Obs_time <= 4)),
                                                     (max_signal, (8 <= Obs_time) & (Obs_time <= 11)),
                                                     (max_signal, (15 <= Obs_time)& (Obs_time <= 19)),
                                                     (0, True))),
                     Expression('signal3', max_signal),
                     Expression('signal4', max_signal*sin(frequency1*Obs_time)),
                     Expression('signal5', max_signal*sin(frequency2*Obs_time)),
                     Expression('signal6', max_signal*sin(frequency3*Obs_time)),
                     Expression('signal7', max_signal * 2.71828**(-1.3863*Obs_time)),
                     Expression('signal8', max_signal * 2.71828**(-0.693147*Obs_time)),
                     Expression('signal9', max_signal * 2.71828**(-0.11552*Obs_time))]


Expression('signal_factor', Obs_signal**n/(Obs_signal**n + T**n))


# def p53_mdm2_loop():
#     Rule('activate_p53', p53(bf=None, state='I') >> p53(bf=None, state='A'), activate_p53_k*signal_factor)
#     # p53 is recycled back into iniactive p53
#     catalyze('mdm2_inactivate_p53', mdm2(bf=None, state='A'), p53(bf=None, state='A') >> p53(bf=None, state='I'), )
#     # mdm2 binds inactive p53 and prevents it from being activated
#     Rule('mdm2_bind_p53I', mdm2(bf=None, state='A') + p53(bf=None, state='I') | p53(bf=1, state='I') % mdm2(bf=1, state='A'), mdm2_bind_p53I_kf, mdm2_bind_p53I_kr)
#     # mdm2 binds active p53 and prevents it from transcribing
#     Rule('mdm2_bind_p53I', mdm2(bf=None, state='A') + p53(bf=None, state='A') | p53(bf=1, state='A') % mdm2(bf=1, state='A'), mdm2_bind_p53A_kf, mdm2_bind_p53A_kr)


for i, x in enumerate(signal_expression):       # this will assign an index number to each of the signals in the signal_expression array
    print  '%s is %d' % (x, i)

for x in signal_expression:                    # loop through the different signals in the signal_expression array and save all the data to a separate array
    #Observable('signal_%d' % (i), signal())
    signal_%d % (i) = [x]
    # p53_mdm2_loop()
    alias_model_components()