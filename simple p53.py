

from pysb import *
from pysb.util import alias_model_components

from earm import shared
from earm import lopez_modules
from earm import albeck_modules
from sympy import sin
from sympy import Piecewise
from pysb.macros import catalyze

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

Parameter('activate_p53_k', 10)             # per hour
Parameter('recycle_p53_k', )                # recycle rate equal to __ rate, the limiting step between degradation and synthesis
Parameter('transcribe_Puma_k', )
Parameter('transcribe_BIM_k', )
Parameter('transcribe_Bid_k', )
Parameter('transcribe_Bax_k', )
Parameter('transcribe_Bak_k', )
Parameter('mdm2_bind_p53I_kf', )
Parameter('mdm2_bind_p53I_kr', )
Parameter('mdm2_bind_p53A_kf', )
Parameter('mdm2_bind_p53A_kr', )
Parameter('n', 4)                           # Hill coefficient
Parameter('T', 75)                          # the signal concentration at half maximal p53 activation = 0.75(maximum concentration of signal
Parameter('frequency1', 1)                   # frequency of signal sine wave in units of per hour
Parameter('frequency2', 3)                   # frequency of signal sine wave in units of per hour
Parameter('frequency3', 6)                   # frequency of signal sine wave in units of per hour
Parameter('max_signal', 1.2)                   # maximum amount of signal throughout the simulation, also the amplitude for the oscillatory signal nd height of pulsed signals
Parameter('time_delay_expression', 0.7)
Parameter('simulated_conc_unit', Expression('simulated_conc_unit', 0.75*max_signal))

Observable('Obs_p53', p53(bf=None, state='A'))
Observable('Obs_signal', signal)
Observable('Obs_time', t)
Observable('time_delay', time_delay_func)

"""________________________________________________________________________________________________________________________"""

"""Signal Component for ODEs and Signal Input Functions"""

Rule('time', None >> t, 0.2)         #time points every 0.2 of an hour

Expression('signal_factor', Obs_signal^n/(Obs_signal^n + T^n))


# signal duration of 1 hour given every 4 hours, with the first signal starting at 1 hour
Expression('signal', Piecewise((0, 0 <= Obs_time < 1)( max_signal, 1 < Obs_time <= 2), (0, 2 < Obs_time < 6),( max_signal, 6 <= Obs_time <= 7), (0, 7 < Obs_time < 11), (max_signal, 11 <= Obs_time <= 12), (0, 12 < Obs_time < 16), (max_signal, 16 <= Obs_time <= 17), (0, 17 < Obs_time < 21), (max_signal, 21 <= Obs_time <= 22), (0, 22 < Obs_time <= 24)))

# maximum signal is given for a duration of 3 hours with four hours in between signals starting at 1 hour
Expression('signal', Piecewise((0, 0 <= Obs_time < 1)( max_signal, 1 < Obs_time <= 4), (0, 4 < Obs_time < 8),( max_signal, 8 <= Obs_time <= 11), (0, 11 < Obs_time < 15), (max_signal, 15 <= Obs_time <= 19), (0, 19 < Obs_time <= 24)))

# Maximum signal for the duration of the simulation
Expression('signal', max_signal)

# sin function for oscillatory signal

Expression('signal', max_signal*sin(frequency1*Obs_time))
Expression('signal', max_signal*sin(frequency2*Obs_time))
Expression('signal', max_signal*sin(frequency3*Obs_time))

# exponential decay of drug signal

Expression('signal', max_signal * 2.71828^(-1.3863*Obs_time))              # for a drug half life of 30 mins
Expression('signal', max_signal * 2.71828^(-0.693147*Obs_time))              # for a drug half life of 60 mins
Expression('signal', max_signal * 2.71828^(-0.11552*Obs_time))              # for a drug half life of 6 hours




"""__________________________________________________________________________________________________________________________"""

# p53 is activated by a signal and then starts transcribing PUMA, Bid, and BIM
# to represent transcription a pool of inactive proteins will be activated
# the system is assumed to be well mixed

Expression('time_delay_func', Obs_time + time_delay_expression)

def p53_mdm2_loop():
    Rule('activate_p53', p53(bf=None, state='I') >> p53(bf=None, state='A'), activate_p53_k*signal_factor)
    # p53 is recycled back into iniactive p53
    catalyze('mdm2_inactivate_p53', mdm2(bf=None, state='A'), p53(bf=None, state='A') >> p53(bf=None, state='I'), )
    # mdm2 binds inactive p53 and prevents it from being activated
    Rule('mdm2_bind_p53I', mdm2(bf=None, state='A') + p53(bf=None, state='I') | p53(bf=1, state='I') % mdm2(bf=1, state='A'), mdm2_bind_p53I_kf, mdm2_bind_p53I_kr)
    # mdm2 binds active p53 and prevents it from transcribing
    Rule('mdm2_bind_p53I', mdm2(bf=None, state='A') + p53(bf=None, state='A') | p53(bf=1, state='A') % mdm2(bf=1, state='A'), mdm2_bind_p53A_kf, mdm2_bind_p53A_kr)


def transcribe():
    # p53 activates PUMA by transcribing it
    Rule('transcribe_Puma', Puma(bf=None, state='I') >> Puma(bf=None, state='A'), transcribe_Puma_k*time_delay)
    # p53 activates BIM by transcribing it
    Rule('transcribe_BIM', BIM(bf=None, state='I') >> BIM(bf=None, state='A'), transcribe_BIM_k*time delay)
    # p53 activates Bid by transcribing it
    Rule('transcribe_Bid', None >> Bid(bf=None, state='U'), transcribe_Bid_k*time_delay)
    # p53 is activated and starts transcribing Bak, and Bax
    Rule('transcribe_Bax', None >> Bax(bf=None, s1=None, s2=None, state='C'), transcribe_Bax_k*Obs_p53)
    Rule('transcribe_Bak', None >> Bak(bf=None, s1=None, s2=None, state='M'), transcribe_Bak_k*Obs_p53)


def BIM_Puma_activate_Bax_and_Bak():
    Rule('BIM_bind_Bax', BIM(bf=None, state='A') + Bax(bf=None, state='C') | Bim(bf=1, state='A') % Bax(bf=1, s1=None, s2=None, state='C'))
    Rule('Puma_bind_Bax', Puma(bf=None, state='A') + Bax(bf=None, state='C') | Puma(bf=1, state='A') % Bax(bf=1, s1=None, s2=None, state='C'))
    Rule('BIM_activate_BaxC', Bim(bf=1, state='A') % Bax(bf=1, state='C') >> Bim(bf=None, state='A') + Bax(bf=None, s1=None, s2=None, state='M'))
    Rule('Puma_activate_BaxC', Puma(bf=1, state='A') % Bax(bf=1, state='C') >> Puma(bf=None, state='A') + Bax(bf=None, s1=None, s2=None state='M'))
    

