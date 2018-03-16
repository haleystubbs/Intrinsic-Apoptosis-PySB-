# import the pysb module and all its methods and functions
from pysb import *
import matplotlib.pyplot as plt
from pysb.bng import *
from pysb.simulator import ScipyOdeSimulator
import pylab as pl

import numpy as np
# instantiate a model

Model()


"""declare monomers"""
Monomer('Bid', ['b', 'S'], {'S': ['u', 't']})
# a protein Bid has two sites, a binding site and a truncation site
# the truncation site has two states: unbound (whole) or truncated
Monomer('C8', ['b'])
# the monomer C8 only has one site, which is the binding site


"""Input the Parameter Values"""
Parameter('kf', 1.0e-07)
Parameter('kr', 1.0e-03)
Parameter('kc', 1.0)


"""Input the Rules"""
Rule('C8_Bid_bind', C8(b=None) + Bid(b=None, S='u') | C8(b=1) % Bid(b=1, S='u'), kf, kr)
# we set up a rule for the reversible binding of C8 to Bid, the + indicates complexation, | means its reversible,
# b is the number of bound units
Rule('tBid_from_C8Bid', C8(b=1) % Bid(b=1, S='u') >> C8(b=None) + Bid(b=None, S='t'), kc)
# here we set up a reaction to represent the bound Bid and C8 complex catalyzing the truncation of Bid


# need to place bounds on the system, so that the ODEs can be solved by integration
# need to set more parameter values and then set them as initial conditions


"""Set up Initial Conditions"""
Parameter('C8_0', 1000)
# set the inital number of C8 molecules at time 0 secs equal to 1000
Parameter('Bid_0', 10000)
Initial(C8(b=None), C8_0)
Initial(Bid(b=None, S='u'), Bid_0)


"""Enter Observables"""
Observable('obsC8', C8(b=None))
Observable('obsBid', Bid(b=None, S='u'))
Observable('obstBid', Bid(b=None, S='t'))

generate_network(model)
generate_equations(model)

t = pl.linspace(0,20000)
simres = ScipyOdeSimulator(model, tspan=t).run()
yout = simres.all

yout['obsBid']
print(yout)

pl.ion()
pl.figure()
pl.plot(t, yout['obsBid'], label="Bid")
pl.plot(t, yout['obstBid'], label="tBid")
pl.plot(t, yout['obsC8'], label="C8")
pl.legend()
pl.xlabel("Time (s)")
pl.ylabel("Molecules/cell")
pl.show()

print(pl)
