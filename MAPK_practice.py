from pysb import *
Model()

Monomer('Raf', ['s','k'], {'s': ['u','p']})
# this creates a monomer named Raf with two sites, k is the catalytic kinase domain and s is a phosphorylation site
#the site s has two states, u is for unphosphorylated, and p is for phosphorylated

Monomer('MEK', ['s218', 's222', 'k'], {'s218': ['u','p'], 's222': ['u', 'p']})
#this creates another monomer for MEK with two phosphorylation sites and a catalytic site, each phosphorylation site has
#two states, one phosphorylated and one unpohosphorylated, the sites are at serine residues 218 and 222

#now we will enter in different parameters, which can be compartment volumes, rate constants, initial conditions or boundaries, etc
#the only two attributes to a parameter are the name and a numerical value
#default value is 0
#can access the parameters by typing model.parameters into the command line
#you need to keep unit consistency on your own because here they are undefined

Parameter('kf', 1e-5)
#here we set the forward reaction rate equal of MEK and Raf binding
Parameter('Raf_0', 7e4)
#here we set an initial concentration for Raf
Parameter('MEK_0', 3e6)
#set an initial protein concentration for MEK

#now we can set up rules
#rules have a name, specify what molecules are reacting, specify how the reactants are transformed into products, and specify rate constants
#rule constructor includes a rule expression (basically just the reaction written out), and two parameters for the rate constants

