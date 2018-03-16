from pysb import *

Model()                                      #create model
Monomer('A')                                 #enter a single monomer called A
Parameter('k',3.0)                           #enter a rate parameter of 3
Rule('synthesize_A', None >> A(), k)         #create a rule that synthesizes the monomer protein at a rate of 3 copies per second

# the above are model definitions, they aren't for specific model usage.
# we will now define time paramteres and a solver and pysb.integrate to set up how we want to run the model we just defined
# every model component automatically assigns itself a variable named the components name, so the model will create a monomer A
# the variable model holds the Model object
# every model must begin with from pysb import * and also the Model()
# first line brings in the python classes needed to define the model, so it lets pysb do its thing
#the second line creates an instance of the model and implicitly assigns the object to the variable model
#Required component declarations: Monomer, Parameter, and Rule
#optional component declarations: observable, and compartment,
# each component is represented by a python class with base class Component
#the component base class sets up a name for the subclass and then a self export functionality that assigns a variable named after the argument
#monomers have a name and a list of sites
#sites can bind other sites or take on a state
#the site list is optional, but most models will need this specified
#after defining the monomer name and a list of all the site names, define a dictionary for the allowable states for the sites
#


from _future_ import print_function
from pysb.simulator import ScipyOdeSimulator
from testing import model



