#!/usr/bin/env Python3

import gillespy2

def Dimerization(parameter_values=None):
    # First call the gillespy2.Model initializer.
    model = gillespy2.Model()

    # Define parameters for the rates of creation and dissociation.
    k_c = gillespy2.Parameter(name='k_c', expression=0.005)
    k_d = gillespy2.Parameter(name='k_d', expression=0.08)
    model.add_parameter([k_c, k_d])

    # Define variables for the molecular species representing M & D.
    m = gillespy2.Species(name='monomer', initial_value=30)
    d = gillespy2.Species(name='dimer',   initial_value=0)
    model.add_species([m, d])

    # The list of reactants and products for a Reaction object are
    # each a Python dictionary in which the dictionary keys are
    # Species objects and the values are stoichiometries of the
    # species in the reaction.
    r_c = gillespy2.Reaction(name="r_creation", rate=k_c,
                             reactants={m:2}, products={d:1})
    r_d = gillespy2.Reaction(name="r_dissociation", rate=k_d,
                             reactants={d:1}, products={m:2})
    model.add_reaction([r_c, r_d])

    # Set the timespan for the simulation.
    tspan = gillespy2.TimeSpan.linspace(t=100, num_points=101)
    model.timespan(tspan)
    return model

model = Dimerization()
results = model.run(number_of_trajectories=10)

import matplotlib.pyplot as plt

for index in range(0, 10):
    trajectory = results[index]
    plt.plot(trajectory['time'], trajectory['monomer'], 'r')
    plt.plot(trajectory['time'], trajectory['dimer'],   'b')
    

# Define the model.

def AutomaticSwitchExample(parameter_values=None):
    # First call the gillespy2.Model initializer.
    model = gillespy2.Model(name="Automatic Switch Example")

    # Define parameters.
    k1 = gillespy2.Parameter(name='k1', expression=3e-4)
    k2 = gillespy2.Parameter(name='k2', expression=0.5e-2)
    k3 = gillespy2.Parameter(name='k3', expression=2e-1)
    model.add_parameter([k1,k2,k3])

    # Define species.
    A = gillespy2.Species(name='A', initial_value=400)
    B = gillespy2.Species(name='B', initial_value=10000)
    C = gillespy2.Species(name='C', initial_value=10000)
    model.add_species([A, B, C])

    # Define reactions.
    r1 = gillespy2.Reaction(name="r1", rate=k1,
                            reactants={A:1,B:1}, products={B:1,C:1})

    r2 = gillespy2.Reaction(name="r2", rate=k2,
                            reactants={B:1}, products={})

    r3 = gillespy2.Reaction(name="r3", rate=k3,
                            reactants={C:1}, products={A:1})

    model.add_reaction([r1,r2,r3])

    # Set the timespan for the simulation.
    tspan = gillespy2.TimeSpan.linspace(t=600, num_points=601)
    model.timespan(tspan)
    return model

# Create an instance of the model object, then run the simulation.

model = AutomaticSwitchExample()
results = model.run(algorithm="Tau-Hybrid")

plt.figure(figsize=(15, 10))
for species in results[0]:
    if species == 'time':
        continue
    plt.plot(results[0]['time'], results[0][species],
             label='{0}'.format(species))
plt.title('Example Hybrid Switching Model')
plt.legend(loc='best') 
    
    

