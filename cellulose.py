from pysb import *
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt

Model()

# Monomers
Monomer('Cellulose')
Monomer('Glucose')
Monomer('Fructose')
Monomer('HMF')
Monomer('Furfural')
Monomer('FA')
Monomer('LA')

# Initial conditions
Parameter('Cellulose_0', 50) # 0.5
Initial(Cellulose(), Cellulose_0)

# Rules
Parameter('k1', 0.01) # 0.211)
# Rule('Cellulose_to_Glucose', Cellulose() >> Cellulose() + Glucose(), k1)
Rule('Cellulose_to_Glucose', Cellulose() >> Glucose(), k1)

Parameter('k2', 0.01) # 0.143)
Rule('Glucose_to_HMF', Glucose() >> HMF(), k2)

Parameter('k3', 0.01) # 1.651)
Rule('HMF_to_LA_FA', HMF() >> LA() + FA(), k3)

Parameter('k4', 0.1) # 1.870)
Rule('Glucose_to_Furfural_FA', Glucose() >> Furfural() + FA(), k4)

Parameter('k5', 0.001) # 1.533)
Rule('Furfural_to_LA', Furfural() >> LA(), k5)

Parameter('kf6', 0.1) # 0.350)
Parameter('kr6', 0.1) # 0.349)
Rule('Glucose_to_Fructose', Glucose() | Fructose(), kf6, kr6)

Parameter('k7', 0.1) # 1.810)
Rule('Fructose_to_HMF', Fructose() >> HMF, k7)

# Observables
Observable('Glucose_Obs', Glucose())
Observable('Fructose_Obs', Fructose())
Observable('HMF_Obs', HMF())
Observable('Furfural_Obs', Furfural())
Observable('FA_Obs', FA())
Observable('LA_Obs', LA())

if __name__ == '__main__':

    # experimental data
    data = np.genfromtxt('cellulose_data_160C.csv', dtype=None, delimiter=',', names=True, encoding="utf_8_sig")
    exp_observables = np.unique([d['observable'] for d in data])

    # Simulation commands
    tspan = np.linspace(0, 16*60, 16*60*10+1)
    # tspan = np.linspace(0, 150, 1501)
    sim = ScipyOdeSimulator(model, tspan, verbose=True)
    output = sim.run()

    plt.figure(constrained_layout=True, figsize=(6.4*1.2, 4.8))
    for obs in model.observables:
        # plot simulation result
        p = plt.plot(tspan, output.observables[obs.name], lw=2, label=obs.name)
        if obs.name in exp_observables:
            # plot experimental data
            t_pts = [d['time'] * 60 for d in data if d['observable'] == obs.name] # convert from hrs to minutes
            avg = [d['average'] for d in data if d['observable'] == obs.name]
            stderr = [d['stderr'] for d in data if d['observable'] == obs.name]
            plt.errorbar(t_pts, avg, stderr, fmt='o', color=p[0].get_color(), capsize=6)
    plt.xlabel('Time (min)')
    plt.ylabel('Concentration (mM)')
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1))

    # ODEs
    for i, ode in enumerate(model.odes):
        print('%s:' % model.species[i], ode)

    plt.show()
