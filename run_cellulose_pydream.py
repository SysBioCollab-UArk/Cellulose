from cellulose import model
from param_calibration import *
from pysb.simulator import ScipyOdeSimulator

solver = ScipyOdeSimulator(model)
sim_protocol = SimulationProtocol(solver)

custom_priors = {'Cellulose_0': ('uniform', 0.3)}

obs_labels = {'Furfural_Obs': 'Furfural', 'HMF_Obs': 'HMF', 'LA_Obs': 'LA'}
exp_data_file = os.path.join('.', 'cellulose_data_160C.csv')

if __name__ == '__main__':

    calibrator = ParameterCalibration(model, exp_data_file, sim_protocol, priors=custom_priors)

    calibrator.run(niterations=50000, nchains=5, obs_labels=obs_labels, plot_results=True,
                   plot_tc_args={'separate_plots': False})
