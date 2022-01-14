###############################################################################
## Setting-up, initial conditions, options ####################################
###############################################################################

def setup():
    InputOptions = ['Multitrial', 'Svoboda']
    # Multitrial uses PSTH data collected by Aguilar & Castro-Alamancos  ("Spatiotemporal Gating of Sensory Inputs in Thalamus during Quiescent and Activated States", J. Neurosci. 2005)
    # Svoboda loads spike trains generated by the matlab model, based on whisker recordings (Svoboda dataset) and barreloids' filtering models
    Input = {'Option' : 'Multitrial',
             'Folder' : 'ReadingData_Aguilar',
             'Filename' : 'psth.dat'}
#    Input = {'Option' : 'Svoboda',
#             'Folder' : 'InstantiatedModel',
#             'Filename' : 'Test_sim_Svoboda_an171923_19-Sep-2021_Thalamic_Spike_Trains.mat'}
    
    Model = {'Folder' : 'InstantiatedModel',
             'PreModel' : 'CMDMs_Test_sim_Svoboda.mat',
             'FullModel' : 'CMDMs_Test_sim_Svoboda_ConData.mat'}
    
    Ntrials = 25          # number of trials (different inputs)
    Nrepetitions = 2     # repetitions of an unique input condition
    ExperimConds = {'Ntrials' : Ntrials,
                    'Nrepetitions' : Nrepetitions}
    
    vr = -60  # [-60, -70, -80] resting potential, also sets the initial condition
    v0 = -70
    u0 = 0    # recovery variable

    dyn_thres = 1     # 0: normal izhikevich model - 1: adaptative threshold (Huang et al, "Adaptive Spike Threshold Enables Robust and Temporally Precise Neuronal Encoding", PLoS Comp. Biol. 2016)
    tau_plas = 120.0  # plasticity - STD
    fr = 1.00         # Fraction of recurrent connections - Useful for setting-up sims

    Settings = {'vr' : vr, 
                'v0' : v0,
                'u0' : u0,
                'dyn_thres' : dyn_thres,
                'tau_plas' : tau_plas,
                'fr' : fr}

    return Input, Model, ExperimConds, Settings
