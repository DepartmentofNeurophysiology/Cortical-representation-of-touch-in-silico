from netpyne import specs

def set_cfg(input,Nbarrels):
    # Simulation options
    simConfig = specs.SimConfig()       # object of class SimConfig to store simulation configuration
    if input['Option'] == 'Multitrial':
        simConfig.duration = 40               # Duration of the simulation, in ms
        simConfig.dt = 0.025                  # Internal integration timestep to use
        simConfig.recordStep = 0.025          # Step size in ms to save data (eg. V traces, LFP, etc)
    elif input['Option'] == 'Svoboda':
        simConfig.duration = 6000             # Duration of the simulation, in ms - set in the Matlab code to generate the SpikeTrain struct
        simConfig.dt = 0.1                    # Internal integration timestep to use
        simConfig.recordStep = 0.1            # Step size in ms to save data (eg. V traces, LFP, etc)
    else:
        print('Sim config set at default conds')
    
    simConfig.recordTraces = {'v_soma':{'sec':'soma','loc':0.5,'var':'v'},
                              'u_soma':{'sec':'soma', 'pointp':'Izhi','var':'u'},
                              'vt_soma':{'sec':'soma', 'pointp':'Izhi','var':'vt'},
                              'synapse_exc_std':{'sec':'soma','loc':0.5,'synMech':'exc','var':'std'}, # STD: just first synapse (called from post-side)
                              'synapse_exc_s':{'sec':'soma','loc':0.5,'synMech':'exc','var':'s'},     # Gating s: just first synapse
                              'synapse_inh_std':{'sec':'soma','loc':0.5,'synMech':'inh','var':'std'}, # STD: just first synapse
                              'synapse_inh_s':{'sec':'soma','loc':0.5,'synMech':'inh','var':'s'}      # Gating s: just first synapse
                              }  # Dict with traces to record


    simConfig.gatherOnlySimData = False
    simConfig.saveDataInclude = ['simData']
    simConfig.saveJson = True

    simConfig.analysis['plotRaster'] = {'include': ['all'], 'saveFig': True}                  # Plot a raster

    # bulk cells in main barrel, one per population
    if Nbarrels == 1:
        simConfig.analysis['plotTraces'] = {'include': [420, 1107, 1621, 1715, 2877, 4046, 4101, 4131, 4190, 4254, 4303, 4350, 4399], 'saveFig': True}  # Plot recorded traces for this list of cells
    elif Nbarrels == 3:
        simConfig.analysis['plotTraces'] = {'include': [1260, 3321, 4863, 5144, 8631, 12138, 12303, 12392, 12570, 12762, 12909, 13050, 13197], 'saveFig': True}  # Plot recorded traces for this list of cells
    
    return simConfig