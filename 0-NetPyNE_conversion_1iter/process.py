import os
import numpy as np
from scipy.io import  loadmat
from netpyne import specs, sim


reading_folder  = 'InstantiatedModel'
datamodel_folder = os.path.join(os.getcwd(),reading_folder)

filename = 'CMDMs_Test_sim_Svoboda_ConData.mat'
reading_filename = os.path.join(datamodel_folder,filename)
inst_model = loadmat(reading_filename, struct_as_record=False, squeeze_me=True)

filename = 'Test_sim_Svoboda_an171923_19-Sep-2021_Thalamic_Spike_Trains.mat'
reading_filename = os.path.join(datamodel_folder,filename)
inst_input = loadmat(reading_filename, struct_as_record=False, squeeze_me=True)

# Data from cells and conns
PMat_IntoAll = inst_model['ConData'].PMat_IntoAll
PMat_AlltoAll = inst_model['ConData'].PMat_AlltoAll
Neuron_Para = inst_model['ConData'].Neuron_Para
Neuron_Vt = inst_model['ConData'].Neuron_Vt
Cellinfo_In = inst_model['ConData'].Cellinfo_In
Cellinfo_All = inst_model['ConData'].Cellinfo_All
NtIn = inst_model['ConData'].NtIn
NtAll = inst_model['ConData'].NtAll

# Data from spike train inputs
SpikeTrainStruct = inst_input['SpikeTrainStruct']
#PSTH1 = SpikeTrainStruct[0].PSTH
#SpikeTimes1 = SpikeTrainStruct[0].SpikeTimes
#SpikeCount1 = SpikeTrainStruct[0].SpikeCount

###############################################################################
## Initial conditions #########################################################
###############################################################################
v0 = -65
vr = -65
u0 = 0

###############################################################################
## Creating populations #######################################################
###############################################################################
Ncells = Cellinfo_All.shape[0]
Npops = int(max(np.unique(Cellinfo_All[0:Ncells,3])))
Pops = {'1':{'Label':'Pyr_SP',  'Layer':'l4',  'Type':'exc', 'Neuron_Vt':{}},
        '2':{'Label':'Pyr_SS',  'Layer':'l4',  'Type':'exc', 'Neuron_Vt':{}},
        '3':{'Label':'Inh_FS',   'Layer':'l4',  'Type':'inh', 'Neuron_Vt':{}},
        '4':{'Label':'Inh_RSNP', 'Layer':'l4',  'Type':'inh', 'Neuron_Vt':{}},
        '5':{'Label':'Pyr',  'Layer':'l23', 'Type':'exc', 'Neuron_Vt':{}},
        '6':{'Label':'Inh_FSBS', 'Layer':'l23', 'Type':'inh', 'Neuron_Vt':{}},   # PV
        '7':{'Label':'Inh_FSCH', 'Layer':'l23', 'Type':'inh', 'Neuron_Vt':{}},   # PV
        '8':{'Label':'Inh_BSPV', 'Layer':'l23', 'Type':'inh', 'Neuron_Vt':{}},   # PV
        '9':{'Label':'Inh_Mar',  'Layer':'l23', 'Type':'inh', 'Neuron_Vt':{}},   # SOM
        '10':{'Label':'Inh_Bit', 'Layer':'l23', 'Type':'inh', 'Neuron_Vt':{}},   # SOM
        '11':{'Label':'Inh_DBC', 'Layer':'l23', 'Type':'inh', 'Neuron_Vt':{}},   # VIP
        '12':{'Label':'Inh_Bip', 'Layer':'l23', 'Type':'inh', 'Neuron_Vt':{}},   # VIP
        '13':{'Label':'Inh_Bip', 'Layer':'l23', 'Type':'inh', 'Neuron_Vt':{}},   # CR
        '14':{'Label':'Inh_SBC', 'Layer':'l23', 'Type':'inh', 'Neuron_Vt':{}},   # CR
        '15':{'Label':'Inh_NG',  'Layer':'l23', 'Type':'inh', 'Neuron_Vt':{}}}   # AC

N = [[n for n in range(Ncells) if Cellinfo_All[n,3]==(ntype+1)] for ntype in range(Npops)]

for ntype in range(len(Pops)):
    Pops[str(ntype+1)]['ids'] = N[ntype]  #ids starting from 0 - differs in 1 from rows in matlab
    if len(N[ntype])==NtAll[ntype]:
        Pops[str(ntype+1)]['Ncells'] = len(N[ntype])
    else:
        print('Check number of cells per population - Pop',ntype+1)
        
    #for field in Neuron_Vt._fieldnames:
    #    Pops[str(ntype+1)]['Neuron_Vt'][field] = getattr(Neuron_Vt, field)[ntype]
    
    # reading position information, in NetPyNE paradigm (y-axis as depth)
    Pops[str(ntype+1)]['x'] = Cellinfo_All[N[ntype],0]
    Pops[str(ntype+1)]['y'] = Cellinfo_All[N[ntype],2]  # here, y is the depth
    Pops[str(ntype+1)]['z'] = Cellinfo_All[N[ntype],1]

    # reading parameters related to individual intrinsic dynamics (Izhikevich model)
    Pops[str(ntype+1)]['vr'] = vr
    Pops[str(ntype+1)]['vt'] = Neuron_Vt.avg[ntype]
    Pops[str(ntype+1)]['a'] = Neuron_Para[N[ntype],0]
    Pops[str(ntype+1)]['b'] = Neuron_Para[N[ntype],1]
    Pops[str(ntype+1)]['c'] = Neuron_Para[N[ntype],2]
    Pops[str(ntype+1)]['d'] = Neuron_Para[N[ntype],3]
    
    Pops[str(ntype+1)]['list'] = [{'x': Pops[str(ntype+1)]['x'][n],
                                   'y': Pops[str(ntype+1)]['y'][n],
                                   'z': Pops[str(ntype+1)]['z'][n],
                                   'params':{'vr': Pops[str(ntype+1)]['vr'],
                                             'vt': Pops[str(ntype+1)]['vt'],
                                             'a': Pops[str(ntype+1)]['a'][n],
                                             'b': Pops[str(ntype+1)]['b'][n],
                                             'c': Pops[str(ntype+1)]['c'][n],
                                             'd': Pops[str(ntype+1)]['d'][n]}
                                   } for n in range(len(N[ntype]))]

# Network parameters
netParams = specs.NetParams()  # object of class NetParams to store the network parameters

## Cell parameters/rules
GenericCell = {'secs': {}}
GenericCell['secs']['soma'] = {'geom': {}, 'pointps': {}}                  # soma params dict
GenericCell['secs']['soma']['geom'] = {'diam': 6.366, 'L': 5.0, 'cm': 1.0} # Area of 100 um2 --> point process current in [mA/cm2]
GenericCell['secs']['soma']['pointps']['Izhi'] = {                         # soma Izhikevich properties
    'mod': 'Izhi2007b',
    'C': 1,
    'k': 0.04,
    'vpeak': 0.0,
    'celltype': 1}
netParams.cellParams['IzhiCell'] = GenericCell

netParams.stimSourceParams['Input'] = {'type': 'IClamp', 'del': 50, 'dur': 500, 'amp': 0.015} # in mA/cm2
for ntype in range(len(Pops)):
    name = Pops[str(ntype+1)]['Label'] + '-' + Pops[str(ntype+1)]['Layer']
    netParams.popParams[name] = {'cellType': 'IzhiCell', 'cellsList': Pops[str(ntype+1)]['list']}
    netParams.stimTargetParams['Input->'+name] = {'source': 'Input', 'sec':'soma', 'loc': 0.5, 'conds': {'pop':name}}


# Simulation options
simConfig = specs.SimConfig()       # object of class SimConfig to store simulation configuration
simConfig.hParams['v_init'] = -65.0
simConfig.duration = 500          # Duration of the simulation, in ms
simConfig.dt = 0.025                # Internal integration timestep to use
simConfig.verbose = False           # Show detailed messages
simConfig.recordTraces = {'V_soma':{'sec':'soma','loc':0.5,'var':'v'}}  # Dict with traces to record
simConfig.recordStep = 0.025            # Step size in ms to save data (eg. V traces, LFP, etc)
simConfig.filename = 'iclamp'         # Set file output name
simConfig.saveDat = True 

simConfig.analysis['plotTraces'] = {'include': [0,1,2], 'saveFig': True}  # Plot recorded traces for this list of cells


# Create network and run simulation
sim.createSimulateAnalyze(netParams = netParams, simConfig = simConfig)
