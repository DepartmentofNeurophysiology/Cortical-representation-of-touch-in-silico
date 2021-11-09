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

Nconns = 1000  # just to develop

###############################################################################
## Creating populations #######################################################
###############################################################################
Ncells = Cellinfo_All.shape[0]
Npops = int(max(np.unique(Cellinfo_All[0:Ncells,3])))
Pops = {'1':{'Label':'Pyr_SP',   'Layer':'l4',  'Type':'exc'},
        '2':{'Label':'Pyr_SS',   'Layer':'l4',  'Type':'exc'},
        '3':{'Label':'Inh_FS',   'Layer':'l4',  'Type':'inh'},
        '4':{'Label':'Inh_RSNP', 'Layer':'l4',  'Type':'inh'},
        '5':{'Label':'Pyr',      'Layer':'l23', 'Type':'exc'},
        '6':{'Label':'Inh_FSBS', 'Layer':'l23', 'Type':'inh'},   # PV
        '7':{'Label':'Inh_FSCH', 'Layer':'l23', 'Type':'inh'},   # PV
        '8':{'Label':'Inh_BSPV', 'Layer':'l23', 'Type':'inh'},   # PV
        '9':{'Label':'Inh_Mar',  'Layer':'l23', 'Type':'inh'},   # SOM
        '10':{'Label':'Inh_Bit', 'Layer':'l23', 'Type':'inh'},   # SOM
        '11':{'Label':'Inh_DBC', 'Layer':'l23', 'Type':'inh'},   # VIP
        '12':{'Label':'Inh_Bip', 'Layer':'l23', 'Type':'inh'},   # VIP
        '13':{'Label':'Inh_Bip', 'Layer':'l23', 'Type':'inh'},   # CR
        '14':{'Label':'Inh_SBC', 'Layer':'l23', 'Type':'inh'},   # CR
        '15':{'Label':'Inh_NG',  'Layer':'l23', 'Type':'inh'}}   # AC

N = [[n for n in range(Ncells) if Cellinfo_All[n,3]==(ntype+1)] for ntype in range(Npops)]

for ntype in range(len(Pops)):
    Pops[str(ntype+1)]['ids'] = N[ntype]  #ids starting from 0 - differs in 1 from rows in matlab
    if len(N[ntype])==NtAll[ntype]:
        Pops[str(ntype+1)]['Ncells'] = len(N[ntype])
    else:
        print('Check number of cells per population - Pop',ntype+1)
        
    # reading position information, in NetPyNE paradigm (y-axis as depth)
    Pops[str(ntype+1)]['x'] = Cellinfo_All[N[ntype],0]
    Pops[str(ntype+1)]['y'] = Cellinfo_All[N[ntype],2]  # here, y is the depth
    Pops[str(ntype+1)]['z'] = Cellinfo_All[N[ntype],1]

    # reading parameters related to individual intrinsic dynamics (Izhikevich model)
    Pops[str(ntype+1)]['vr'] = vr
    Pops[str(ntype+1)]['u0'] = u0
    Pops[str(ntype+1)]['vt0'] = Neuron_Vt.avg[ntype]
    Pops[str(ntype+1)]['a'] = Neuron_Para[N[ntype],0]
    Pops[str(ntype+1)]['b'] = Neuron_Para[N[ntype],1]
    Pops[str(ntype+1)]['c'] = Neuron_Para[N[ntype],2]
    Pops[str(ntype+1)]['d'] = Neuron_Para[N[ntype],3]
    Pops[str(ntype+1)]['alpha'] = Neuron_Vt.alpha[ntype]
    Pops[str(ntype+1)]['Vi'] = Neuron_Vt.Vi[ntype]
    Pops[str(ntype+1)]['Vmin'] = Neuron_Vt.Vmin[ntype]
    Pops[str(ntype+1)]['Ka'] = Neuron_Vt.Ka[ntype]
    Pops[str(ntype+1)]['Ki'] = Neuron_Vt.Ki[ntype]
    Pops[str(ntype+1)]['tau'] = Neuron_Vt.tau[ntype]
    
    Pops[str(ntype+1)]['list'] = [{'x': Pops[str(ntype+1)]['x'][n],
                                   'y': Pops[str(ntype+1)]['y'][n],
                                   'z': Pops[str(ntype+1)]['z'][n],
                                   'params':{'vr': Pops[str(ntype+1)]['vr'],
                                             'u0': Pops[str(ntype+1)]['u0'],
                                             'vt0': Pops[str(ntype+1)]['vt0'],
                                             'a': Pops[str(ntype+1)]['a'][n],
                                             'b': Pops[str(ntype+1)]['b'][n],
                                             'c': Pops[str(ntype+1)]['c'][n],
                                             'd': Pops[str(ntype+1)]['d'][n],
                                             'alpha': Pops[str(ntype+1)]['alpha'],
                                             'Vi': Pops[str(ntype+1)]['Vi'],
                                             'Vmin': Pops[str(ntype+1)]['Vmin'],
                                             'Ka': Pops[str(ntype+1)]['Ka'],
                                             'Ki': Pops[str(ntype+1)]['Ki'],
                                             'tau': Pops[str(ntype+1)]['tau']}
                                   } for n in range(len(N[ntype]))]

# Network parameters
netParams = specs.NetParams()  # object of class NetParams to store the network parameters

## Cell parameters/rules
GenericCell = {'secs': {}}
GenericCell['secs']['soma'] = {'geom': {}, 'pointps': {}}                  # soma params dict
GenericCell['secs']['soma']['geom'] = {'diam': 6.366, 'L': 5.0, 'cm': 1.0} # Area of 100 um2 --> point process current in [mA/cm2]
GenericCell['secs']['soma']['pointps']['Izhi'] = {                         # soma Izhikevich properties
    'mod': 'Izhi2007b_dyn_thr',
    'C': 1,
    'k': 0.04,
    'vpeak': 10.0,
    'celltype': 1}
netParams.cellParams['IzhiCell'] = GenericCell

for ntype in range(len(Pops)):
    name = Pops[str(ntype+1)]['Label'] + '-' + Pops[str(ntype+1)]['Layer']
    netParams.popParams[name] = {'cellType': 'IzhiCell', 'cellsList': Pops[str(ntype+1)]['list']}

# Assembling conns
netParams.synMechParams['exc'] = {'mod': 'FluctExp2Syn', 'tau_rise': 1.0, 'tau_fall': 2.0, 'cn': 4.0, 'type': 1}
netParams.synMechParams['inh'] = {'mod': 'FluctExp2Syn', 'tau_rise': 1.0, 'tau_fall': 2.0, 'cn': 4.0, 'type': -1}

# Just an input instead of the thalamic inputs
netParams.stimSourceParams['Input'] = {'type': 'NetStim', 'rate': 1000, 'noise': 0.3} # in mA/cm2
for ntype in range(len(Pops)):
    name = Pops[str(ntype+1)]['Label'] + '-' + Pops[str(ntype+1)]['Layer']
    netParams.stimTargetParams['Input->'+name] = {'source': 'Input', 'sec':'soma', 'loc': 0.5, 'conds': {'pop':name}, 'weight': 0.005, 'synMech': 'exc'}


# Connections List: [Pre, Post]
# List of excitatory connections
connListExc = [ [PMat_AlltoAll.preCell[nc][1],PMat_AlltoAll.preCell[nc][0]] for nc in range(len(PMat_AlltoAll.preCell)) if PMat_AlltoAll.preCell[nc][2]==1]
connListExc_gID = [ [PMat_AlltoAll.preCell[nc][1]-1,PMat_AlltoAll.preCell[nc][0]-1] for nc in range(len(PMat_AlltoAll.preCell)) if PMat_AlltoAll.preCell[nc][2]==1]
# weights (nA or equiv. mA/cm2 if thought of as distributed)
weightListExc = [ 0.001*PMat_AlltoAll.Am[nc] for nc in range(len(PMat_AlltoAll.preCell)) if PMat_AlltoAll.preCell[nc][2]==1]
CVListExc =     [ 0.001*PMat_AlltoAll.CV[nc] for nc in range(len(PMat_AlltoAll.preCell)) if PMat_AlltoAll.preCell[nc][2]==1]
#delays
delayListExc = [ PMat_AlltoAll.Delay[nc] for nc in range(len(PMat_AlltoAll.preCell)) if PMat_AlltoAll.preCell[nc][2]==1]
# reduced (just for dev.)
connListExc_gID2 = connListExc_gID[:Nconns]
weightListExc2 = weightListExc[:Nconns]
CVListExc2 = CVListExc[:Nconns]
delayListExc2 = delayListExc[:Nconns]


# List of inhibitory connections
connListInh = [ [PMat_AlltoAll.preCell[nc][1],PMat_AlltoAll.preCell[nc][0]] for nc in range(len(PMat_AlltoAll.preCell)) if PMat_AlltoAll.preCell[nc][2]==-1]
connListInh_gID = [ [PMat_AlltoAll.preCell[nc][1]-1,PMat_AlltoAll.preCell[nc][0]-1] for nc in range(len(PMat_AlltoAll.preCell)) if PMat_AlltoAll.preCell[nc][2]==-1]
# weights (nA or equiv. mA/cm2 if thought of as distributed)
weightListInh = [ 0.001*PMat_AlltoAll.Am[nc] for nc in range(len(PMat_AlltoAll.preCell)) if PMat_AlltoAll.preCell[nc][2]==-1]
CVListInh =     [ 0.001*PMat_AlltoAll.CV[nc] for nc in range(len(PMat_AlltoAll.preCell)) if PMat_AlltoAll.preCell[nc][2]==-1]
#delays
delayListInh = [ PMat_AlltoAll.Delay[nc] for nc in range(len(PMat_AlltoAll.preCell)) if PMat_AlltoAll.preCell[nc][2]==-1]
# reduced (just for dev.)
connListInh_gID2 = connListInh_gID[:Nconns]
weightListInh2 = weightListInh[:Nconns]
CVListInh2 = CVListInh[:Nconns]
delayListInh2 = delayListInh[:Nconns]


netParams.connParams['All->All_exc'] = {
        'preConds': {'cellType': 'IzhiCell'},   # conditions of presyn cells
        'postConds': {'cellType': 'IzhiCell'},  # conditions of postsyn cells
        'connList': connListExc_gID2,            # list of conns
        'weight': weightListExc2,                # synaptic weight
        'delay': delayListExc2,                  # transmission delay (ms)
        'synMech': 'exc'}

netParams.connParams['All->All_inh'] = {
        'preConds': {'cellType': 'IzhiCell'},   # conditions of presyn cells
        'postConds': {'cellType': 'IzhiCell'},  # conditions of postsyn cells
        'connList': connListInh_gID2,            # list of conns
        'weight': weightListInh2,                # synaptic weight
        'delay': delayListInh2,                  # transmission delay (ms)
        'synMech': 'inh'}

# Simulation options
simConfig = specs.SimConfig()       # object of class SimConfig to store simulation configuration
simConfig.hParams['v_init'] = -65.0
simConfig.duration = 500          # Duration of the simulation, in ms
simConfig.dt = 0.025                # Internal integration timestep to use
simConfig.verboseLevel = 1           # Show detailed messages
simConfig.recordTraces = {'V_soma1':{'sec':'soma','loc':0.5,'var':'v'}}  # Dict with traces to record
simConfig.recordStep = 0.025            # Step size in ms to save data (eg. V traces, LFP, etc)
simConfig.filename = 'sim1'         # Set file output name
simConfig.saveDat = True 

simConfig.analysis['plotTraces'] = {'include': [0,1,2,3,4], 'saveFig': True}  # Plot recorded traces for this list of cells
simConfig.analysis['plotRaster'] = {'include': [0,1,2,3,4], 'saveFig': True}                  # Plot a raster

# Create the network with a generic synaptic mechanism
sim.createSimulateAnalyze(netParams = netParams, simConfig = simConfig)

# Modify connections (synaptic properties of each connection)
for nconn in range(len(connListExc_gID2)):
    Pre_Gid = connListExc_gID2[nconn][0]
    Post_Gid = connListExc_gID2[nconn][1]
    #print('Nconn: ',nconn,' - Pre: ',Pre_Gid,' - Post: ',Post_Gid)
    list_preGids = [sim.net.cells[Post_Gid].conns[nc]['preGid'] for nc in range(len(sim.net.cells[Post_Gid].conns))]
    try:
        conn_index = list_preGids.index(Pre_Gid)
    except:
        print("Error in connections pre-post Ids")
        
    # Setting individual parameters
    Pop_pre = [nc for nc in range(len(N)) if (Pre_Gid in N[nc])==True][0]
    sim.net.cells[Post_Gid].conns[conn_index]['hObj'].syn().tau_rise = PMat_AlltoAll.Trise[Post_Gid][Pop_pre]
    sim.net.cells[Post_Gid].conns[conn_index]['hObj'].syn().tau_fall = PMat_AlltoAll.Tfall[Post_Gid][Pop_pre]
    sim.net.cells[Post_Gid].conns[conn_index]['hObj'].syn().cn = PMat_AlltoAll.CN[Post_Gid][Pop_pre]
    sim.net.cells[Post_Gid].conns[conn_index]['hObj'].syn().mean_amp = weightListExc2[nconn]
    sim.net.cells[Post_Gid].conns[conn_index]['hObj'].syn().cv = CVListExc2[nconn]
    
for nconn in range(len(connListInh_gID2)):
    Pre_Gid = connListInh_gID2[nconn][0]
    Post_Gid = connListInh_gID2[nconn][1]
    #print('Nconn: ',nconn,' - Pre: ',Pre_Gid,' - Post: ',Post_Gid)
    list_preGids = [sim.net.cells[Post_Gid].conns[nc]['preGid'] for nc in range(len(sim.net.cells[Post_Gid].conns))]
    try:
        conn_index = list_preGids.index(Pre_Gid)
    except:
        print("Error in connections pre-post Ids")
        
    # Setting individual parameters
    Pop_pre = [nc for nc in range(len(N)) if (Pre_Gid in N[nc])==True][0]
    sim.net.cells[Post_Gid].conns[conn_index]['hObj'].syn().tau_rise = PMat_AlltoAll.Trise[Post_Gid][Pop_pre]
    sim.net.cells[Post_Gid].conns[conn_index]['hObj'].syn().tau_fall = PMat_AlltoAll.Tfall[Post_Gid][Pop_pre]
    sim.net.cells[Post_Gid].conns[conn_index]['hObj'].syn().cn = PMat_AlltoAll.CN[Post_Gid][Pop_pre]
    sim.net.cells[Post_Gid].conns[conn_index]['hObj'].syn().mean_amp = weightListInh2[nconn]
    sim.net.cells[Post_Gid].conns[conn_index]['hObj'].syn().cv = CVListInh2[nconn]
    
# adding flag to watch
conn_to_watch = 0
Pre_Gid = connListExc_gID2[conn_to_watch][0]
Post_Gid = connListExc_gID2[conn_to_watch][1]
list_preGids = [sim.net.cells[Post_Gid].conns[nc]['preGid'] for nc in range(len(sim.net.cells[Post_Gid].conns))]
try:
    conn_index = list_preGids.index(Pre_Gid)
except:
    print("Error in connections pre-post Ids")

sim.net.cells[Post_Gid].conns[conn_index]['hObj'].syn().flag_print = 1
    
# simulate again
sim.cfg.filename = 'sim2'
sim.simulate()
sim.analyze()
