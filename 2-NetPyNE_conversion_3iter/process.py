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

###############################################################################
## Thalamic Inputs ############################################################
###############################################################################
Nbarrels = len(SpikeTrainStruct)
n , Ntrials = np.shape(SpikeTrainStruct[0].SpikeTimes)
Nthalamic = sum([len(SpikeTrainStruct[n].SpikeTimes) for n in range(Nbarrels)])


ntrial = 1
SpikeTimes_perBarrel = []
for nbarrel in range(Nbarrels):
    SpikeTimes_perBarrel.append([ [SpikeTrainStruct[nbarrel].SpikeTimes[n,ntrial-1]] if isinstance(SpikeTrainStruct[nbarrel].SpikeTimes[n,ntrial-1],int) else SpikeTrainStruct[nbarrel].SpikeTimes[n,ntrial-1].tolist() for n in range(len(SpikeTrainStruct[nbarrel].SpikeTimes)) ])
SpikeTimes = np.concatenate(SpikeTimes_perBarrel).ravel()


###############################################################################
## Initial conditions #########################################################
###############################################################################
v0 = -65 # not used
vr = -65
u0 = 0

tau_plas = 120.0

###############################################################################
## Formatting data ############################################################
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


###############################################################################
## NetPyNE specifications #####################################################
###############################################################################

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

# Population parameters - First, we define thalamic cells, then cortical. This impacts on individual gIDs read from Matlab connections
# Population corresponding to thalamic cells
netParams.popParams['artificial'] = {'cellModel': 'VecStim', 'numCells': Nthalamic, 'spkTimes': SpikeTimes, 'xRange': [-0.01,0.01],'yRange': [0,0.01],'zRange': [-0.01,0.01]}

# Populations in cortex
for ntype in range(len(Pops)):
    name = Pops[str(ntype+1)]['Label'] + '-' + Pops[str(ntype+1)]['Layer']
    netParams.popParams[name] = {'cellType': 'IzhiCell', 'cellsList': Pops[str(ntype+1)]['list']}


# Defining a generic synapse
netParams.synMechParams['exc'] = {'mod': 'FluctExp2Syn', 'tau_rise': 1.0, 'tau_fall': 2.0, 'cn': 4.0, 'type': 1}
netParams.synMechParams['inh'] = {'mod': 'FluctExp2Syn', 'tau_rise': 1.0, 'tau_fall': 2.0, 'cn': 4.0, 'type': -1}


# Assembling conns - Connections List: [Pre, Post], gIDs are relative to the regions (thalamus or cortex). This is consistent to the connParam rule set by "connList"
# Thalamus -> Cortex (IntoAll)
connListExc_ThtoAll = [ [PMat_IntoAll.preCell[nc][1],PMat_IntoAll.preCell[nc][0]] for nc in range(len(PMat_IntoAll.preCell)) if PMat_IntoAll.preCell[nc][2]==1]
connListExc_ThtoAll_gID = [ [PMat_IntoAll.preCell[nc][1]-1,PMat_IntoAll.preCell[nc][0]-1] for nc in range(len(PMat_IntoAll.preCell)) if PMat_IntoAll.preCell[nc][2]==1]
# weights (nA or equiv. mA/cm2 if thought of as distributed)
weightListExc_ThtoAll = [ 0.001*PMat_IntoAll.Am[nc] for nc in range(len(PMat_IntoAll.preCell)) if PMat_IntoAll.preCell[nc][2]==1]
CVListExc_ThtoAll =     [ 0.001*PMat_IntoAll.CV[nc] for nc in range(len(PMat_IntoAll.preCell)) if PMat_IntoAll.preCell[nc][2]==1]
#delays
delayListExc_ThtoAll = [ PMat_IntoAll.Delay[nc] for nc in range(len(PMat_IntoAll.preCell)) if PMat_IntoAll.preCell[nc][2]==1]
# short-term learning & failure
STDListExc_ThtoAll = [ PMat_IntoAll.STD[nc] for nc in range(len(PMat_IntoAll.preCell)) if PMat_IntoAll.preCell[nc][2]==1]
FailListExc_ThtoAll = [ PMat_IntoAll.Fail[nc] for nc in range(len(PMat_IntoAll.preCell)) if PMat_IntoAll.preCell[nc][2]==1]
PlasListExc_ThtoAll = [ PMat_IntoAll.Plas[nc] for nc in range(len(PMat_IntoAll.preCell)) if PMat_IntoAll.preCell[nc][2]==1]

if len(connListExc_ThtoAll) != len(PMat_IntoAll.preCell):
    print("Please, consider also inhibitory projections from thalamus to cortex\n")


# Cortex -> Cortex (AlltoAll)
# List of excitatory connections
connListExc = [ [PMat_AlltoAll.preCell[nc][1],PMat_AlltoAll.preCell[nc][0]] for nc in range(len(PMat_AlltoAll.preCell)) if PMat_AlltoAll.preCell[nc][2]==1]
connListExc_gID = [ [PMat_AlltoAll.preCell[nc][1]-1,PMat_AlltoAll.preCell[nc][0]-1] for nc in range(len(PMat_AlltoAll.preCell)) if PMat_AlltoAll.preCell[nc][2]==1]
# weights (nA or equiv. mA/cm2 if thought of as distributed)
weightListExc = [ 0.001*PMat_AlltoAll.Am[nc] for nc in range(len(PMat_AlltoAll.preCell)) if PMat_AlltoAll.preCell[nc][2]==1]
CVListExc =     [ 0.001*PMat_AlltoAll.CV[nc] for nc in range(len(PMat_AlltoAll.preCell)) if PMat_AlltoAll.preCell[nc][2]==1]
#delays
delayListExc = [ PMat_AlltoAll.Delay[nc] for nc in range(len(PMat_AlltoAll.preCell)) if PMat_AlltoAll.preCell[nc][2]==1]
# short-term learning & failure
STDListExc = [ PMat_AlltoAll.STD[nc] for nc in range(len(PMat_AlltoAll.preCell)) if PMat_AlltoAll.preCell[nc][2]==1]
FailListExc = [ PMat_AlltoAll.Fail[nc] for nc in range(len(PMat_AlltoAll.preCell)) if PMat_AlltoAll.preCell[nc][2]==1]
PlasListExc = [ PMat_AlltoAll.Plas[nc] for nc in range(len(PMat_AlltoAll.preCell)) if PMat_AlltoAll.preCell[nc][2]==1]

# List of inhibitory connections
connListInh = [ [PMat_AlltoAll.preCell[nc][1],PMat_AlltoAll.preCell[nc][0]] for nc in range(len(PMat_AlltoAll.preCell)) if PMat_AlltoAll.preCell[nc][2]==-1]
connListInh_gID = [ [PMat_AlltoAll.preCell[nc][1]-1,PMat_AlltoAll.preCell[nc][0]-1] for nc in range(len(PMat_AlltoAll.preCell)) if PMat_AlltoAll.preCell[nc][2]==-1]
# weights (nA or equiv. mA/cm2 if thought of as distributed)
weightListInh = [ 0.001*PMat_AlltoAll.Am[nc] for nc in range(len(PMat_AlltoAll.preCell)) if PMat_AlltoAll.preCell[nc][2]==-1]
CVListInh =     [ 0.001*PMat_AlltoAll.CV[nc] for nc in range(len(PMat_AlltoAll.preCell)) if PMat_AlltoAll.preCell[nc][2]==-1]
#delays
delayListInh = [ PMat_AlltoAll.Delay[nc] for nc in range(len(PMat_AlltoAll.preCell)) if PMat_AlltoAll.preCell[nc][2]==-1]
# short-term learning & failure
STDListInh = [ PMat_AlltoAll.STD[nc] for nc in range(len(PMat_AlltoAll.preCell)) if PMat_AlltoAll.preCell[nc][2]==-1]
FailListInh = [ PMat_AlltoAll.Fail[nc] for nc in range(len(PMat_AlltoAll.preCell)) if PMat_AlltoAll.preCell[nc][2]==-1]
PlasListInh = [ PMat_AlltoAll.Plas[nc] for nc in range(len(PMat_AlltoAll.preCell)) if PMat_AlltoAll.preCell[nc][2]==-1]

## Connections by "connList" uses IDs relative to the "preConds" and "postConds"
netParams.connParams['Thalamus->All_exc'] = {
        'preConds': {'pop': 'artificial'},       # conditions of presyn cells
        'postConds': {'cellType': 'IzhiCell'},   # conditions of postsyn cells
        'connList': connListExc_ThtoAll_gID,    # list of conns
        'weight': weightListExc_ThtoAll,        # synaptic weight
        'delay': delayListExc_ThtoAll,          # transmission delay (ms)
        'synMech': 'exc'}

netParams.connParams['All->All_exc'] = {
        'preConds': {'cellType': 'IzhiCell'},   # conditions of presyn cells
        'postConds': {'cellType': 'IzhiCell'},  # conditions of postsyn cells
        'connList': connListExc_gID,            # list of conns
        'weight': weightListExc,                # synaptic weight
        'delay': delayListExc,                  # transmission delay (ms)
        'synMech': 'exc'}

netParams.connParams['All->All_inh'] = {
        'preConds': {'cellType': 'IzhiCell'},   # conditions of presyn cells
        'postConds': {'cellType': 'IzhiCell'},  # conditions of postsyn cells
        'connList': connListInh_gID,            # list of conns
        'weight': weightListInh,                # synaptic weight
        'delay': delayListInh,                  # transmission delay (ms)
        'synMech': 'inh'}


# # # Just an input instead of the thalamic inputs
# netParams.stimSourceParams['Input'] = {'type': 'NetStim', 'rate': 1000, 'noise': 0.3} # in mA/cm2
# for ntype in range(len(Pops)):
#     name = Pops[str(ntype+1)]['Label'] + '-' + Pops[str(ntype+1)]['Layer']
#     netParams.stimTargetParams['Input->'+name] = {'source': 'Input', 'sec':'soma', 'loc': 0.5, 'conds': {'pop':name}, 'weight': 0.005, 'synMech': 'exc'}


# Simulation options
simConfig = specs.SimConfig()       # object of class SimConfig to store simulation configuration
simConfig.duration = 6000           # Duration of the simulation, in ms
simConfig.dt = 0.1                  # Internal integration timestep to use
simConfig.verbose = False           # Show detailed messages
simConfig.recordTraces = {'V_soma1':{'sec':'soma','loc':0.5,'var':'v'}}  # Dict with traces to record
simConfig.recordStep = 0.1            # Step size in ms to save data (eg. V traces, LFP, etc)
simConfig.filename = 'sim1'         # Set file output name
simConfig.printRunTime = 0.1

simConfig.analysis['plotRaster'] = {'include': ['all'], 'saveFig': True}                  # Plot a raster
simConfig.analysis['plotTraces'] = {'include': [600,4818,4819,4820,9036], 'saveFig': True}  # Plot recorded traces for this list of cells

# Create the network with a generic synaptic mechanism
sim.create(netParams = netParams, simConfig = simConfig)

###############################################################################
## Modifying connections, first gIDs need to be consistent. It is assumed that populations are loaded by regions (first thalamus, then cortex. Or viceversa. But not different pops in cortex interleaved with thalamus) 
print("Modifying network ...\n")
Ncells_allin = len(sim.net.cells)
FirstThalamus_Gid = min([n for n in range(Ncells_allin) if 'VecStim' in list(sim.net.cells[n].tags.values())])
FirstCortex_Gid = min([n for n in range(Ncells_allin) if 'IzhiCell' in list(sim.net.cells[n].tags.values())])


# Modify connections (synaptic properties of each connection)
list_preGids = []
for npost in range(Ncells_allin):
    list_preGids.append( [sim.net.cells[npost].conns[nc]['preGid'] for nc in range(len(sim.net.cells[npost].conns))] )

# Thalamus to cortex
for nconn in range(len(connListExc_ThtoAll_gID)):
    Pre_Gid_rel = connListExc_ThtoAll_gID[nconn][0]
    Post_Gid_rel = connListExc_ThtoAll_gID[nconn][1]
    Pre_Gid = FirstThalamus_Gid + Pre_Gid_rel
    Post_Gid = FirstCortex_Gid + Post_Gid_rel
    try:
        conn_index = list_preGids[Post_Gid].index(Pre_Gid)
    except:
        print("Error in connections pre-post Ids")
        
    # Setting individual parameters
    sim.net.cells[Post_Gid].conns[conn_index]['hObj'].syn().tau_rise = PMat_IntoAll.Trise[Post_Gid_rel]
    sim.net.cells[Post_Gid].conns[conn_index]['hObj'].syn().tau_fall = PMat_IntoAll.Tfall[Post_Gid_rel]
    sim.net.cells[Post_Gid].conns[conn_index]['hObj'].syn().cn = PMat_IntoAll.CN[Post_Gid_rel]
    sim.net.cells[Post_Gid].conns[conn_index]['hObj'].syn().mean_amp = weightListExc_ThtoAll[nconn]
    sim.net.cells[Post_Gid].conns[conn_index]['hObj'].syn().cv = CVListExc_ThtoAll[nconn]
    sim.net.cells[Post_Gid].conns[conn_index]['hObj'].syn().std0 = STDListExc_ThtoAll[nconn]
    sim.net.cells[Post_Gid].conns[conn_index]['hObj'].syn().pf = FailListExc_ThtoAll[nconn]
    sim.net.cells[Post_Gid].conns[conn_index]['hObj'].syn().plas = PlasListExc_ThtoAll[nconn]
    sim.net.cells[Post_Gid].conns[conn_index]['hObj'].syn().tau_plas = tau_plas
    if nconn%100000==0: print("Thalamus->Cortex: ",nconn, "\n")

# Cortex to cortex
for nconn in range(len(connListExc_gID)):
    Pre_Gid_rel = connListExc_gID[nconn][0]
    Post_Gid_rel = connListExc_gID[nconn][1]
    Pre_Gid = FirstCortex_Gid + Pre_Gid_rel
    Post_Gid = FirstCortex_Gid + Post_Gid_rel
    try:
        conn_index = list_preGids[Post_Gid].index(Pre_Gid)
    except:
        print("Error in connections pre-post Ids")
        
    # Setting individual parameters
    Pop_pre = [nc for nc in range(len(N)) if (Pre_Gid_rel in N[nc])==True][0]
    sim.net.cells[Post_Gid].conns[conn_index]['hObj'].syn().tau_rise = PMat_AlltoAll.Trise[Post_Gid_rel][Pop_pre]
    sim.net.cells[Post_Gid].conns[conn_index]['hObj'].syn().tau_fall = PMat_AlltoAll.Tfall[Post_Gid_rel][Pop_pre]
    sim.net.cells[Post_Gid].conns[conn_index]['hObj'].syn().cn = PMat_AlltoAll.CN[Post_Gid_rel][Pop_pre]
    sim.net.cells[Post_Gid].conns[conn_index]['hObj'].syn().mean_amp = weightListExc[nconn]
    sim.net.cells[Post_Gid].conns[conn_index]['hObj'].syn().cv = CVListExc[nconn]
    sim.net.cells[Post_Gid].conns[conn_index]['hObj'].syn().std0 = STDListExc[nconn]
    sim.net.cells[Post_Gid].conns[conn_index]['hObj'].syn().pf = FailListExc[nconn]
    sim.net.cells[Post_Gid].conns[conn_index]['hObj'].syn().plas = PlasListExc[nconn]
    sim.net.cells[Post_Gid].conns[conn_index]['hObj'].syn().tau_plas = tau_plas
    if nconn%100000==0: print("Exc Cortex->Cortex: ",nconn, "\n")
    
for nconn in range(len(connListInh_gID)):
    Pre_Gid_rel = connListInh_gID[nconn][0]
    Post_Gid_rel = connListInh_gID[nconn][1]    
    Pre_Gid = FirstCortex_Gid + Pre_Gid_rel
    Post_Gid = FirstCortex_Gid + Post_Gid_rel
    try:
        conn_index = list_preGids[Post_Gid].index(Pre_Gid)
    except:
        print("Error in connections pre-post Ids")
        
    # Setting individual parameters
    Pop_pre = [nc for nc in range(len(N)) if (Pre_Gid_rel in N[nc])==True][0]
    sim.net.cells[Post_Gid].conns[conn_index]['hObj'].syn().tau_rise = PMat_AlltoAll.Trise[Post_Gid_rel][Pop_pre]
    sim.net.cells[Post_Gid].conns[conn_index]['hObj'].syn().tau_fall = PMat_AlltoAll.Tfall[Post_Gid_rel][Pop_pre]
    sim.net.cells[Post_Gid].conns[conn_index]['hObj'].syn().cn = PMat_AlltoAll.CN[Post_Gid_rel][Pop_pre]
    sim.net.cells[Post_Gid].conns[conn_index]['hObj'].syn().mean_amp = weightListInh[nconn]
    sim.net.cells[Post_Gid].conns[conn_index]['hObj'].syn().cv = CVListInh[nconn]
    sim.net.cells[Post_Gid].conns[conn_index]['hObj'].syn().std0 = STDListInh[nconn]
    sim.net.cells[Post_Gid].conns[conn_index]['hObj'].syn().pf = FailListInh[nconn]
    sim.net.cells[Post_Gid].conns[conn_index]['hObj'].syn().plas = PlasListInh[nconn]
    sim.net.cells[Post_Gid].conns[conn_index]['hObj'].syn().tau_plas = tau_plas
    if nconn%100000==0: print("Inh Cortex->Cortex: ",nconn, "\n")

print("Finish to set the network ...\n")
        
# simulate again
sim.cfg.filename = 'sim2'
sim.simulate()
sim.analyze()
