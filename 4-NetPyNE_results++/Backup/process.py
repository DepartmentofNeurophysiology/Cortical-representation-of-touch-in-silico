import os
import sys
import json
import numpy as np
import pandas as pd
from random import random
from scipy.io import  loadmat
from netpyne import specs, sim
from neuron import h

###############################################################################
## Initial conditions and properties ##########################################
###############################################################################
vr = -60
u0 = 0
tau_plas = 120.0
fr = 1.00  ## Fraction of recurrent connections - Useful for setting-up sims ###

###############################################################################
## Input selection ############################################################
###############################################################################
Input = 'Multitrial'
#Input = 'Svoboda'
Nrepetitions = 1    # repetitions of an unique input condition
Ntrials = 10

###############################################################################
## Reading data ###############################################################
###############################################################################
reading_folder  = 'InstantiatedModel'
datamodel_folder = os.path.join(os.getcwd(),reading_folder)

# Data about barrel structure (among others)
filename = 'CMDMs_Test_sim_Svoboda.mat'
reading_filename = os.path.join(datamodel_folder,filename)
inst_premodel = loadmat(reading_filename, struct_as_record=False, squeeze_me=True)

# Structured data about cells and connections
filename = 'CMDMs_Test_sim_Svoboda_ConData.mat'
reading_filename = os.path.join(datamodel_folder,filename)
inst_model = loadmat(reading_filename, struct_as_record=False, squeeze_me=True)

# Psth and spike trains from Svoboda dataset
filename = 'Test_sim_Svoboda_an171923_19-Sep-2021_Thalamic_Spike_Trains.mat'
reading_filename = os.path.join(datamodel_folder,filename)
inst_input = loadmat(reading_filename, struct_as_record=False, squeeze_me=True)

# Psth data from Aguilar's paper
reading_folder  = 'ReadingData_Aguilar'
datamodel_folder = os.path.join(os.getcwd(),reading_folder)
pw_psth = pd.read_csv(os.path.join(datamodel_folder,'PrincipalWhisker_psth.dat'), sep=" ", header=None)
sw_psth = pd.read_csv(os.path.join(datamodel_folder,'SurroundWhisker_psth.dat'), sep=" ", header=None)

###############################################################################
## Conditioning data ##########################################################
###############################################################################
# Data about barrels
BarrelStruct = inst_premodel['barrelstruct']

# Data from cells and conns
PMat_IntoAll = inst_model['ConData'].PMat_IntoAll
PMat_AlltoAll = inst_model['ConData'].PMat_AlltoAll
Neuron_Para = inst_model['ConData'].Neuron_Para
Neuron_Vt = inst_model['ConData'].Neuron_Vt
Cellinfo_In = inst_model['ConData'].Cellinfo_In
Cellinfo_All = inst_model['ConData'].Cellinfo_All
NtIn = inst_model['ConData'].NtIn
NtAll = inst_model['ConData'].NtAll

# Data from spike train inputs. Curated from Svoboda dataset
SpikeTrainStruct = inst_input['SpikeTrainStruct']


###############################################################################
## Assembling barrels #########################################################
###############################################################################
Nbarrels = len(BarrelStruct)
Barrels = {}
for nb in range(Nbarrels):
    Barrels['Barrel'+str(nb+1)] = {
        'Nthalamic': BarrelStruct[nb].Nthalamic,
        'Xpos': BarrelStruct[nb].xpos,
        'Ypos': BarrelStruct[nb].ypos
        }
    if BarrelStruct[nb].mainbarrel == 1:
        Barrels['Barrel'+str(nb+1)]['Type'] = 'Principal'
    elif BarrelStruct[nb].mainbarrel == 2:
        Barrels['Barrel'+str(nb+1)]['Type'] = 'Surround'
    else:
        print('Barrel category not recognized')
        
Nthalamic = sum(Barrels['Barrel'+str(nb+1)]['Nthalamic'] for nb in range(Nbarrels))

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

    # reading barrel
    Pops[str(ntype+1)]['Barrel'] = Cellinfo_All[N[ntype],4]

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

netParams.defaultThreshold = 0.0

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
netParams.popParams['artificial'] = {'cellModel': 'VecStim', 'numCells': Nthalamic, 'spkTimes': [1]*Nthalamic, 'xRange': [-0.01,0.01],'yRange': [0,0.01],'zRange': [-0.01,0.01]}

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
idx_e = []
idx_i = []
for nc in range(len(PMat_AlltoAll.preCell)):
    if random() < fr:
        if PMat_AlltoAll.preCell[nc][2]==1:
            idx_e.append(nc)
        elif PMat_AlltoAll.preCell[nc][2]==-1:
            idx_i.append(nc)
        else:
            print("Connection without specification of exc/inh")

# List of excitatory connections
connListExc     = [ [PMat_AlltoAll.preCell[idx_e[nc]][1],PMat_AlltoAll.preCell[idx_e[nc]][0]] for nc in range(len(idx_e)) ]
connListExc_gID = [ [PMat_AlltoAll.preCell[idx_e[nc]][1]-1,PMat_AlltoAll.preCell[idx_e[nc]][0]-1] for nc in range(len(idx_e)) ]
weightListExc   = [ 0.001*PMat_AlltoAll.Am[idx_e[nc]] for nc in range(len(idx_e)) ]
CVListExc       = [ 0.001*PMat_AlltoAll.CV[idx_e[nc]] for nc in range(len(idx_e)) ]
delayListExc    = [ PMat_AlltoAll.Delay[idx_e[nc]] for nc in range(len(idx_e)) ]
STDListExc      = [ PMat_AlltoAll.STD[idx_e[nc]] for nc in range(len(idx_e)) ]
FailListExc     = [ PMat_AlltoAll.Fail[idx_e[nc]] for nc in range(len(idx_e)) ]
PlasListExc     = [ PMat_AlltoAll.Plas[idx_e[nc]] for nc in range(len(idx_e)) ]

# List of inhibitory connections
connListInh     = [ [PMat_AlltoAll.preCell[idx_i[nc]][1],PMat_AlltoAll.preCell[idx_i[nc]][0]] for nc in range(len(idx_i)) ]
connListInh_gID = [ [PMat_AlltoAll.preCell[idx_i[nc]][1]-1,PMat_AlltoAll.preCell[idx_i[nc]][0]-1] for nc in range(len(idx_i)) ]
weightListInh   = [ 0.001*PMat_AlltoAll.Am[idx_i[nc]] for nc in range(len(idx_i)) ]
CVListInh       = [ 0.001*PMat_AlltoAll.CV[idx_i[nc]] for nc in range(len(idx_i)) ]
delayListInh    = [ PMat_AlltoAll.Delay[idx_i[nc]] for nc in range(len(idx_i)) ]
STDListInh      = [ PMat_AlltoAll.STD[idx_i[nc]] for nc in range(len(idx_i)) ]
FailListInh     = [ PMat_AlltoAll.Fail[idx_i[nc]] for nc in range(len(idx_i)) ]
PlasListInh     = [ PMat_AlltoAll.Plas[idx_i[nc]] for nc in range(len(idx_i)) ]


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


# Simulation options
simConfig = specs.SimConfig()       # object of class SimConfig to store simulation configuration
simConfig.duration = 40             # Duration of the simulation, in ms
simConfig.dt = 0.025                # Internal integration timestep to use
simConfig.verboseLevel = 1          # Show detailed messages
simConfig.recordTraces = {'v_soma':{'sec':'soma','loc':0.5,'var':'v'},
                          'u_soma':{'sec':'soma', 'pointp':'Izhi','var':'u'},
                          'vt_soma':{'sec':'soma', 'pointp':'Izhi','var':'vt'},
                          'synapse_exc_std':{'sec':'soma','loc':0.5,'synMech':'exc','var':'std'},
                          'synapse_exc_s':{'sec':'soma','loc':0.5,'synMech':'exc','var':'s'},
                          'synapse_inh_std':{'sec':'soma','loc':0.5,'synMech':'inh','var':'std'},
                          'synapse_inh_s':{'sec':'soma','loc':0.5,'synMech':'inh','var':'s'}
                          }  # Dict with traces to record


simConfig.recordStep = 0.025        # Step size in ms to save data (eg. V traces, LFP, etc)
simConfig.gatherOnlySimData = False
simConfig.saveDataInclude = ['simData']
simConfig.saveJson = True

simConfig.analysis['plotRaster'] = {'include': ['all'], 'saveFig': True}                  # Plot a raster
simConfig.analysis['plotTraces'] = {'include': [1260, 3321, 4863, 5144, 8631, 12138, 12303, 12392, 12570, 12762, 12909, 13050, 13197], 'saveFig': True}  # Plot recorded traces for this list of cells

# Create the network with a generic synaptic mechanism
sim.create(netParams = netParams, simConfig = simConfig, output=True)

###############################################################################
# Modifying connections, first gIDs need to be consistent. It is assumed that populations are loaded by regions (first thalamus, then cortex. Or viceversa. But not different pops in cortex interleaved with thalamus) 
if sim.rank == 0: print("Modifying network ...\n")

Ncells_allin = len(sim.net.cells)
FirstThalamus_Gid = 0
FirstCortex_Gid = 600

# Modify connections (synaptic properties of each connection)
list_preGids = []
for npost in range(Ncells_allin):
    list_preGids.append( [sim.net.cells[npost].conns[nc]['preGid'] for nc in range(len(sim.net.cells[npost].conns))] )

PostGids_nhost = [sim.net.cells[npost].gid for npost in range(Ncells_allin)]

# with open("kk"+str(sim.rank)+".txt","w") as f:
#     for npost in range(Ncells_allin):
#         print(sim.net.cells[npost].gid, list_preGids[npost], file=f)

# with open("kkhuate"+str(sim.rank)+".txt","w") as f:
#     print(PostGids_nhost,file=f)
    
# Thalamus to cortex
for nconn in range(len(connListExc_ThtoAll_gID)):

    if nconn%100000 == 0 and sim.rank == 0: print("Thalamus->Cortex: ",nconn, "\n")

    Pre_Gid_rel = connListExc_ThtoAll_gID[nconn][0]
    Post_Gid_rel = connListExc_ThtoAll_gID[nconn][1]
    Pre_Gid = FirstThalamus_Gid + Pre_Gid_rel
    Post_Gid = FirstCortex_Gid + Post_Gid_rel
    if Post_Gid in PostGids_nhost:
        indexPost_Gid = PostGids_nhost.index(Post_Gid)
        try:
            conn_index = list_preGids[indexPost_Gid].index(Pre_Gid)
        except:
            print("Error in connections pre-post Ids")
        
        # Setting individual parameters
        sim.net.cells[indexPost_Gid].conns[conn_index]['hObj'].syn().tau_rise = PMat_IntoAll.Trise[Post_Gid_rel]
        sim.net.cells[indexPost_Gid].conns[conn_index]['hObj'].syn().tau_fall = PMat_IntoAll.Tfall[Post_Gid_rel]
        sim.net.cells[indexPost_Gid].conns[conn_index]['hObj'].syn().cn = PMat_IntoAll.CN[Post_Gid_rel]
        sim.net.cells[indexPost_Gid].conns[conn_index]['hObj'].syn().mean_amp = weightListExc_ThtoAll[nconn]
        sim.net.cells[indexPost_Gid].conns[conn_index]['hObj'].syn().cv = CVListExc_ThtoAll[nconn]
        sim.net.cells[indexPost_Gid].conns[conn_index]['hObj'].syn().std0 = STDListExc_ThtoAll[nconn]
        sim.net.cells[indexPost_Gid].conns[conn_index]['hObj'].syn().pf = FailListExc_ThtoAll[nconn]
        sim.net.cells[indexPost_Gid].conns[conn_index]['hObj'].syn().plas = PlasListExc_ThtoAll[nconn]
        sim.net.cells[indexPost_Gid].conns[conn_index]['hObj'].syn().tau_plas = tau_plas

# Cortex to cortex
# Excitatory
for nconn in range(len(connListExc_gID)):
    if nconn%100000 == 0 and sim.rank == 0: print("Exc Cortex->Cortex: ",nconn, "\n")

    Pre_Gid_rel = connListExc_gID[nconn][0]
    Post_Gid_rel = connListExc_gID[nconn][1]
    Pre_Gid = FirstCortex_Gid + Pre_Gid_rel
    Post_Gid = FirstCortex_Gid + Post_Gid_rel
    if Post_Gid in PostGids_nhost:
        indexPost_Gid = PostGids_nhost.index(Post_Gid)
        try:
            conn_index = list_preGids[indexPost_Gid].index(Pre_Gid)
        except:
            print("Error in connections pre-post Ids")
        
    # Setting individual parameters
    Pop_pre = [nc for nc in range(len(N)) if (Pre_Gid_rel in N[nc])==True][0]
    sim.net.cells[indexPost_Gid].conns[conn_index]['hObj'].syn().tau_rise = PMat_AlltoAll.Trise[Post_Gid_rel][Pop_pre]
    sim.net.cells[indexPost_Gid].conns[conn_index]['hObj'].syn().tau_fall = PMat_AlltoAll.Tfall[Post_Gid_rel][Pop_pre]
    sim.net.cells[indexPost_Gid].conns[conn_index]['hObj'].syn().cn = PMat_AlltoAll.CN[Post_Gid_rel][Pop_pre]
    sim.net.cells[indexPost_Gid].conns[conn_index]['hObj'].syn().mean_amp = weightListExc[nconn]
    sim.net.cells[indexPost_Gid].conns[conn_index]['hObj'].syn().cv = CVListExc[nconn]
    sim.net.cells[indexPost_Gid].conns[conn_index]['hObj'].syn().std0 = STDListExc[nconn]
    sim.net.cells[indexPost_Gid].conns[conn_index]['hObj'].syn().pf = FailListExc[nconn]
    sim.net.cells[indexPost_Gid].conns[conn_index]['hObj'].syn().plas = PlasListExc[nconn]
    sim.net.cells[indexPost_Gid].conns[conn_index]['hObj'].syn().tau_plas = tau_plas
    
    
# Inhibitory
for nconn in range(len(connListInh_gID)):
    if nconn%100000 == 0 and sim.rank == 0: print("Inh Cortex->Cortex: ",nconn, "\n")

    Pre_Gid_rel = connListInh_gID[nconn][0]
    Post_Gid_rel = connListInh_gID[nconn][1]
    Pre_Gid = FirstCortex_Gid + Pre_Gid_rel
    Post_Gid = FirstCortex_Gid + Post_Gid_rel
    if Post_Gid in PostGids_nhost:
        indexPost_Gid = PostGids_nhost.index(Post_Gid)
        try:
            conn_index = list_preGids[indexPost_Gid].index(Pre_Gid)
        except:
            print("Error in connections pre-post Ids")

    # Setting individual parameters
    Pop_pre = [nc for nc in range(len(N)) if (Pre_Gid_rel in N[nc])==True][0]
    sim.net.cells[indexPost_Gid].conns[conn_index]['hObj'].syn().tau_rise = PMat_AlltoAll.Trise[Post_Gid_rel][Pop_pre]
    sim.net.cells[indexPost_Gid].conns[conn_index]['hObj'].syn().tau_fall = PMat_AlltoAll.Tfall[Post_Gid_rel][Pop_pre]
    sim.net.cells[indexPost_Gid].conns[conn_index]['hObj'].syn().cn = PMat_AlltoAll.CN[Post_Gid_rel][Pop_pre]
    sim.net.cells[indexPost_Gid].conns[conn_index]['hObj'].syn().mean_amp = weightListInh[nconn]
    sim.net.cells[indexPost_Gid].conns[conn_index]['hObj'].syn().cv = CVListInh[nconn]
    sim.net.cells[indexPost_Gid].conns[conn_index]['hObj'].syn().std0 = STDListInh[nconn]
    sim.net.cells[indexPost_Gid].conns[conn_index]['hObj'].syn().pf = FailListInh[nconn]
    sim.net.cells[indexPost_Gid].conns[conn_index]['hObj'].syn().plas = PlasListInh[nconn]
    sim.net.cells[indexPost_Gid].conns[conn_index]['hObj'].syn().tau_plas = tau_plas

if sim.rank == 0: print("Finish to set the network ...\n")
        
for ntrial in range(Ntrials):
    if sim.rank == 0: print("Trial: ",ntrial+1)
    ###########################################################################
    ## Setting-up stimulus: Thalamic Inputs ###################################
    ###########################################################################
    if Input == 'Multitrial':
        binpsth = 1
        SpikeTimes_perBarrel = []
        for nb in range(Nbarrels):
            SpikeTimes_inBarrels = []
            for nn in range(Barrels['Barrel'+str(nb+1)]['Nthalamic']):
                SpikeTms = []
                for nt in range(len(pw_psth)):
                    if Barrels['Barrel'+str(nb+1)]['Type'] == 'Principal':
                        pspike = pw_psth[1][nt] * binpsth / 1000.0
                    elif Barrels['Barrel'+str(nb+1)]['Type'] == 'Surround':
                        pspike = sw_psth[1][nt] * binpsth / 1000.0
                    else:
                        print("Barrel type not supported")
                    if random() < pspike:
                        SpikeTms.append(nt*binpsth)
                SpikeTimes_inBarrels.append(SpikeTms)
            SpikeTimes_perBarrel.append(SpikeTimes_inBarrels)
        SpikeTimes = np.concatenate(SpikeTimes_perBarrel).ravel()

    elif Input == 'Svoboda':
        SpikeTimes_perBarrel = []
        for nbarrel in range(Nbarrels):
            SpikeTimes_perBarrel.append([ [SpikeTrainStruct[nbarrel].SpikeTimes[n,ntrial]] if isinstance(SpikeTrainStruct[nbarrel].SpikeTimes[n,ntrial],int) else SpikeTrainStruct[nbarrel].SpikeTimes[n,ntrial].tolist() for n in range(len(SpikeTrainStruct[nbarrel].SpikeTimes)) ])
        SpikeTimes = np.concatenate(SpikeTimes_perBarrel).ravel()

    # Set spike times in VecStims
    for nc in range(Nthalamic):
        if nc in PostGids_nhost:
            index = PostGids_nhost.index(nc)
            sim.net.cells[index].hPointp.play(h.Vector().from_python(SpikeTimes[nc]))
    
    simConfig.filename = 'sim'+str(ntrial+1)         # Set file output name
    sim.simulate()
    sim.analyze()