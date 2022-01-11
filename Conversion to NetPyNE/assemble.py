# Data about barrels
def define_barrels(premodel):
    import numpy as np
    
    Nbarrels = premodel['Nbarrel']
    if Nbarrels != 1:
        BarrelStruct = premodel['barrelstruct']
    else:
        BarrelStruct = np.array([premodel['barrelstruct']])

    Barrels = {}
    for nb in range(Nbarrels):   ## one dimensional barrels (Nby = 1) - otherwise, it may change the access to BarrelStruct[nx][ny], for example, with proper loops on nx, ny
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

    return Nbarrels, Nthalamic, Barrels

# Data about cells and conns
def define_pops(model, settings):
    
    vr = settings['vr']
    v0 = settings['v0']
    u0 = settings['u0']
    dyn_thres = settings['dyn_thres']

    Neuron_Para = model['ConData'].Neuron_Para
    Neuron_Vt = model['ConData'].Neuron_Vt
    Cellinfo_In = model['ConData'].Cellinfo_In
    Cellinfo_All = model['ConData'].Cellinfo_All
    NtIn = model['ConData'].NtIn
    NtAll = model['ConData'].NtAll

    Ncells = Cellinfo_All.shape[0]
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
    Npops = len(Pops)

    N = [[n for n in range(Ncells) if Cellinfo_All[n,3]==(ntype+1)] for ntype in range(Npops)]

    for ntype in range(Npops):
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
        Pops[str(ntype+1)]['v0'] = v0
        Pops[str(ntype+1)]['u0'] = u0
        Pops[str(ntype+1)]['vt0'] = Neuron_Vt.avg[ntype]     # NOT USED IN THE ADAPTATIVE THRESHOLD MODEL
        Pops[str(ntype+1)]['a'] = Neuron_Para[N[ntype],0]
        Pops[str(ntype+1)]['b'] = Neuron_Para[N[ntype],1]
        Pops[str(ntype+1)]['c'] = Neuron_Para[N[ntype],2]
        Pops[str(ntype+1)]['d'] = Neuron_Para[N[ntype],3]
        if dyn_thres == 1:
            Pops[str(ntype+1)]['alpha'] = Neuron_Vt.alpha[ntype]
            Pops[str(ntype+1)]['Vi'] = Neuron_Vt.Vi[ntype]
            Pops[str(ntype+1)]['Vmin'] = Neuron_Vt.Vmin[ntype]
            Pops[str(ntype+1)]['Ka'] = Neuron_Vt.Ka[ntype]
            Pops[str(ntype+1)]['Ki'] = Neuron_Vt.Ki[ntype]
            Pops[str(ntype+1)]['tau'] = Neuron_Vt.tau[ntype]
    
        # Setting-up list, based on above properties, for instantiate individual neurons in NetPyNE
        Pops[str(ntype+1)]['list'] = [{'x': Pops[str(ntype+1)]['x'][n],
                                       'y': Pops[str(ntype+1)]['y'][n],
                                       'z': Pops[str(ntype+1)]['z'][n],
                                       'params':{'vr': Pops[str(ntype+1)]['vr'],
                                                 'v0': Pops[str(ntype+1)]['v0'],
                                                 'u0': Pops[str(ntype+1)]['u0'],
                                                 'vt0': Pops[str(ntype+1)]['vt0'],
                                                 'dyn_th': dyn_thres,
                                                 'a': Pops[str(ntype+1)]['a'][n],
                                                 'b': Pops[str(ntype+1)]['b'][n],
                                                 'c': Pops[str(ntype+1)]['c'][n],
                                                 'd': Pops[str(ntype+1)]['d'][n]}
                                       } for n in range(len(N[ntype]))]

        # adding properties for dynamical threshold evolution
        if dyn_thres == 1:
            for n in range(len(N[ntype])):
                Pops[str(ntype+1)]['list'][n]['params']['alpha'] = Pops[str(ntype+1)]['alpha']
                Pops[str(ntype+1)]['list'][n]['params']['Vi'] = Pops[str(ntype+1)]['Vi']
                Pops[str(ntype+1)]['list'][n]['params']['Vmin'] = Pops[str(ntype+1)]['Vmin']
                Pops[str(ntype+1)]['list'][n]['params']['Ka'] = Pops[str(ntype+1)]['Ka']
                Pops[str(ntype+1)]['list'][n]['params']['Ki'] = Pops[str(ntype+1)]['Ki']
                Pops[str(ntype+1)]['list'][n]['params']['tau'] = Pops[str(ntype+1)]['tau']


    return Ncells, Npops, N, Pops

def define_conns(model, settings):
    from random import random

    fr = settings['fr']
    
    PMat_IntoAll = model['ConData'].PMat_IntoAll
    PMat_AlltoAll = model['ConData'].PMat_AlltoAll

    # Assembling conns - Connections List: [Pre, Post], gIDs are relative to the regions (thalamus or cortex). This is consistent to the connParam rule set by "connList"
    # Thalamus -> Cortex (IntoAll)
    Exc_ThtoAll = {}
    Exc_ThtoAll['connList'] = [ [PMat_IntoAll.preCell[nc][1],PMat_IntoAll.preCell[nc][0]] for nc in range(len(PMat_IntoAll.preCell)) if PMat_IntoAll.preCell[nc][2]==1]
    Exc_ThtoAll['connList_gID'] = [ [PMat_IntoAll.preCell[nc][1]-1,PMat_IntoAll.preCell[nc][0]-1] for nc in range(len(PMat_IntoAll.preCell)) if PMat_IntoAll.preCell[nc][2]==1]
    # weights (nA or equiv. mA/cm2 if thought of as distributed)
    Exc_ThtoAll['weightList'] = [ 0.001*PMat_IntoAll.Am[nc] for nc in range(len(PMat_IntoAll.preCell)) if PMat_IntoAll.preCell[nc][2]==1]
    Exc_ThtoAll['CVList'] = [ 0.001*PMat_IntoAll.CV[nc] for nc in range(len(PMat_IntoAll.preCell)) if PMat_IntoAll.preCell[nc][2]==1]
    #delays
    Exc_ThtoAll['delayList'] = [ PMat_IntoAll.Delay[nc] for nc in range(len(PMat_IntoAll.preCell)) if PMat_IntoAll.preCell[nc][2]==1]
    # short-term learning & failure
    Exc_ThtoAll['STDList'] = [ PMat_IntoAll.STD[nc] for nc in range(len(PMat_IntoAll.preCell)) if PMat_IntoAll.preCell[nc][2]==1]
    Exc_ThtoAll['FailList'] = [ PMat_IntoAll.Fail[nc] for nc in range(len(PMat_IntoAll.preCell)) if PMat_IntoAll.preCell[nc][2]==1]
    Exc_ThtoAll['PlasList'] = [ PMat_IntoAll.Plas[nc] for nc in range(len(PMat_IntoAll.preCell)) if PMat_IntoAll.preCell[nc][2]==1]
    # temporal characteristics
    Exc_ThtoAll['Trise'] = [ PMat_IntoAll.Trise[PMat_IntoAll.preCell[nc][0]-1] for nc in range(len(PMat_IntoAll.preCell)) if PMat_IntoAll.preCell[nc][2]==1]
    Exc_ThtoAll['Tfall'] = [ PMat_IntoAll.Tfall[PMat_IntoAll.preCell[nc][0]-1] for nc in range(len(PMat_IntoAll.preCell)) if PMat_IntoAll.preCell[nc][2]==1]
    Exc_ThtoAll['CN'] = [ PMat_IntoAll.CN[PMat_IntoAll.preCell[nc][0]-1] for nc in range(len(PMat_IntoAll.preCell)) if PMat_IntoAll.preCell[nc][2]==1]
    

    if len(Exc_ThtoAll['connList']) != len(PMat_IntoAll.preCell):
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
    Exc_AlltoAll = {}
    Exc_AlltoAll['connList']     = [ [PMat_AlltoAll.preCell[idx_e[nc]][1],PMat_AlltoAll.preCell[idx_e[nc]][0]] for nc in range(len(idx_e)) ]
    Exc_AlltoAll['connList_gID'] = [ [PMat_AlltoAll.preCell[idx_e[nc]][1]-1,PMat_AlltoAll.preCell[idx_e[nc]][0]-1] for nc in range(len(idx_e)) ]
    Exc_AlltoAll['weightList']   = [ 0.001*PMat_AlltoAll.Am[idx_e[nc]] for nc in range(len(idx_e)) ]
    Exc_AlltoAll['CVList']       = [ 0.001*PMat_AlltoAll.CV[idx_e[nc]] for nc in range(len(idx_e)) ]
    Exc_AlltoAll['delayList']    = [ PMat_AlltoAll.Delay[idx_e[nc]] for nc in range(len(idx_e)) ]
    Exc_AlltoAll['STDList']      = [ PMat_AlltoAll.STD[idx_e[nc]] for nc in range(len(idx_e)) ]
    Exc_AlltoAll['FailList']     = [ PMat_AlltoAll.Fail[idx_e[nc]] for nc in range(len(idx_e)) ]
    Exc_AlltoAll['PlasList']     = [ PMat_AlltoAll.Plas[idx_e[nc]] for nc in range(len(idx_e)) ]
    Exc_AlltoAll['Trise']        = [ PMat_AlltoAll.Trise[PMat_AlltoAll.preCell[idx_e[nc]][0]-1][PMat_AlltoAll.preID[idx_e[nc]][1]-1] for nc in range(len(idx_e)) ]
    Exc_AlltoAll['Tfall']        = [ PMat_AlltoAll.Tfall[PMat_AlltoAll.preCell[idx_e[nc]][0]-1][PMat_AlltoAll.preID[idx_e[nc]][1]-1] for nc in range(len(idx_e)) ]
    Exc_AlltoAll['CN']           = [ PMat_AlltoAll.CN[PMat_AlltoAll.preCell[idx_e[nc]][0]-1][PMat_AlltoAll.preID[idx_e[nc]][1]-1] for nc in range(len(idx_e)) ]

    # List of inhibitory connections
    Inh_AlltoAll = {}
    Inh_AlltoAll['connList']     = [ [PMat_AlltoAll.preCell[idx_i[nc]][1],PMat_AlltoAll.preCell[idx_i[nc]][0]] for nc in range(len(idx_i)) ]
    Inh_AlltoAll['connList_gID'] = [ [PMat_AlltoAll.preCell[idx_i[nc]][1]-1,PMat_AlltoAll.preCell[idx_i[nc]][0]-1] for nc in range(len(idx_i)) ]
    Inh_AlltoAll['weightList']   = [ 0.001*PMat_AlltoAll.Am[idx_i[nc]] for nc in range(len(idx_i)) ]
    Inh_AlltoAll['CVList']       = [ 0.001*PMat_AlltoAll.CV[idx_i[nc]] for nc in range(len(idx_i)) ]
    Inh_AlltoAll['delayList']    = [ PMat_AlltoAll.Delay[idx_i[nc]] for nc in range(len(idx_i)) ]
    Inh_AlltoAll['STDList']      = [ PMat_AlltoAll.STD[idx_i[nc]] for nc in range(len(idx_i)) ]
    Inh_AlltoAll['FailList']     = [ PMat_AlltoAll.Fail[idx_i[nc]] for nc in range(len(idx_i)) ]
    Inh_AlltoAll['PlasList']     = [ PMat_AlltoAll.Plas[idx_i[nc]] for nc in range(len(idx_i)) ]
    Inh_AlltoAll['Trise']        = [ PMat_AlltoAll.Trise[PMat_AlltoAll.preCell[idx_i[nc]][0]-1][PMat_AlltoAll.preID[idx_i[nc]][1]-1] for nc in range(len(idx_i)) ]
    Inh_AlltoAll['Tfall']        = [ PMat_AlltoAll.Tfall[PMat_AlltoAll.preCell[idx_i[nc]][0]-1][PMat_AlltoAll.preID[idx_i[nc]][1]-1] for nc in range(len(idx_i)) ]
    Inh_AlltoAll['CN']           = [ PMat_AlltoAll.CN[PMat_AlltoAll.preCell[idx_i[nc]][0]-1][PMat_AlltoAll.preID[idx_i[nc]][1]-1] for nc in range(len(idx_i)) ]

    return Exc_ThtoAll, Exc_AlltoAll, Inh_AlltoAll