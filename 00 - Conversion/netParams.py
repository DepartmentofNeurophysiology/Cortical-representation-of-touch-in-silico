from netpyne import specs

def set_netParams(Nin, Pops, Exc_ThtoAll, Exc_AlltoAll, Inh_AlltoAll):
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
    netParams.popParams['artificial'] = {'cellModel': 'VecStim', 'numCells': Nin, 'spkTimes': [1]*Nin, 'xRange': [-0.01,0.01],'yRange': [0,0.01],'zRange': [-0.01,0.01]}

    # Populations in cortex
    for ntype in range(len(Pops)):
        name = Pops[str(ntype+1)]['Label'] + '-' + Pops[str(ntype+1)]['Layer']
        netParams.popParams[name] = {'cellType': 'IzhiCell', 'cellsList': Pops[str(ntype+1)]['list']}

    # Defining a generic synapse
    netParams.synMechParams['exc'] = {'mod': 'FluctExp2Syn', 'tau_rise': 1.0, 'tau_fall': 2.0, 'cn': 4.0, 'type': 1}
    netParams.synMechParams['inh'] = {'mod': 'FluctExp2Syn', 'tau_rise': 1.0, 'tau_fall': 2.0, 'cn': 4.0, 'type': -1}

    ## Connections by "connList" uses IDs relative to the "preConds" and "postConds"
    netParams.connParams['Thalamus->All_exc'] = {
            'preConds': {'pop': 'artificial'},          # conditions of presyn cells
            'postConds': {'cellType': 'IzhiCell'},      # conditions of postsyn cells
            'connList': Exc_ThtoAll['connList_gID'],    # list of conns
            'weight': Exc_ThtoAll['weightList'],        # synaptic weight
            'delay': Exc_ThtoAll['delayList'],          # transmission delay (ms)
            'synMech': 'exc'}

    netParams.connParams['All->All_exc'] = {
            'preConds': {'cellType': 'IzhiCell'},       # conditions of presyn cells
            'postConds': {'cellType': 'IzhiCell'},      # conditions of postsyn cells
            'connList': Exc_AlltoAll['connList_gID'],   # list of conns
            'weight': Exc_AlltoAll['weightList'],       # synaptic weight
            'delay': Exc_AlltoAll['delayList'],         # transmission delay (ms)
            'synMech': 'exc'}

    netParams.connParams['All->All_inh'] = {
            'preConds': {'cellType': 'IzhiCell'},       # conditions of presyn cells
            'postConds': {'cellType': 'IzhiCell'},      # conditions of postsyn cells
            'connList': Inh_AlltoAll['connList_gID'],   # list of conns
            'weight': Inh_AlltoAll['weightList'],       # synaptic weight
            'delay': Inh_AlltoAll['delayList'],         # transmission delay (ms)
            'synMech': 'inh'}

    return netParams