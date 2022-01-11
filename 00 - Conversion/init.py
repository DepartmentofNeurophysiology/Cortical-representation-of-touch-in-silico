import settings, reading, assemble
import netParams, cfg, netModify, multitrial
from netpyne import sim

## Setting-up, initial conditions, options
Input, Model, ExperimConds, Settings = settings.setup()

## Reading input, model data
inst_input = reading.read_input(Input)
inst_premodel, inst_model = reading.read_data(Model)

## Conditioning data & assembling barrels, pops, conns
Nbarrels, Nthalamic, Barrels = assemble.define_barrels(inst_premodel)
Ncells, Npops, N, Pops = assemble.define_pops(inst_model, Settings)
Exc_ThtoAll, Exc_AlltoAll, Inh_AlltoAll  = assemble.define_conns(inst_model, Settings)

## NetPyNE specifications
netParams = netParams.set_netParams(Nthalamic, Pops, Exc_ThtoAll, Exc_AlltoAll, Inh_AlltoAll)
simConfig = cfg.set_cfg(Input,Nbarrels)

# Create the network with a generic synaptic mechanism
sim.create(netParams = netParams, simConfig = simConfig, output=True)

# Modifying the instantiated network with particularized synaptic information
netModify.update_net(Nthalamic,Ncells,Exc_ThtoAll, Exc_AlltoAll, Inh_AlltoAll, Settings)

# Running multitrial experiment - plotting raster and traces (set at simConfig)- saving
multitrial.run_multitrial(ExperimConds, Input, inst_input, Nbarrels, Nthalamic, Barrels, simConfig)