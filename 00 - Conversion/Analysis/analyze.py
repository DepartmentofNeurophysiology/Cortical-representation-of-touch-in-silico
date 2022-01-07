import numpy as np
import settings, reading, assemble, collect_stats


## Setting-up, initial conditions, options
Input, Model, ExperimConds, Settings = settings.setup()

## Reading input, model data
inst_input = reading.read_input(Input)
inst_premodel, inst_model = reading.read_data(Model)

## Conditioning data & assembling barrels and pops
Nbarrels, Nthalamic, Barrels = assemble.define_barrels(inst_premodel)
Ncells, Npops, N, Pops = assemble.define_pops(inst_model, Settings)

## Reading simulation output data
vr = [60, 70, 80]
spkt,spkid = reading.read_outputdata(ExperimConds,vr)

MeanStatistics_per_condition = collect_stats.meanstat(ExperimConds,vr,Barrels,Nthalamic,Pops,N,spkid)

# Detailed analysis
# Mean firing per stimulus
binSize = 0.5
bins_list = np.arange(-0.25, 5.25, binSize)

# Excitatory/inhibitory populations
popslabels = [ ['1', '2'],
                ['3', '4'],
                ['5'],
                ['6', '7', '8', '9', '10', '11', '12', '13', '14', '15'] ]

MeanStatistics_per_condition_GroupedPops = {}

# Principal barrel
nbarrel = 1
name = [ 'ExcL4_mainbarrel',
          'InhL4_mainbarrel',
          'ExcL23_mainbarrel',
          'InhL23_mainbarrel' ]
MeanStatistics_per_condition_GroupedPops['Principal'] = collect_stats.meanstat_grouped(MeanStatistics_per_condition,vr,nbarrel,popslabels,name,binSize,bins_list)

# Secondary barrels
nbarrel = [0,2]
name = [ 'ExcL4_secbarrels',
          'InhL4_secbarrels',
          'ExcL23_secbarrels',
          'InhL23_secbarrels' ]
MeanStatistics_per_condition_GroupedPops['Secondary'] = collect_stats.meanstat_grouped(MeanStatistics_per_condition,vr,nbarrel,popslabels,name,binSize,bins_list)

# Post-stimulus time histogram
timeRange = [-0.5, 40.5]
binSize = 0.5

# Excitatory/inhibitory populations
popslabels = [ ['1', '2'],
               ['3', '4'],
               ['5'],
               ['6', '7', '8', '9', '10', '11', '12', '13', '14', '15'] ]

Psth_per_condition_GroupedPops = {}

# Principal barrel
nbarrel = 1
name = [ 'ExcL4_mainbarrel',
         'InhL4_mainbarrel',
         'ExcL23_mainbarrel',
         'InhL23_mainbarrel' ]
Psth_per_condition_GroupedPops['Principal'] = collect_stats.psth_grouped(ExperimConds,vr,Barrels,Nthalamic,Pops,N,spkt,spkid,nbarrel,popslabels,name,timeRange,binSize)

nbarrel = [0,2]
name = [ 'ExcL4_secbarrels',
          'InhL4_secbarrels',
          'ExcL23_secbarrels',
          'InhL23_secbarrels' ]
Psth_per_condition_GroupedPops['Secondary'] = collect_stats.psth_grouped(ExperimConds,vr,Barrels,Nthalamic,Pops,N,spkt,spkid,nbarrel,popslabels,name,timeRange,binSize)


