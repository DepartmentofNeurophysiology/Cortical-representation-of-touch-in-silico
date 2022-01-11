from netpyne import sim
import random
import numpy as np
from neuron import h

def run_multitrial(Conds, Input, inst_input, Nbarrels, Nthalamic, Barrels, simConfig):
    Ntrials = Conds['Ntrials']
    Nrepetitions = Conds['Nrepetitions']

    if sim.rank == 0: print("Finish to set the network ...\n")

    Nlocalcells = len(sim.net.cells)
    PostGids_nhost = [sim.net.cells[npost].gid for npost in range(Nlocalcells)]
    
    for ntrial in range(Ntrials):
        if sim.rank == 0: print("Trial: ",ntrial+1)

        ## Setting-up stimulus: Thalamic Inputs
        if Input['Option'] == 'Multitrial':
            binpsth =  inst_input[0].values.tolist()[1]-inst_input[0].values.tolist()[0]
            SpikeTimes_perBarrel = []
            for nb in range(Nbarrels):
                SpikeTimes_inBarrels = []
                for nn in range(Barrels['Barrel'+str(nb+1)]['Nthalamic']):
                    SpikeTms = []
                    for nt in range(len(inst_input)):
                        if Barrels['Barrel'+str(nb+1)]['Type'] == 'Principal':
                            pspike = inst_input[1].values.tolist()[nt] * binpsth / 1000.0
                        elif Barrels['Barrel'+str(nb+1)]['Type'] == 'Surround':
                            pspike = inst_input[2].values.tolist()[nt] * binpsth / 1000.0
                        else:
                            print("Barrel type not supported")
                        if random.random() < pspike:
                            SpikeTms.append(nt*binpsth)
                    SpikeTimes_inBarrels.append(SpikeTms)
                SpikeTimes_perBarrel.append(SpikeTimes_inBarrels)
            SpikeTimes = np.concatenate(SpikeTimes_perBarrel).ravel()

        elif Input['Option'] == 'Svoboda':
            # Check for number of trials
            Ntrials2 = len(inst_input['SpikeTrainStruct'][0].SpikeTimes[0])
            if Ntrials2 < Ntrials and sim.rank == 0:
                print('Spike trains from Svoboda dataset are not enough to complete Ntrials')
 
            # ntrial beyond the number of trials provided by SpikeTrainStruct are randomly repeated
            ntrial2 = ntrial
            if ntrial >= Ntrials2:
                ntrial2 = random.randint(0,Ntrials2-1)

            SpikeTimes_perBarrel = []
            for nbarrel in range(Nbarrels):
                SpikeTimes_perBarrel.append([ [inst_input['SpikeTrainStruct'][nbarrel].SpikeTimes[n,ntrial2]] if isinstance(inst_input['SpikeTrainStruct'][nbarrel].SpikeTimes[n,ntrial2],int) else inst_input['SpikeTrainStruct'][nbarrel].SpikeTimes[n,ntrial2].tolist() for n in range(len(inst_input['SpikeTrainStruct'][nbarrel].SpikeTimes)) ])
            SpikeTimes = np.concatenate(SpikeTimes_perBarrel).ravel()

        # Set spike times in VecStims
        for nc in range(Nthalamic):
            if nc in PostGids_nhost:
                index = PostGids_nhost.index(nc)
                sim.net.cells[index].hPointp.play(h.Vector().from_python(SpikeTimes[nc]))
    
        for nrep in range(Nrepetitions):
            if sim.rank == 0: print("   Repetition: ",nrep+1)
            simConfig.filename = 'sim'+str(ntrial+1)    # Set file output name
            if Nrepetitions>1:
                simConfig.filename = simConfig.filename + '_' + str(nrep+1)    # Set file output name
            sim.simulate()
            sim.analyze()

