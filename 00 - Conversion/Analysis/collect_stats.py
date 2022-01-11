import numpy as np

def meanstat(Conds, vr, Barrels, Nthalamic, Pops, N, spkid):

    Ntrials = Conds['Ntrials']
    Nrepetitions = Conds['Nrepetitions']
    Nbarrels = len(Barrels)
    
    MeanStatistics_per_condition = []
    for nvr in range(len(vr)):
        Statistics_per_barrel = []
        for nb in range(Nbarrels):
            Statistics = {}
            # Thalamic population
            Statistics['Thalamic'] = []
            for nn in range(Barrels['Barrel'+str(nb+1)]['Nthalamic']):
                spike_count = []
                for ntrial in range(Ntrials):
                    spkid_trial = spkid[str(vr[nvr])+'mV'][ntrial*Nrepetitions]   # only first repetition
                    spike_count.append(sum([1 for ns in range(len(spkid_trial)) if spkid_trial[ns]==(nb*Nthalamic/Nbarrels)+nn]))
                Statistics['Thalamic'].append(np.mean(spike_count))
                 
            # Cortical populations
            for ntype in range(len(Pops)):
                Statistics[str(ntype+1)] = []
                for nn in range(len(N[ntype])):
                    gid = Nthalamic + N[ntype][nn]
                    if Pops[str(ntype+1)]['Barrel'][nn]==(nb+1):
                        spike_count = []
                        for ntrial in range(Ntrials):
                            spkid_trial = spkid[str(vr[nvr])+'mV'][ntrial*Nrepetitions]   # only first repetition
                            spike_count.append(sum([1 for ns in range(len(spkid_trial)) if spkid_trial[ns]==gid]))
                        Statistics[str(ntype+1)].append(np.mean(spike_count))

            Statistics_per_barrel.append(Statistics)
        MeanStatistics_per_condition.append(Statistics_per_barrel)

    return MeanStatistics_per_condition

def meanstat_grouped(Statistics,vr,nbarrel,popslabels,name,binSize,bins_list):
    import matplotlib as plt

    MeanStatistics = {}
    histoT = bins_list[:-1]+binSize/2
    for nl in range(len(popslabels)):
        MeanStatistics_per_condition = []
        for nvr in range(len(vr)):
            if isinstance(nbarrel,int):
                sample = [val for elem in [Statistics[nvr][nbarrel][label] for label in popslabels[nl]] for val in elem]
            elif isinstance(nbarrel,list):
                sample = []
                for nb in range(len(nbarrel)):
                    sample.append([val for elem in [Statistics[nvr][nbarrel[nb]][label] for label in popslabels[nl]] for val in elem])
                sample = [val for elem in sample for val in elem]
            else:
                print('Data format not supported')
        
            histo = np.histogram(sample,bins_list)
            MeanStatistics_per_condition.append(histo[0]/(binSize*len(sample)))

            if nvr == 0: color='black'
            elif nvr == 1: color='blue'
            else: color = 'red'

            plt.pyplot.plot(histoT,MeanStatistics_per_condition[nvr],color=color)

        plt.pyplot.legend(loc='upper right')
        plt.pyplot.suptitle(name[nl])
        plt.pyplot.xlabel('Mean Spike Counts per Stimulus')
        plt.pyplot.ylabel('Probability density')
        plt.pyplot.savefig(name[nl]+'.png',bbox_inches='tight')
        plt.pyplot.show()
        
        MeanStatistics[name[nl]] = MeanStatistics_per_condition

    return MeanStatistics

def psth_grouped(Conds, vr, Barrels, Nthalamic, Pops, N, spkt, spkid, nbarrel, popslabels, name, timeRange, binSize):
    from builtins import zip
    import matplotlib as plt

    Ntrials = Conds['Ntrials']
    Nrepetitions = Conds['Nrepetitions']
    Nbarrels = len(Barrels)
    
    bins = np.arange(timeRange[0], timeRange[1], binSize)
    histo = np.histogram([],bins)
    histoT = histo[1][:-1]+binSize/2
    
    Cells_ids = {}
    if isinstance(nbarrel,int):
        # Thalamic cells
        Cells_ids['Thalamic'] = np.arange(nbarrel * Nthalamic/Nbarrels,(nbarrel+1) * Nthalamic/Nbarrels).tolist()
        # Cortical cells
        for ntype in range(len(Pops)):
            Cells_ids[str(ntype+1)] = [(Nthalamic + Pops[str(ntype+1)]['ids'][n]) for n in range(len(N[ntype])) if Pops[str(ntype+1)]['Barrel'][n] == (nbarrel+1)]
        # Grouped pops
        for nl in range(len(popslabels)):
            Cells_ids[name[nl]] = []
            for n in range(len(popslabels[nl])):
                Cells_ids[name[nl]].append(Cells_ids[popslabels[nl][n]])
            Cells_ids[name[nl]] = [val for elem in Cells_ids[name[nl]] for val in elem]

    elif isinstance(nbarrel,list):
        # Thalamic cells
        Cells_ids['Thalamic'] = [val for elem in [np.arange(nb * Nthalamic/Nbarrels,(nb+1) * Nthalamic/Nbarrels).tolist() for nb in nbarrel] for val in elem]
        # Cortical cells
        for ntype in range(len(Pops)):
            Cells_ids[str(ntype+1)] = [ val for elem in [ [(Nthalamic + Pops[str(ntype+1)]['ids'][n]) for n in range(len(N[ntype])) if Pops[str(ntype+1)]['Barrel'][n] == (nb+1)] for nb in nbarrel] for val in elem]
        # Grouped pops
        for nl in range(len(popslabels)):
            Cells_ids[name[nl]] = []
            for n in range(len(popslabels[nl])):
                Cells_ids[name[nl]].append(Cells_ids[popslabels[nl][n]])
            Cells_ids[name[nl]] = [val for elem in Cells_ids[name[nl]] for val in elem]

    else:
        print('Data format not supported')
    
    # Creating psth
    spkt_multitrial = {}
    spkid_multitrial = {}
    for nvr in range(len(vr)):
        spkt_multitrial[str(vr[nvr])+'mV'] = []
        spkid_multitrial[str(vr[nvr])+'mV'] = []
        for ntrial in range(Ntrials):
            spkt_multitrial[str(vr[nvr])+'mV'].append(spkt[str(vr[nvr])+'mV'][ntrial*Nrepetitions])   # only first repetition
            spkid_multitrial[str(vr[nvr])+'mV'].append(spkid[str(vr[nvr])+'mV'][ntrial*Nrepetitions])   # only first repetition
        spkt_multitrial[str(vr[nvr])+'mV'] = [ val for elem in spkt_multitrial[str(vr[nvr])+'mV'] for val in elem ]
        spkid_multitrial[str(vr[nvr])+'mV'] = [ val for elem in spkid_multitrial[str(vr[nvr])+'mV'] for val in elem ]

    histoCount = {}
    for nvr in range(len(vr)):
        
        spkt_ = spkt_multitrial[str(vr[nvr])+'mV']
        spkid_ = spkid_multitrial[str(vr[nvr])+'mV']
        histoCount[str(vr[nvr])+'mV'] = {}
        
        # Thalamic input
        spkinds,spkts = list(zip(*[(spkgid,spkt) for spkgid,spkt in zip(spkid_,spkt_) if spkgid in Cells_ids['Thalamic']]))
        histoCount[str(vr[nvr])+'mV']['thalamic'], bin_edges = np.histogram(spkts, bins)
        histoCount[str(vr[nvr])+'mV']['thalamic'] = histoCount[str(vr[nvr])+'mV']['thalamic']*(1000.0 / binSize) / (len(Cells_ids['Thalamic']*Ntrials))

        # Cortical populations
        for ntype in range(len(Pops)):
            if len(N[ntype]) > 0:
                try:
                    spkinds,spkts = list(zip(*[(spkgid,spkt) for spkgid,spkt in zip(spkid_,spkt_) if spkgid in Cells_ids[str(ntype+1)]]))
                except:
                    spkinds,spkts = [],[]
                histoCount[str(vr[nvr])+'mV'][str(ntype+1)], bin_edges = np.histogram(spkts, bins)
                histoCount[str(vr[nvr])+'mV'][str(ntype+1)] = histoCount[str(vr[nvr])+'mV'][str(ntype+1)]*(1000.0 / binSize) / (len(Cells_ids[str(ntype+1)]*Ntrials))

        # Grouped pops
        for nl in range(len(popslabels)):
            try:
                spkinds,spkts = list(zip(*[(spkgid,spkt) for spkgid,spkt in zip(spkid_,spkt_) if spkgid in Cells_ids[name[nl]]]))
            except:
                spkinds,spkts = [],[]
            histoCount[str(vr[nvr])+'mV'][name[nl]], bin_edges = np.histogram(spkts, bins)
            histoCount[str(vr[nvr])+'mV'][name[nl]] = histoCount[str(vr[nvr])+'mV'][name[nl]]*(1000.0 / binSize) / (len(Cells_ids[name[nl]]*Ntrials))

    for nl in range(len(popslabels)):
        for nvr in range(len(vr)):
            if nvr == 0: color='black'
            elif nvr == 1: color='blue'
            else: color = 'red'

            plt.pyplot.plot(histoT,histoCount[str(vr[nvr])+'mV'][name[nl]],color = color)
        plt.pyplot.title(name[nl])
        plt.pyplot.xticks([0,10,20,30,40])
        plt.pyplot.yticks([0,300,600])
        plt.pyplot.xlim([0,40])
        plt.pyplot.ylim([0,700])
        plt.pyplot.ylabel('Firing rate (Hz)')

        plt.pyplot.savefig(name[nl]+'_psth.png',bbox_inches='tight')
        plt.pyplot.show()

    return histoCount
