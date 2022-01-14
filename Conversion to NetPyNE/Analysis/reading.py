import os, json

def read_data(model):
    from scipy.io import  loadmat

    folder = model['Folder']
    data_premodel = model['PreModel']
    data_full = model['FullModel']
    
    reading_folder  = folder
    datamodel_folder = os.path.join(os.getcwd(),reading_folder)

    # Data about barrel structure (among others)
    filename = data_premodel
    reading_filename = os.path.join(datamodel_folder,filename)
    inst_premodel = loadmat(reading_filename, struct_as_record=False, squeeze_me=True)

    # Structured data about cells and connections
    filename = data_full
    reading_filename = os.path.join(datamodel_folder,filename)
    inst_model = loadmat(reading_filename, struct_as_record=False, squeeze_me=True)

    return inst_premodel, inst_model

def read_input(Input):
    
    option = Input['Option']
    folder = Input['Folder']
    filename = Input['Filename']

    reading_folder  = folder
    datamodel_folder = os.path.join(os.getcwd(),reading_folder)

    if option == 'Svoboda':
        # Psth and spike trains from Svoboda dataset
        from scipy.io import  loadmat
        inst_input = loadmat(os.path.join(datamodel_folder,filename), struct_as_record=False, squeeze_me=True)

    elif option == 'Multitrial':
        # Psth data from Aguilar's paper
        import pandas as pd
        inst_input = pd.read_csv(os.path.join(datamodel_folder,filename), sep=" ", header=None)
    else:
        print('Input option not available')
        
    return inst_input

def read_outputdata(Conds,vr):
    Ntrials = Conds['Ntrials']
    Nrepetitions = Conds['Nrepetitions']
    
    foldername_parent = os.getcwd()

    spkt = {}
    spkid = {}

    for nvr in range(len(vr)):
        foldername = str(vr[nvr])+'mV'
        os.chdir(foldername)

        spkt[str(vr[nvr])+'mV'] = []
        spkid[str(vr[nvr])+'mV'] = []
    
        for ntrial in range(Ntrials):
            for nrep in range(Nrepetitions):

                filename = 'sim'+str(ntrial+1)    # Set file input name
                if Nrepetitions>1:
                    filename = filename + '_' + str(nrep+1)    # Set file input name
            
                filename_ = open(filename+'.json')
                data = json.load(filename_)

                spkt[str(vr[nvr])+'mV'].append(data['simData']['spkt'])
                spkid[str(vr[nvr])+'mV'].append(data['simData']['spkid'])

        os.chdir(foldername_parent)

    return spkt, spkid