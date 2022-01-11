import os

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
