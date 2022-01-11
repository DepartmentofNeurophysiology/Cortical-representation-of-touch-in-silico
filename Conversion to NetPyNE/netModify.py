from netpyne import sim

def update_net(Nth, Ncx, Exc_ThtoAll, Exc_AlltoAll, Inh_AlltoAll, settings):
    if sim.rank == 0: print("Modifying network ...\n")

    tau_plas = settings['tau_plas']
    FirstThalamus_Gid, FirstCortex_Gid = 0, Nth

    Nlocalcells = len(sim.net.cells)

    # Modify connections (synaptic properties of each connection)
    list_preGids = []
    for npost in range(Nlocalcells):
        list_preGids.append( [sim.net.cells[npost].conns[nc]['preGid'] for nc in range(len(sim.net.cells[npost].conns))] )

    PostGids_nhost = [sim.net.cells[npost].gid for npost in range(Nlocalcells)]

    with open("Conns"+str(sim.rank)+".txt","w") as f:
        for npost in range(Nlocalcells):
            print(sim.net.cells[npost].gid, list_preGids[npost], file=f)

    with open("ListPost"+str(sim.rank)+".txt","w") as f:
        print(PostGids_nhost,file=f)
        
    
    # Thalamus to cortex
    for nconn in range(len(Exc_ThtoAll['connList_gID'])):

        if nconn%100000 == 0 and sim.rank == 0: print("Thalamus->Cortex: ",nconn, "\n")

        Pre_Gid_rel = Exc_ThtoAll['connList_gID'][nconn][0]
        Post_Gid_rel = Exc_ThtoAll['connList_gID'][nconn][1]
        Pre_Gid = FirstThalamus_Gid + Pre_Gid_rel
        Post_Gid = FirstCortex_Gid + Post_Gid_rel
        if Post_Gid in PostGids_nhost:
            indexPost_Gid = PostGids_nhost.index(Post_Gid)
            try:
                conn_index = list_preGids[indexPost_Gid].index(Pre_Gid)
            except:
                print("Error in connections pre-post Ids")
        
            # Setting individual parameters
            sim.net.cells[indexPost_Gid].conns[conn_index]['hObj'].syn().tau_rise = Exc_ThtoAll['Trise'][nconn]
            sim.net.cells[indexPost_Gid].conns[conn_index]['hObj'].syn().tau_fall = Exc_ThtoAll['Tfall'][nconn]
            sim.net.cells[indexPost_Gid].conns[conn_index]['hObj'].syn().cn = Exc_ThtoAll['CN'][nconn]
            sim.net.cells[indexPost_Gid].conns[conn_index]['hObj'].syn().mean_amp = Exc_ThtoAll['weightList'][nconn]
            sim.net.cells[indexPost_Gid].conns[conn_index]['hObj'].syn().cv = Exc_ThtoAll['CVList'][nconn]
            sim.net.cells[indexPost_Gid].conns[conn_index]['hObj'].syn().std0 = Exc_ThtoAll['STDList'][nconn]
            sim.net.cells[indexPost_Gid].conns[conn_index]['hObj'].syn().pf = Exc_ThtoAll['FailList'][nconn]
            sim.net.cells[indexPost_Gid].conns[conn_index]['hObj'].syn().plas = Exc_ThtoAll['PlasList'][nconn]
            sim.net.cells[indexPost_Gid].conns[conn_index]['hObj'].syn().tau_plas = tau_plas
    
    # Cortex to cortex
    # Excitatory
    for nconn in range(len(Exc_AlltoAll['connList_gID'])):
        if nconn%100000 == 0 and sim.rank == 0: print("Exc Cortex->Cortex: ",nconn, "\n")

        Pre_Gid_rel = Exc_AlltoAll['connList_gID'][nconn][0]
        Post_Gid_rel = Exc_AlltoAll['connList_gID'][nconn][1]
        Pre_Gid = FirstCortex_Gid + Pre_Gid_rel
        Post_Gid = FirstCortex_Gid + Post_Gid_rel
        if Post_Gid in PostGids_nhost:
            indexPost_Gid = PostGids_nhost.index(Post_Gid)
            try:
                conn_index = list_preGids[indexPost_Gid].index(Pre_Gid)
            except:
                print("Error in connections pre-post Ids")
        
        # Setting individual parameters
        sim.net.cells[indexPost_Gid].conns[conn_index]['hObj'].syn().tau_rise = Exc_AlltoAll['Trise'][nconn]
        sim.net.cells[indexPost_Gid].conns[conn_index]['hObj'].syn().tau_fall = Exc_AlltoAll['Tfall'][nconn]
        sim.net.cells[indexPost_Gid].conns[conn_index]['hObj'].syn().cn = Exc_AlltoAll['CN'][nconn]
        sim.net.cells[indexPost_Gid].conns[conn_index]['hObj'].syn().mean_amp = Exc_AlltoAll['weightList'][nconn]
        sim.net.cells[indexPost_Gid].conns[conn_index]['hObj'].syn().cv = Exc_AlltoAll['CVList'][nconn]
        sim.net.cells[indexPost_Gid].conns[conn_index]['hObj'].syn().std0 = Exc_AlltoAll['STDList'][nconn]
        sim.net.cells[indexPost_Gid].conns[conn_index]['hObj'].syn().pf = Exc_AlltoAll['FailList'][nconn]
        sim.net.cells[indexPost_Gid].conns[conn_index]['hObj'].syn().plas = Exc_AlltoAll['PlasList'][nconn]
        sim.net.cells[indexPost_Gid].conns[conn_index]['hObj'].syn().tau_plas = tau_plas

    # Inhibitory
    for nconn in range(len(Inh_AlltoAll['connList_gID'])):
        if nconn%100000 == 0 and sim.rank == 0: print("Inh Cortex->Cortex: ",nconn, "\n")

        Pre_Gid_rel = Inh_AlltoAll['connList_gID'][nconn][0]
        Post_Gid_rel = Inh_AlltoAll['connList_gID'][nconn][1]
        Pre_Gid = FirstCortex_Gid + Pre_Gid_rel
        Post_Gid = FirstCortex_Gid + Post_Gid_rel
        if Post_Gid in PostGids_nhost:
            indexPost_Gid = PostGids_nhost.index(Post_Gid)
            try:
                conn_index = list_preGids[indexPost_Gid].index(Pre_Gid)
            except:
                print("Error in connections pre-post Ids")
        
        # Setting individual parameters
        sim.net.cells[indexPost_Gid].conns[conn_index]['hObj'].syn().tau_rise = Inh_AlltoAll['Trise'][nconn]
        sim.net.cells[indexPost_Gid].conns[conn_index]['hObj'].syn().tau_fall = Inh_AlltoAll['Tfall'][nconn]
        sim.net.cells[indexPost_Gid].conns[conn_index]['hObj'].syn().cn = Inh_AlltoAll['CN'][nconn]
        sim.net.cells[indexPost_Gid].conns[conn_index]['hObj'].syn().mean_amp = Inh_AlltoAll['weightList'][nconn]
        sim.net.cells[indexPost_Gid].conns[conn_index]['hObj'].syn().cv = Inh_AlltoAll['CVList'][nconn]
        sim.net.cells[indexPost_Gid].conns[conn_index]['hObj'].syn().std0 = Inh_AlltoAll['STDList'][nconn]
        sim.net.cells[indexPost_Gid].conns[conn_index]['hObj'].syn().pf = Inh_AlltoAll['FailList'][nconn]
        sim.net.cells[indexPost_Gid].conns[conn_index]['hObj'].syn().plas = Inh_AlltoAll['PlasList'][nconn]
        sim.net.cells[indexPost_Gid].conns[conn_index]['hObj'].syn().tau_plas = tau_plas
