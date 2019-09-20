function CM_old=conn_reduction(CM_old, lim_pre, lim_post)
%modifies the connectivity matrix based on the maximum number of pre/post
%synpatic input could exists
warning off;
%modify number of contacts one presynaptic neuron could have (divergence rate)
% keyboard
a=sum(CM_old,1);
b=lim_pre';
ind=find(a-b>0);%actual number exceeds the limitation
if isempty(find(ind, 1))==0
    %randomly chose the contact to remove
    for lp=1:length(ind)
        all=randperm(a(ind(lp)));
        N_red=floor(a(ind(lp))-b(ind(lp)));%number of connections need to be reduced
        old=CM_old(:,ind(lp));%original connectivity matrix
        ind_C=find(old);%originally connected pairs
        old(ind_C(all(1:N_red)))=0;%remove those contacts
        CM_old(:,ind(lp))=old;
    end
end

%now we do not need to worry about this
%modify number of contacts one postsynaptic neuron could have (convergence rate)
a=sum(CM_old,2);
b=lim_post;
ind=find(a-b>0);%actual number exceeds the limitation
if isempty(find(ind, 1))==0
    %randomly chose the contact to remove
    for lp=1:length(ind)
        all=randperm(a(ind(lp)));
        N_red=a(ind(lp))-b(ind(lp));%number of connections need to be reduced
        old=CM_old(ind(lp),:);%original connectivity matrix
        ind_C=find(old);%originally connected pairs
        old(ind_C(all(1:N_red)))=0;%remove those contacts
        CM_old(ind(lp),:)=old;
    end
end
% 