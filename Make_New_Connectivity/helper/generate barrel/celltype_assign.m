function [M_parameter,cellinfo,tag]=celltype_assign(cellinfo,Ep,layer)
%assign cell type to cells, can be either layer 4 or layer 2/3 cells

% suffle cells
idx = randperm(size(cellinfo, 1));
cellinfo = cellinfo(idx, :);

layer=lower(layer);
switch layer
    case 'l4'
        
        if size(cellinfo,2)==3
            Ne=round(length(sum(cellinfo,2))*Ep);
            Ni=length(sum(cellinfo,2))-Ne;
            label=zeros(size(cellinfo,1),1);
            
            nl4=Ne+Ni;
            S_P=round(0.32*Ne);
            S_S=Ne-S_P;
            NeN=[S_P,S_S];
            
            FS=round(0.5*Ni);
            RSNP=Ni-FS;
            NiN=[FS,RSNP];
            
            N_cell=[NeN,NiN];
            for lp=1:length(N_cell)
                ind=[0,cumsum(N_cell)];
                label(ind(lp)+1:ind(lp+1))=lp;
            end
            cellinfo=[cellinfo,label];
        else
            Ne=length(find(cellinfo(:,4)==1))+length(find(cellinfo(:,4)==2));
            NiN=[length(find(cellinfo(:,4)==3)),length(find(cellinfo(:,4)==4))];
        end
        
        a=[0.01+0.02*rand(Ne,1);0.08+0.04*rand(NiN(1),1); 0.01+0.02*rand(NiN(2),1)];
        b=[0.3+0.05*rand(Ne,1);0.35-0.05*rand(NiN(1),1); 0.3+0.1*rand(NiN(2),1)];
        c=[-60+10*(rand(Ne,1)).^2; -65+10*rand(NiN(1),1); -65+10*rand(NiN(2),1)];
        d=[10+6*(rand(Ne,1)); 12+6*rand(NiN(1),1);8+4*rand(NiN(2),1)];
        
        M_parameter=[a,b,c,d];
        
        tag={'S_P','S_S','FS','RSNP'};
%         keyboard
    case 'l23'
        if size(cellinfo,2)==3
            Ne=round(size(cellinfo,1)*Ep);
            Ni=size(cellinfo,1)-Ne;
            
            Ne2=round(size(cellinfo,1)*Ep); Ni2=size(cellinfo,1)-Ne2;
            DBC_VIP=round(Ni2*0.15);
            Bip_VIP=round(Ni2*0.00);
            SBC_CR=round(Ni2*0.14);
            Bip_CR=round(Ni2*0.08);
            BS_PV=round(Ni2*0.25);
            CH_PV=round(Ni2*0.01);
            Mar_SOM=round(Ni2*0.15);
            Res_PV=round(Ni2*0.13);
            Bit_SOM=round(Ni2*0.00);
            NGC_AC = Ni2-sum([DBC_VIP,Bip_VIP,SBC_CR,Bip_CR,BS_PV,CH_PV,Mar_SOM,Res_PV,Bit_SOM]);
            
            NiN=[BS_PV,CH_PV,Res_PV,Mar_SOM,Bit_SOM,DBC_VIP,Bip_VIP,Bip_CR,SBC_CR,NGC_AC];
%             clear BS_PV CH_PV Res_PV Mar_SOM Bit_SOM DBC_VIP Bip_VIP Bip_CR SBC_CR NGC_AC
            NiN_types = {'BS_PV', 'CH_PV', 'Res_PV', 'Mar_SOM', 'Bit_SOM', 'DBC_VIP', ...
                'Bip_VIP', 'Bip_CR', 'SBC_CR', 'NGC_AC'};
            
            label=zeros(size(cellinfo,1),1);
            N_cell=[Ne,NiN];
            for lp=1:length(N_cell)
                ind=[0,cumsum(N_cell)];
                label(ind(lp)+1:ind(lp+1))=lp;
            end
            cellinfo=[cellinfo,label];
            
        else 
            Ne=length(find(cellinfo(:,4)==1));
            for lp=2:max(cellinfo(:,4))
                NiN(lp-1)=length(find(cellinfo(:,4)==lp));
            end
            BS_PV=NiN(1);
            CH_PV=NiN(2);
            Res_PV=NiN(3);
            Mar_SOM=NiN(4);
            Bit_SOM=NiN(5);
            DBC_VIP=NiN(6);
            Bip_VIP=NiN(7);
            Bip_CR=NiN(8);
            SBC_CR=NiN(9);
            NGC_AC=NiN(10);
        end
        re=rand(Ne,1); % re and ri is used to give some varibility to spiking pattern of neuron model
        ri=rand(sum(NiN),1); 

        a=[0.01+0.02*rand(Ne,1);0.08+0.04*rand(BS_PV+CH_PV,1); 0.01+0.02*rand(Res_PV,1); 0.01+0.03*rand(Mar_SOM,1).^1.5;0.01+0.02*rand(Bit_SOM+DBC_VIP+Bip_VIP+Bip_CR,1);0.01+0.03*rand(SBC_CR,1); 0.01+0.02*rand(NGC_AC,1)]; %
        b=[0.5+0.05*rand(Ne,1);0.9-0.05*rand(BS_PV+CH_PV,1); 0.5+0.1*rand(Res_PV,1).^1.5; 0.5+0.05*rand(Mar_SOM+Bit_SOM,1);0.55+0.05*rand(DBC_VIP+Bip_VIP+Bip_CR,1); 0.5*ones(SBC_CR,1);0.55+0.05*rand(NGC_AC,1).^2];
        c=[-65+15*re.^2; -65-10*rand(BS_PV+CH_PV,1); -60+10*rand(Res_PV,1).^0.9; -70+15*rand(Mar_SOM,1).^1.5;-70+15*rand(DBC_VIP+Bip_VIP+Bip_CR+Bit_SOM,1).^0.8;-65+10*rand(SBC_CR,1);-90*ones(NGC_AC,1)];
        d=[8+6*re; 10+6*rand(BS_PV+CH_PV,1);4+4*rand(Res_PV,1); 12-8*rand(Mar_SOM+Bit_SOM,1).^0.7; 4+6*rand(DBC_VIP+Bip_VIP+Bip_CR,1).^0.7;8+6*rand(SBC_CR,1).^0.7;14+4*rand(NGC_AC,1)];
        
        M_parameter=[a,b,c,d];
        
        % move NG type of cell to near the border of L2
        idx = find(cellinfo(:,4) == 11);
        idx = idx(1:round(0.7*length(idx)));
        cellinfo(idx, 3) = abs(10 + 50*randn(length(idx), 1));
        % move CR cells towards L2
        idx = find(cellinfo(:,4) == 10 & cellinfo(:,3) > 200);
        idx = idx(1:round(0.5*length(idx)));
        cellinfo(idx, 3) = cellinfo(idx, 3) - 200.*(cellinfo(idx, 3))/400;
        % adjust excitatory neuron distribution
        idx = find(cellinfo(:,4) == 1 & cellinfo(:,3) < 10);
        idx = idx(1:35);
        cellinfo(idx, 3) = max(cellinfo(:,3)) - cellinfo(idx, 3);
        idx = find(cellinfo(:,4) == 1 & cellinfo(:,3) < 20 & cellinfo(:,3) > 10);
        idx = idx(1:27);
        cellinfo(idx, 3) = max(cellinfo(:,3)) - cellinfo(idx, 3);
        idx = find(cellinfo(:,4) == 1 & cellinfo(:,3) < 30 & cellinfo(:,3) > 20);
        idx = idx(1:15);
        cellinfo(idx, 3) = max(cellinfo(:,3)) - cellinfo(idx, 3);
        idx = find(cellinfo(:,4) == 1 & cellinfo(:,3) < 40 & cellinfo(:,3) > 30);
        idx = idx(1:10);
        cellinfo(idx, 3) = max(cellinfo(:,3)) - cellinfo(idx, 3);
        
%         % move FsPV towards L3
%         idx = find(cellinfo(:, 4) == 2 & cellinfo(:,3) < 200);
%         idx = idx(round(0.1*length(idx)));
%         cellinfo(idx, 3) = cellinfo(idx, 3) + 100*rand(length(idx), 1);
        
        tag={'Pyr','FSBS','FSCH','BSPV','Mar','Bit','DBC','Bip','Bip','SBC','NG'};
end