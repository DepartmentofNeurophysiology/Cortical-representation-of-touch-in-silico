function CM = simConn_ThtoL23_ThtoAll(ADO_Mat, postIND, postType, AxonRatio, DendRatio, DMs)
% ADO_Mat: axon-dendritic overlap matrix
% postIND: indexing to seperate different type of post synaptic cells
% postType: numbers specify the type of post-synaptic cell
% keyboard
% connectivity matrix need to be calculated for each type of post synaptic
% cell
for i = 1:length(postIND) - 1
    temp = ADO_Mat(postIND(i)+1:postIND(i+1), :);

    ADO_temp = temp;
    DM = DMs(postIND(i)+1:postIND(i+1), :);
    DR = DendRatio(postIND(i)+1:postIND(i+1));

    switch postType(i)
        case(1) % post-synaptic cell is pyramidal
            CM(postIND(i)+1:postIND(i+1), :) = simConn_ThtoL23_ThtoPyr(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
        case{2,3} % post-synaptic cell is fast spiking PV+
            CM(postIND(i)+1:postIND(i+1), :) = simConn_ThtoL23_ThtoFsPV(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
        case(4) % post-synaptic cell is bursting PV+
            CM(postIND(i)+1:postIND(i+1), :) = simConn_ThtoL23_ThtoRsPV(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
        case{5, 6} % post-synaptic cell is MarSOM
            CM(postIND(i)+1:postIND(i+1), :) = simConn_ThtoL23_ThtoMarSOM(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
        case{7, 8} % post-synaptic cell is BipVIP
            CM(postIND(i)+1:postIND(i+1), :) = simConn_ThtoL23_ThtoBipVIP(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
        case(9) % post-synaptic cell is BipCR
            CM(postIND(i)+1:postIND(i+1), :) = simConn_ThtoL23_ThtoBipCR(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
        case(10) % post-synaptic cell in SbcCR
            CM(postIND(i)+1:postIND(i+1), :) = simConn_ThtoL23_ThtoSbcCR(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
        case(11) % post-synaptic cell in NG cells
            CM(postIND(i)+1:postIND(i+1), :) = simConn_ThtoL23_ThtoNG(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
    end
end

%% nested function
    function CMs = simConn_ThtoL23_ThtoPyr(ADO_temp, AR, DR)
        % map axon-dendrite overlapping index into connectivity matrix
        % the value is normalized to meausre average connection probability
        % from experiment; for Pyr-Pyr connection, within 100 micron p =
        % 0.72
        % also total number of connection is controlled by the axon and
        % dendrite ratio in the model region
        % assuming a total convengenc/divergence rate of 140 for total
        % connection
        
        pConn = ADO_temp/4.2454e-12*0.6;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMs = zeros(size(pConn));
        CMs(RI < pConn) = 1;
        
        CMs = conn_reduction(CMs, 120*AR', 35*DR');
    end


    function CMs = simConn_ThtoL23_ThtoFsPV(ADO_temp, AR, DR)
        % for Pyr-FsPV connection, within 50 micron p =
        % 0.18
        % assuming a total convengenc rate of 700 and divergence rate of 90
        
        pConn = ADO_temp/4.0471e-12*0.8;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMs = zeros(size(pConn));
        CMs(RI < pConn) = 1;
        
        CMs = conn_reduction(CMs, 60*AR', 45*DR');
    end


    function CMs = simConn_ThtoL23_ThtoRsPV(ADO_temp, AR, DR)
        % for Pyr-FsPV connection, within 100 micron p =
        % 0.5
        % assuming a total convengenc rate of 240 and divergence rate of 40
        
        pConn = ADO_temp/3.8930e-12*0.8;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMs = zeros(size(pConn));
        CMs(RI < pConn) = 1;
        
        CMs = conn_reduction(CMs, 60*AR', 35*DR');
    end


    function CMs = simConn_ThtoL23_ThtoMarSOM(ADO_temp, AR, DR)
        % for Pyr-MarSOM connection, no connection
        % assuming a total convengenc rate of 700 and divergence rate of 90
%         
        pConn = ADO_temp/5.9268e-12*0.4;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMs = zeros(size(ADO_temp));
        CMs(RI < pConn) = 1;
        
        CMs = conn_reduction(CMs, 40*AR', 30*DR');
    end
  
  function CMs = simConn_ThtoL23_ThtoBipVIP(ADO_temp, AR, DR)
    % for Pyr-MarSOM connection, no connection
    % assuming a total convengenc rate of 700 and divergence rate of 90
    %
    pConn = ADO_temp/4.9268e-12*0.8;
    % convert pConn into binary connectivity matrix
    RI = rand(size(pConn));
    CMs = zeros(size(ADO_temp));
    CMs(RI < pConn) = 1;
    
    CMs = conn_reduction(CMs, 40*AR', 40*DR');
  end

  function CMs = simConn_ThtoL23_ThtoBipCR(ADO_temp, AR, DR)
    % for Pyr-MarSOM connection, no connection
    % assuming a total convengenc rate of 700 and divergence rate of 90
    %
    pConn = ADO_temp/4.9268e-12*0.8;
    % convert pConn into binary connectivity matrix
    RI = rand(size(pConn));
    CMs = zeros(size(ADO_temp));
    CMs(RI < pConn) = 1;
    
    CMs = conn_reduction(CMs, 40*AR', 40*DR');
  end

  function CMs = simConn_ThtoL23_ThtoSbcCR(ADO_temp, AR, DR)
    % for Pyr-MarSOM connection, no connection
    % assuming a total convengenc rate of 700 and divergence rate of 90
    %
    pConn = ADO_temp/4.9268e-12*0.8;
    % convert pConn into binary connectivity matrix
    RI = rand(size(pConn));
    CMs = zeros(size(ADO_temp));
    CMs(RI < pConn) = 1;
    
    CMs = conn_reduction(CMs, 40*AR', 40*DR');
  end

  function CMs = simConn_ThtoL23_ThtoNG(ADO_temp, AR, DR)
    % for Pyr-MarSOM connection, no connection
    % assuming a total convengenc rate of 700 and divergence rate of 90
    %
    pConn = ADO_temp/4.9268e-12*0.8;
    % convert pConn into binary connectivity matrix
    RI = rand(size(pConn));
    CMs = zeros(size(ADO_temp));
    CMs(RI < pConn) = 1;
    
    CMs = conn_reduction(CMs, 40*AR', 40*DR');
  end

end

% %%
% %check connection probability as a function of distance
% bin = 0:25:500;
% N_all = [0, cumsum(NiN2)];
% Pconn = {};
% for post = 1:length(N_all) - 1
%     Pconn{post,1} = PConn_pairs(CM_modified(N_all(post)+1:N_all(post+1),:), ...
%         bin, DM(N_all(post)+1:N_all(post+1), :));
% end