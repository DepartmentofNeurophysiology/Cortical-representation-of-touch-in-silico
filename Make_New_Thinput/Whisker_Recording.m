classdef Whisker_Recording < matlab.mixin.Copyable
    % Transform deflection angles, base angles, curvatures and grid
    % displacements into one another, given the assumption that the whisker
    % bends circularly (i.e. no deformations)
    % Use transfrom('what_to_transform_to'), it will either calculate it or give you
    % the quantities that are missing
    % Possible transformations:
    % * length_whisker
    % * grid_displacement
    % * base_angle
    % * deflection_angle
    % * curvature
    % * whisker_tip_position
    
    properties
        length_whisker                  % in mm
        distance_grid                   % in mm
        binsize
        grid_displacement_recordings    % Ntrace cell array with each an Ntimex1 array or NtracexNtime array with recordings 
        grid_displacement_unit          % should be mm
        grid_displacement_baseline      % baseline (rest) value
        grid_displacement_absrel        % whether recording is absolute or relative to baseline (cellstring or string)
        deflection_angle_recordings     % Ntrace cell array with each an Ntimex1 array or NtracexNtime array with recordings
        deflection_angle_unit           % 'degree' or 'radian'
        deflection_angle_baseline  
        deflection_angle_absrel  
        base_angle_recordings           % Ntrace cell array with each an Ntimex1 array or NtracexNtime array with recordings
        base_angle_unit                 % 'degree' or 'radian'
        base_angle_baseline
        base_angle_absrel
        curvature_recordings            % Ntrace cell array with each an Ntimex1 array or NtracexNtime array with recordings
        curvature_unit                  % should be 'mm-1'
        curvature_baseline
        curvature_absrel
        whisker_tip_position_recordings % Ntrace cell array with each an Ntimex2 array or NtracexNtimex2 array with recordings (NB 1: s (horizontal), 2: x (vertical))
        whisker_tip_position_unit       % Should be 1 value, same for both!
        whisker_tip_position_baseline 
        whisker_tip_position_absrel     % Should be 1 value, same for both!
        
        

    end
    
    properties (Dependent, Hidden)

    end
    
    methods
        %% Generator
        function b = Whisker_Recording()
        end
        
        %% Transformations
        function b = transform(b,outputfieldname, delete_original)
            if nargin == 2
                delete_original = 0;
            end
            if strcmp(outputfieldname, 'length_whisker')
                outputrec = outputfieldname;
            else
                outputrec = [outputfieldname '_recordings'];
            end
            if ~isempty(b.(outputrec)) 
                if strcmp(outputfieldname, 'length_whisker')
                    okornot = 1;
                else
                    [~, ~, okornot] = check_unit(b, outputfieldname);
                end
                if okornot
                    disp('No calculation needed')
                end
            elseif ~isempty(b.distance_grid)
                disp('Assuming fixed grid distance, so grid displacement experiment')
                    inputstruct.distance_grid = b.distance_grid;
                    if ~isempty(b.length_whisker)
                        inputstruct.length_whisker = b.length_whisker;
                        if ~isempty(b.base_angle_recordings)
                            inputfieldname = 'base_angle';
                        elseif ~isempty(b.deflection_angle_recordings)
                            inputfieldname = 'deflection_angle';
                        elseif ~isempty(b.curvature_recordings)
                            inputfieldname = 'curvature';
                        elseif ~isempty(b.grid_displacement_recordings)
                            inputfieldname = 'grid_displacement';
                        else
                            error('Two sets of recordings or 1 set and the whisker length are needed to do this calculation');
                        end
                        [~, ~, okornot] = check_unit(b, inputfieldname);
                        if okornot
                            disp(['Calculating ' outputfieldname ' from ' inputfieldname ' and whisker length']);
                            function_handle = ['Whisker_Recording.' inputfieldname '_to_' outputfieldname '_length_whisker_fixed_grid'];
                            b = calculate_recordings(b, inputfieldname, outputfieldname, inputstruct, function_handle, delete_original);
                        end
                    else
                        if ~isempty(b.base_angle_recordings)
                            inputfieldname1 = 'base_angle';
                            if ~isempty(b.deflection_angle_recordings)
                                inputfieldname2 = 'deflection_angle';
                            elseif ~isempty(b.curvature_recordings)
                                inputfieldname2 = 'curvature';
                            elseif ~isempty(b.grid_displacement_recordings)
                                inputfieldname2 = 'grid_displacement';
                            elseif strcmp(outputfieldname, 'grid_displacement')
                                disp(['No whisker length given, no deflection angle or curvature recordings given, approximating by ' outputfieldname ' = distance grid * tan' inputfieldname1 ])
                                [~, ~, okornot] = check_unit(b, inputfieldname1);
                                if okornot
                                    b = base_angle_and_grid_displacement_approximation(b, inputfieldname1, outputfieldname);
                                end
                                return
                            else
                                error('Give 2 sets of recordings, or 1 set and whisker length')
                            end
                            
                        elseif ~isempty(b.deflection_angle_recordings)
                            inputfieldname1 = 'deflection_angle';
                            if ~isempty(b.curvature_recordings)
                                inputfieldname2 = 'curvature';
                            elseif ~isempty(b.grid_displacement_recordings)
                                inputfieldname2 = 'grid_displacement';
                            else
                                error('Give 2 sets of recordings, or 1 set and whisker length')
                            end
                        elseif ~isempty(b.curvature_recordings)
                            inputfieldname1 = 'curvature';
                            if ~isempty(b.grid_displacement_recordings)
                                inputfieldname2 = 'grid_displacement';
                            end
                        elseif strcmp(outputfieldname, 'base_angle') && ~isempty(b.grid_displacement_recordings)
                            inputfieldname = 'grid_displacement';
                            [~, ~, okornot] = check_unit(b, inputfieldname);
                            disp(['No whisker length given, no deflection angle or curvature recordings given, approximating by ' outputfieldname ' = atan(' inputfieldname '/distance_grid)'])
                            if okornot 
                                b = base_angle_and_grid_displacement_approximation(b, inputfieldname, outputfieldname);
                            end
                            return
                        else
                            error('Two sets of recordings or 1 set and the whisker length are needed to do this calculation');
                        end
                        [~, ~, okornot1] = check_unit(b, inputfieldname1);
                        [~, ~, okornot2] = check_unit(b, inputfieldname2);
                        if okornot1 && okornot2
                            disp(['Calculating ' outputfieldname ' from ' inputfieldname1 ' and ' inputfieldname2]) 
                            function_handle = ['Whisker_Recording.' inputfieldname1 '_to_' outputfieldname '_' inputfieldname2 '_fixed_grid'];
                            if strcmp(outputfieldname, 'length_whisker')
                                b = calculate_length_whisker(b, inputfieldname1, inputfieldname2, inputstruct, function_handle);
                            else
                                b = calculate_recordings_2to1(b, inputfieldname1, inputfieldname2, outputfieldname, inputstruct, function_handle, delete_original);
                            end
                        end
                    end
            elseif ~isempty(b.length_whisker)
                % D = lw*c/2
                inputstruct.length_whisker = b.length_whisker; 
                if strcmp(outputfieldname, 'whisker_tip_position')
                    if ~isempty(b.(outputrec)) 
                        [~, ~, okornot] = check_unit(b, outputfieldname);
                        if okornot
                            disp('No calculation needed')
                        end
                    elseif ~isempty(b.base_angle_recordings)
                        inputfieldname1 = 'base_angle';
                        if  ~isempty(b.deflection_angle_recordings)
                            inputfieldname2 = 'deflection_angle';
                        elseif  ~isempty(b.curvature_recordings)
                            inputfieldname2 = 'curvature';
                        else
                            error('Deflection angles or curvatures are needed for this calculation')
                        end
                        [~, ~, okornot1] = check_unit(b, inputfieldname1);
                        [~, ~, okornot2] = check_unit(b, inputfieldname2);
                        if okornot1 && okornot2
                            disp(['Calculating ' outputfieldname ' from ' inputfieldname1 ' and ' inputfieldname2]) 
                            function_handle = ['Whisker_Recording.' inputfieldname1 '_to_' outputfieldname '_' inputfieldname2 '_length_whisker'];
                            b = calculate_recordings_2to1(b, inputfieldname1, inputfieldname2, outputfieldname, inputstruct, function_handle, delete_original);
                        end
                    else
                        error('Base angles are needed for this calculation')
                    end
                elseif strcmp(outputfieldname, 'length_whisker') || strcmp(outputfieldname, 'grid_displacement')
                    error('This calculation is not possible')
                else
                    if ~isempty(b.deflection_angle_recordings) && strcmp(outputfieldname, 'curvature')
                        inputfieldname = 'deflection_angle';
                    elseif ~isempty(b.curvature_recordings) && strcmp(outputfieldname, 'deflection_angle')
                        inputfieldname = 'curvature';
                    elseif ~isempty(b.whisker_tip_position_recordings)
                        if ~isempty(b.deflection_angle_recordings)
                            ip = input('Do you want to calculate the base angle using deflection angles or whisker length? (da / wl)','s');
                            if strcmp(ip, 'wl')
                                inputfieldname = 'whisker_tip_position';
                            elseif strcmp(ip, 'da')
                                inputfieldname1 = 'deflection_angle';
                                inputfieldname2 = 'whisker_tip_position';
                                inputstruct = [];
                                [~, ~, okornot1] = check_unit(b, inputfieldname1);
                                [~, ~, okornot2] = check_unit(b, inputfieldname2);
                                if okornot1 && okornot2
                                    disp(['Calculating ' outputfieldname ' from ' inputfieldname1 ' and ' inputfieldname2]);
                                    function_handle = ['Whisker_Recording.' inputfieldname1 '_to_' outputfieldname '_' inputfieldname2];
                                    b = calculate_recordings_2to1(b, inputfieldname1, inputfieldname2, outputfieldname, inputstruct, function_handle, delete_original);
                                end  
                                return
                            else
                                error('Please choose whisker length or deflection angle')
                            end
                        else
                            inputfieldname = 'whisker_tip_position';
                        end
                    else
                        error('This calculation is only possible for a fixed grid')
                    end
                    [~, ~, okornot] = check_unit(b, inputfieldname);
                    if okornot
                        disp(['Calculating ' outputfieldname ' from ' inputfieldname ' and whisker length']);
                        function_handle = ['Whisker_Recording.' inputfieldname '_to_' outputfieldname '_length_whisker'];
                        b = calculate_recordings(b, inputfieldname, outputfieldname, inputstruct, function_handle, delete_original);
                    end
                end
            elseif (strcmp(outputfieldname, 'length_whisker') && ~isempty(b.curvature_recordings)) && ~isempty(b.deflection_angle_recordings)
                % lw = 2c/D
                inputfieldname1 = 'deflection_angle';
                inputfieldname2 = 'curvature';
                inputstruct = [];
                [~, ~, okornot1] = check_unit(b, inputfieldname1);
                [~, ~, okornot2] = check_unit(b, inputfieldname2);
                if okornot1 && okornot2
                    disp(['Calculating ' outputfieldname ' from ' inputfieldname1 ' and ' inputfieldname2]);
                    function_handle = ['Whisker_Recording.' inputfieldname1 '_to_' outputfieldname '_' inputfieldname2 '_fixed_grid'];
                    b = calculate_length_whisker(b, inputfieldname1, inputfieldname2, outputfieldname, inputstruct, function_handle);
                end
            elseif (~isempty((b.whisker_tip_position_recordings)) && ~isempty(b.deflection_angle_recordings)) && strcmp(outputfieldname, 'base_angle')
                inputfieldname1 = 'deflection_angle';
                inputfieldname2 = 'whisker_tip_position';
                inputstruct = [];
                [~, ~, okornot1] = check_unit(b, inputfieldname1);
                [~, ~, okornot2] = check_unit(b, inputfieldname2);
                if okornot1 && okornot2
                    disp(['Calculating ' outputfieldname ' from ' inputfieldname1 ' and ' inputfieldname2]);
                    function_handle = ['Whisker_Recording.' inputfieldname1 '_to_' outputfieldname '_' inputfieldname2];
                    b = calculate_recordings_2to1(b, inputfieldname1, inputfieldname2, outputfieldname, inputstruct, function_handle, delete_original);
                end                
            else
                disp('This calculation is not possible')
            end
        end
                
        %% Absolute and relative
        function b = abs_to_rel(b, fieldname)
            disp(['Making ' fieldname ' relative'])
            rec    = [fieldname '_recordings'];
            baslin = [fieldname '_baseline'];
            absrel = [fieldname '_absrel'];
            [~, Ntrace] = size(b.(rec));
            
            if ~isempty(b.(baslin))
                for nt = 1:Ntrace
                    b.(rec){nt} = b.(rec){nt} - b.(baslin).*ones(size(b.(rec){nt}));
                end
                b.(absrel) = 'rel';
            else
                error(['Please give baseline value for ' fieldname])
            end
        end
        
        function b = rel_to_abs(b, fieldname)
            disp(['Making ' fieldname ' absolute'])
            rec    = [fieldname '_recordings'];
            baslin = [fieldname '_baseline'];
            absrel = [fieldname '_absrel'];
            [~, Ntrace] = size(b.(rec));
            if ~isempty(b.(baslin))
                for nt = 1:Ntrace
                    b.(rec){nt} = b.(rec){nt} + b.(baslin).*ones(size(b.(rec){nt}));
                end
                b.(absrel) = 'abs';            
            else
                error(['Please give baseline value for ' fieldname])
            end
        end
        
        %% Units
        function [unit, unit_target, okornot] = check_unit(b, fieldname)
            rec  = [fieldname '_recordings'];
            unit = [fieldname '_unit'];
            baslin = [fieldname '_baseline'];
            % Check if recordings are in a proper cell array
            if iscell(b.(rec))
                [X, Ntrace] = size(b.(rec));
                if (Ntrace == 1) && (X>1)
                    b.(rec) = b.(rec)';
                    [~, Ntrace] = size(b.(rec));
                end
            else
                ND = ndims(b.(rec));
                if ND == 1
                    Ntrace = 1;
                    Ntime = 1;
                elseif ND == 2
                    [Ntrace, Ntime] = size(b.(rec));
                elseif ND == 3
                    [Ntrace, Ntime, ~] = size(b.(rec));
                end
                
                % Put into cells
                temp = cell(1,Ntrace);
                for nt = 1:Ntrace
                    if ND == 1 || ND == 2
                        temp{nt} = b.(rec)(nt,:)';
                    else
                        temp{nt} = squeeze(b.(rec)(nt,:, :));
                    end
                end
                b.(rec) = temp;
            end
            
            % Check need to transpose
            for nt = 1:Ntrace
                [Ntime, Ndim] = size(b.(rec){nt});
                if Ndim>Ntime
                    ip = input(['Number of dimensions of ', fieldname, ' = ' num2str(Ndim) ', this larger than number of time steps(' num2str(Ntime) '). Transpose? (y/n)'] , 's');
                    if strcmp(ip, 'y')
                        b.(rec){nt} = b.(rec){nt}';
                        Ntt = Ntime;
                        Ntime = Ndim;
                        Ndim = Ntt;
                    end
                end
            end
            
            % Check units
            unit_target = Whisker_Recording.set_unitname(fieldname);
            [~, ndimtarget] = size(unit_target);
            if ~isempty(b.(unit))
                b.(unit) = cellstr(b.(unit));
                okornot = nan*ones(1,ndimtarget);
                for nd = 1:ndimtarget
                    if strcmp(b.(unit){nd}, unit_target{nd})
                        % Everything ok
                        okornot(nd) = 1;
                    elseif strcmp(b.(unit){nd}, 'degree') && strcmp(unit_target{nd}, 'radian')
                        disp('Recordings are in degrees, convert to radians')
                        for nt = 1:Ntrace
                            b.(rec){nt}     = Whisker_Recording.degree_to_radian(b.(rec){nt});
                        end
                        b.(baslin)  = Whisker_Recording.degree_to_radian(b.(baslin));
                        b.(unit){nd} = 'radian';
                        okornot(nd) = 1;
                    else
                        okornot(nd) = 0;
                        error(['Please give ' fieldname ' recordings in ' unit_target{nd}])
                    end
                end
            else
                disp([ 'Warning: ' fieldname ' unit is empty'])
                okornot = 1;
            end
          
        end
        
        %% Calculate new recordings
        function b = calculate_length_whisker(b, inputfieldname1, inputfieldname2, inputstruct, function_handle)
            inputrec1     = [inputfieldname1 '_recordings'];
            inputrec2     = [inputfieldname2 '_recordings'];
            
            % Make sure both recordings are absolute
            for n = 1:2    
                if n==1
                    fn = inputfieldname1;
                elseif n==2
                    fn = inputfieldname2;
                end
                inputbaslin  = [fn '_baseline']; % base angle
                inputabsrel  = [fn '_absrel'];
                if strcmp(b.(inputabsrel), 'rel')
                    if ~isempty(b.(inputbaslin))
                        b = rel_to_abs(b, inputfieldname);
                    else
                        error(['Please give baseline or absolute value for ' fn])
                    end
                end
            end    
            
            % Calculate length whisker
            if iscell(b.(inputrec1))
                [X1, Ntrace1] = size(b.(inputrec1));
                if (Ntrace1 == 1) && (X1>1)
                    b.(inputrec1) = b.(inputrec1)';
                    [~, Ntrace1] = size(b.(inputrec1));
                end
            else
                [Ntrace1, Ntime1] = size(b.(inputrec1));
                ip1temp = cell(1,Ntrace1);
                for nt = 1:Ntrace1
                    ip1temp{nt} = b.(inputrec1)(nt,:)';
                end
                b.(inputrec1) = ip1temp;
            end
            if iscell(b.(inputrec2))
                [X2, Ntrace2] = size(b.(inputrec2));
                if (Ntrace2 == 1) && (X2>1)
                    b.(inputrec2) = b.(inputrec2)';
                    [~, Ntrace2] = size(b.(inputrec2));
                end
            else
                [Ntrace2, Ntime2] = size(b.(inputrec2));
                ip2temp = cell(1,Ntrace2);
                for nt = 1:Ntrace2
                    ip2temp{nt} = b.(inputrec2)(nt,:)';
                end
                b.(inputrec2) = ip2temp;
            end

            if ~(Ntrace1 == Ntrace2)
                error(['Make sure ' inputfieldname1 ' recordings and ' inputfieldname2 ' recordings have same number of trials'])
            end
            wl = nan*ones(Ntrace1,1);
            for nt = 1:Ntrace1
                [Ntime1, ~] = size(b.(inputrec1){nt});
                [Ntime2, ~] = size(b.(inputrec2){nt});
                if ~(Ntime1 == Ntime2)
                    error(['Make sure ' inputfieldname1 ' recordings and ' inputfieldname2 ' recordings have same number of time steps'])
                end
                wl(nt) = feval(function_handle,b.(inputrec1){nt}, b.(inputrec2){nt}, inputstruct);
            end
            f = figure;
            plot(1:Ntrace1, wl, 'o')
            xlabel('Recording nr')
            ylabel('Whisker length')
            xlim([0 Ntrace1+1])
            ylim([0 ceil(max(wl))])
            ip = input('Calculated whisker length ok, take mean over recordings? (y/n)','s');
            close (f)
            if strcmp(ip,'y')
                b.length_whisker = nanmean(wl);
            else
                b.length_whisker = nan;
            end
        end
        
        function b = calculate_recordings(b, inputfieldname, outputfieldname, inputstruct, function_handle, delete_original)
            if nargin == 5
                delete_original = 0;
            end
            % Find out whether source recordings are relative or absolute,
            % calculate target 'recording' and set relevant properteis (unit, absrel, baseline) 
            inputbaslin  = [inputfieldname '_baseline'];
            inputabsrel  = [inputfieldname '_absrel']; 
            inputrec     = [inputfieldname '_recordings'];
            inputunit    = [inputfieldname '_unit'];
            outputabsrel = [outputfieldname '_absrel']; 
            outputunit   = [outputfieldname,'_unit'];
            
            % new units
            unitname = Whisker_Recording.set_unitname(outputfieldname);
            b.(outputunit)   = unitname;
            
            b.(outputabsrel) = b.(inputabsrel);
            
            % calculate recordings and baseline if it exists
            if strcmp(b.(inputabsrel), 'abs')
                % Absolute  
                b  = absolute_calculation(b, inputfieldname, outputfieldname, inputstruct, function_handle);
            elseif strcmp(b.(inputabsrel), 'rel')
                % Relative 
                if ~isempty(b.(inputbaslin)) 
                    b = relative_calculation(b, inputfieldname, outputfieldname, inputstruct, function_handle);
                else
                    error(['Please specify ' inputfieldname ' baseline'])
                end
            else
                error(['Please specify whether ' inputfieldname ' recordings are absolute or relative to a baseline'])
            end
            if delete_original
                disp('Deleting original recordings')
                b.(inputbaslin) = [];
                b.(inputabsrel) = [];
                b.(inputrec)    = [];
                b.(inputunit)   = [];
            end
        end
        
        function b = calculate_recordings_2to1(b, inputfieldname1, inputfieldname2, outputfieldname, inputstruct, function_handle, delete_original)
            if nargin == 6
                delete_original = 0;
            end
            
            inputbaslin1  = [inputfieldname1 '_baseline']; 
            inputabsrel1  = [inputfieldname1 '_absrel'];
            inputrec1     = [inputfieldname1 '_recordings'];
            inputunit1    = [inputfieldname1 '_unit'];
            inputbaslin2  = [inputfieldname2 '_baseline']; 
            inputabsrel2  = [inputfieldname2 '_absrel'];
            inputrec2     = [inputfieldname2 '_recordings'];
            inputunit2    = [inputfieldname1 '_unit'];
            outputabsrel  = [outputfieldname '_absrel']; 
            outputunit    = [outputfieldname,'_unit'];
            outputrec     = [outputfieldname '_recordings'];
            
            % Make both absolute or both relative
            if ~strcmp(b.(inputabsrel1), b.(inputabsrel2))
                if strcmp('rel', b.(inputabsrel2))
                    if ~isempty(b.(inputbaslin1))
                        b = abs_to_rel(b, inputfieldname1);
                    elseif ~isempty(b.deflection_angle_baseline)
                        b = rel_to_abs(b, inputfieldname2);
                    else
                        error([inputfieldname1 ' or ' inputfieldname1 ' baseline value needed'])
                    end
                elseif strcmp('abs', b.(inputabsrel2))
                    if ~isempty(b.(inputbaslin1))
                        b = rel_to_abs(b,inputfieldname1);
                    elseif ~isempty(b.(inputbaslin2))
                        b = abs_to_rel(b, inputfieldname2);
                    else
                        error('Baseline value needed')
                    end
                end
            end
            
            % Make output recordings and settings
            if iscell(b.(inputrec1))
                [~, Ntrace1] = size(b.(inputrec1));
            else
                [Ntrace1, Ntime1] = size(b.(inputrec1));
                ip1temp = cell(1,Ntrace1);
                for nt = 1:Ntrace1
                    ip1temp{nt} = b.(inputrec1)(nt,:)';
                end
                b.(inputrec1) = ip1temp;
            end
            if iscell(b.(inputrec2))
                [~, Ntrace2] = size(b.(inputrec2));
            else
                [Ntrace2, Ntime2] = size(b.(inputrec2));
                ip2temp = cell(1,Ntrace2);
                for nt = 1:Ntrace2
                    ip2temp{nt} = b.(inputrec2)(nt,:)';
                end
                b.(inputrec2) = ip2temp;
            end

            if ~(Ntrace1 == Ntrace2)
                error(['Make sure ' inputfieldname1 ' recordings and ' inputfieldname2 ' recordings have same number of trials'])
            end

            unitname = Whisker_Recording.set_unitname(outputfieldname);
            b.(outputabsrel) = b.(inputabsrel2);                
            
            b.(outputunit)= unitname;
            b.(outputrec) = cell(1,Ntrace1);
            if ~isempty(b.(inputbaslin1)) && ~isempty(b.(inputbaslin2))
                b.(outputbaslin) = feval(function_handle,b.(inputbaslin1),b.(inputbaslin2), inputstruct);
            end
            if strcmp(b.(inputabsrel2), 'abs')
                % Absolute calculation
                for nt = 1:Ntrace1
                    [Ntime1, ~] = size(b.(inputrec1){nt});
                    [Ntime2, ~] = size(b.(inputrec2){nt});
                    if ~(Ntime1 == Ntime2)
                        error(['Make sure ' inputfieldname1 ' recordings and ' inputfieldname2 ' recordings have same number of time steps'])
                    end
                    b.(outputrec){nt} = feval(function_handle, b.(inputrec1){nt},b.(inputrec2){nt}, inputstruct);
                end

            elseif strcmp(b.(inputabsrel2), 'rel')
                % Relative calculation
                if (~isempty(b.(inputbaslin2)) && ~isempty(b.(inputbaslin1))) && ~isempty(b.(outputbaslin))
                    for nt = 1:Ntrace1
                        % Absolute input
                        input1_abs = b.(inputrec1){nt}+b.(inputbaslin1);
                        input2_abs = b.(inputrec2){nt}+b.(inputbaslin2);
                        % Absolute output
                        output_abs = feval(function_handle, input1_abs, input2_abs, inputstruct);
                        % Relative output
                        b.(outputrec){nt} = output_abs - b.(outputbaslin);
                    end
                else
                    error(['Specify ' inputfieldname1 ' and ' inputfieldname1 ' baseline'])
                end
            else
                error('Please specify whether recordings are absolute or relative to a baseline')
            end
                        
            
            if delete_original
                disp('Deleting original recordings')
                b.(inputbaslin1) = [];
                b.(inputabsrel1) = [];
                b.(inputrec1)    = [];
                b.(inputunit1)   = [];
                b.(inputbaslin2) = [];
                b.(inputabsrel2) = [];
                b.(inputrec2)    = [];
                b.(inputunit2)   = [];
            end
        end

        function b = absolute_calculation(b, inputfieldname, outputfieldname, inputstruct, function_handle)
            % Calculate absolute target 'recording'
            inputrec        = [inputfieldname  '_recordings'];
            inputbaslin     = [inputfieldname  '_baseline'];
            outputrec       = [outputfieldname '_recordings'];
            outputbaslin    = [outputfieldname '_baseline'];

            if iscell(b.(inputrec))
                [X1, Ntrace] = size(b.(inputrec));
                if (Ntrace == 1) && (X1>1)
                    b.(inputrec) = b.(inputrec)';
                    [~, Ntrace] = size(b.(inputrec));
                end
            else
                [Ntrace, Ntime] = size(b.(inputrec));
                iptemp = cell(1,Ntrace);
                for nt = 1:Ntrace
                    iptemp{nt} = b.(inputrec)(nt,:)';
                end
                b.(inputrec) = iptemp;
            end
            
            % Calculate recordings
            for nt = 1:Ntrace
                try
                    b.(outputrec){nt} = feval(function_handle,b.(inputrec){nt}, inputstruct);
                catch
                    keyboard
                end
            end
            
            % Calculate baseline
            if ~isempty(b.(inputbaslin))
                b.(outputbaslin) = feval(function_handle,b.(inputbaslin), inputstruct);
            end

        end
        
        function b = relative_calculation(b, inputfieldname, outputfieldname, inputstruct, function_handle)
            % Calculate relative target 'recording'
            inputrec        = [inputfieldname  '_recordings'];
            inputbaslin     = [inputfieldname  '_baseline'];
            outputrec       = [outputfieldname '_recordings'];
            outputbaslin    = [outputfieldname '_baseline'];
            
            if iscell(b.(inputrec))
                [X1, Ntrace] = size(b.(inputrec));
                if (Ntrace == 1) && (X1>1)
                    b.(inputrec) = b.(inputrec)';
                    [~, Ntrace] = size(b.(inputrec));
                end
            else
                [Ntrace, Ntime] = size(b.(inputrec));
                iptemp = cell(1,Ntrace);
                for nt = 1:Ntrace
                    iptemp{nt} = b.(inputrec)(nt,:)';
                end
                b.(inputrec) = iptemp;
            end
            
            % Calculate baseline
            b.(outputbaslin)= feval(function_handle,b.(inputbaslin), inputstruct);
            
            % Calculate recordings
            for nt = 1:Ntrace     
                % Absolute input
                input_abs = b.(inputrec){nt}+b.(inputbaslin);
                % Absolute output
                output_abs = feval(function_handle,input_abs, inputstruct);
                % Relative output
                b.(outputrec){nt} = output_abs - b.(outputbaslin); 
            end
        end
        
        %% Helper functions 
        function b = base_angle_and_grid_displacement_approximation(b, inputfieldname, outputfieldname, delete_original)
            if nargin == 3
                delete_original = 0;
            end
            changed = 0;
            
            % x = s tan alpha or delta x = s delta alpha
            inputbaslin  = [inputfieldname '_baseline'];
            inputabsrel  = [inputfieldname '_absrel']; 
            
           
            
            
            inputstruct.distance_grid  = b.distance_grid;
            
            
            if (~isempty(b.(inputbaslin)) && strcmp(b.(inputabsrel), 'rel') ) || strcmp(b.(inputabsrel), 'abs')
                function_handle = ['Whisker_Recording.' inputfieldname '_to_' outputfieldname '_fixed_grid'];
            else
                 % change absrel for calculation
                if ~strcmp(b.(inputabsrel), 'abs')
                    disp('Ignoring that no baseline is known')
                    b.(inputabsrel) = 'abs';
                    changed = 1;
                end
                function_handle = ['Whisker_Recording.' inputfieldname '_to_' outputfieldname '_rel_nobaseline'];
            end
            b = calculate_recordings(b, inputfieldname, outputfieldname, inputstruct, function_handle, delete_original);
            
            % change absrel back
            if changed
                b.(inputabsrel) = 'rel';
                b.(inputabsrel) = 'abs';
            end
            
        end 
        
        %% Plot
        function b = plot(b, fieldname, traces)  
            hold all
            rec        = [fieldname  '_recordings'];
            if ~isempty(b.(rec))
                if iscell(b.(rec))
                    [X1, Ntrace] = size(b.(rec));
                    if (Ntrace == 1) && (X1>1)
                        b.(rec) = b.(rec)';
                        [~, Ntrace] = size(b.(rec));
                    end
                    if nargin == 2
                        traces = 1:Ntrace;
                    end
                    for nt = traces
                        [Ntime, Ndim] = size(b.(rec){nt});
                        if Ndim>Ntime
                            ip = input(['Number of dimensions of ', fieldname, ' = ' num2str(Ndim) ', this larger than number of time steps(' num2str(Ntime) '). Transpose? (y/n)'] , 's');

                            if strcmp(ip, 'y')
                                b.(rec){nt} = b.(rec){nt}';
                                Ntt = Ntime;
                                Ntime = Ndim;
                                Ndim = Ntt;
                            end
                        end

                        for nd = 1:Ndim
                            if ~isempty(b.binsize)
                                plot((1:Ntime)*b.binsize, b.(rec){nt}(:,nd))
                            else
                                plot(b.(rec){nt}(:,nd))
                            end
                        end
                    end
                else
                    ND = ndims(b.(rec));
                    if ND == 1
                        error('Nothing to plot')
                    elseif ND == 2
                        [Ntrace, Ntime] = size(b.(rec));
                    elseif ND == 3
                        [Ntrace, Ntime, ~] = size(b.(rec));
                    end
                    if Ntrace>Ntime
                        ip = input(['Number of recordings of ', fieldname, ' = ' num2str(Ntrace) ', this larger than number of time steps(' num2str(Ntime) '). Transpose? (y/n)'] , 's');

                        if strcmp(ip, 'y')
                            b.(rec) = b.(rec)';
                            Ntt = Ntime;
                            Ntime = Ntrace;
                            Ntrace = Ntt;
                        end
                    end
                    if nargin == 2
                        traces = 1:Ntrace;
                    end
                    for nt = traces
                        if ND == 2
                            wtp = b.(rec)(nt,:);
                        elseif ND == 3
                            wtp = squeeze(b.(rec)(nt,:,:));
                        end
                        if ~isempty(b.binsize)
                            plot((1:Ntime)*b.binsize, wtp)
                        else
                            plot(wtp)
                        end
                    end
                end
            else
                error(['Recording ' fieldname ' is empty, nothing to plot'])
            end
        end

        
    end
    
    
    methods (Static)
    %% Calculations
        % NB Note that not everything is internally consistent, due to
        % 1) Numerical approximations & assumptions ther
        % 2) Deformations that are ignored (whisker is not circular)
        % 3) Therefore: the recorded combination base angle and deflection
        % angle/curvature is not always consistent 
        % So the order in which calculations are done matters!
        %% To grid displacement
        
        % from base angle and whisker length (numerical approximation)
        function grid_displacement = base_angle_to_grid_displacement_length_whisker_fixed_grid(base_angle, inputstruct)
%             disp('Calculating grid displacement from base angle and whisker length')
            % NB the following function uses fixed s and lw
            deflection_angle = Whisker_Recording.base_angle_to_deflection_angle_length_whisker_fixed_grid(base_angle, inputstruct);
            grid_displacement = Whisker_Recording.base_angle_to_grid_displacement_deflection_angle_fixed_grid(base_angle, deflection_angle, inputstruct);
        end
        
        % from deflection angle and whisker length
        function grid_displacement = deflection_angle_to_grid_displacement_length_whisker_fixed_grid(deflection_angle, inputstruct)
%             disp('Calculating grid displacement from deflection angle and whisker length')
            % Absolute deflection angles 
            disp('Assuming positive grid displacement (2 solutions possible from deflection angle and whisker length)')
            grid_displacement = sqrt((inputstruct.length_whisker.*sin(deflection_angle)./deflection_angle).^2-inputstruct.distance_grid^2);
            grid_displacement(deflection_angle == 0) = sqrt((inputstruct.length_whisker).^2-inputstruct.distance_grid^2);
        end
        
        % from curvature and whisker length
        function grid_displacement = curvature_to_grid_displacement_length_whisker_fixed_grid(curvature, inputstruct)
%             disp('Calculating grid displacement from curvature and whisker length')
            % Absolute curvature
            grid_displacement = sqrt((2*sin(inputstruct.length_whisker.*curvature/2)./curvature).^2-inputstruct.distance_grid^2);
        end
                
        % from base angle and deflection angle
        function grid_displacement = base_angle_to_grid_displacement_deflection_angle_fixed_grid(base_angle, deflection_angle, inputstruct)
%             disp('Calculating grid displacement from base angle and deflection angle')
            % Absolute deflection angles and base angles
            grid_displacement = inputstruct.distance_grid*tan(deflection_angle+base_angle);
            disp('Note that in a freely whisking case, a whisker can make movements that are impossible in a grid experiment')
            lw = inputstruct.distance_grid*deflection_angle./(sin(deflection_angle).*cos(deflection_angle + base_angle));
            f = figure;
            subplot(2,1,1)
            plot(lw , '.')
            title('Whisker length (mm)')
            subplot(2,1,2)
            plot(base_angle)
            hold all
            plot(deflection_angle)
            legend('base angle','deflection angle')
            disp('Check if fixed whisker length is not violated')
            ip = input('Everything ok? (y/n)', 's');
            if strcmp(ip, 'y')
            else
                error('Constant whisker length violated. Try transforming to whisker tip position ')
            end
            close(f)
        end
  
        % from base angle and curvature
        function grid_displacement = base_angle_to_grid_displacement_curvature_fixed_grid(base_angle, curvature, inputstruct)
%             disp('Calculating grid displacement from base angle and curvature')
            deflection_angle = Whisker_Recording.base_angle_to_deflection_angle_curvature_fixed_grid(base_angle, curvature, inputstruct);
            grid_displacement = Whisker_Recording.base_angle_to_grid_displacement_deflection_angle_fixed_grid(base_angle, deflection_angle, inputstruct);
        end
        
        % from deflection angle and curvature
        function grid_displacement = deflection_angle_to_grid_displacement_curvature(deflection_angle, curvature, inputstruct)
%             disp('Calculating grid displacement from deflection angle and curvature')
            disp('Note that in a freely whisking case, a whisker can make movements that are impossible in a grid experiment')
            grid_displacement = sqrt((2*sin(deflection_angle)./curvature).^2-inputstruct.distance_grid^2);
            lw = 2*deflection_angle./curvature;
            f = figure;
            subplot(2,1,1)
            plot(lw , '.')
            title('Whisker length (mm)')
            ip = input('Constant whisker length ok, not violated? (y/n)', 's');
            if strcmp(ip, 'y')
            else
                error('Constant whisker length violated. Try transforming to whisker tip position ')
            end
            close(f)
        end
        
        % from base angle and grid distance (approximation)
        function grid_displacement = base_angle_to_grid_displacement_fixed_grid(base_angle, inputstruct)
%             disp('Estimating grid displacement from base angle')
            disp('Ignoring curvature')
            % Absolute base angles
            grid_displacement = inputstruct.distance_grid*tan(base_angle);
        end
  
        function grid_displacement = base_angle_to_grid_displacement_rel_nobaseline(base_angle, inputstruct)
%             disp('Estimating grid displacement from base angle')
            disp('Ignoring curvature, approximation without baseline (small base angle): Delta x = s Delta alpha')
            grid_displacement = inputstruct.distance_grid*base_angle;
        end
        
        
        %% To deflection angle
        % from base angle and whisker length (numerical approximation)
        function deflection_angle = base_angle_to_deflection_angle_length_whisker_fixed_grid(base_angle, inputstruct)
%             disp('Estimating deflection angle from base angle and whisker length')
            % s/lw = sin(D)cos(alpha+D)/D
            disp('The deflection angle cannot be solved analytically, doing a numerical approximation')
            ip = input('Assuming deflection angle values before (mostly negative) or after (mostly positive) maximum? (+/-)', 's');
            disp('Using a precision of 0.0001 rad')
            
            deflection_angle = nan*ones(size(base_angle));
            dd = 0.0001; % D precision
            
            % Make D values
            dn = -5:dd:5;
            k = inputstruct.distance_grid / inputstruct.length_whisker;
            for n = 1:length(base_angle)
                a = base_angle(n);
                sindcosadd = (sin(dn).*cos(a+dn))./dn;
                
                % find (global) maximum
                [~, maxsindcosaddi] = max(sindcosadd);
                if strcmp(ip, '+')
                    % use part after maximum
                    sindcosadd = sindcosadd(maxsindcosaddi:end);
                    dn = dn(maxsindcosaddi:end);
                    % until minimum
                    [~, minsindcosaddi] = min(sindcosadd);
                    sindcosadd = sindcosadd(1:minsindcosaddi);
                    dn = dn(1:minsindcosaddi);
                elseif strcmp(ip, '-')
                    % use part before maximum
                    sindcosadd = sindcosadd(1:maxsindcosaddi);
                    dn = dn(1:maxsindcosaddi);
                    % from minimum
                    [~, minsindcosaddi] = min(sindcosadd);
                    sindcosadd = sindcosadd(minsindcosaddi:end);
                    dn = dn(minsindcosaddi:end);
                end
                
                % Find D so that sindcosad = curvature(n)
                [~, mindi] = min(abs(sindcosadd-k));
                
                deflection_angle(n) = dn(mindi);
            end
        end
        
        % from curvature and whisker length
        function deflection_angle = curvature_to_deflection_angle_length_whisker(curvature, inputstruct)
%             disp('Calculating deflection angle from curvature and whisker length')
            % because this is linear, relative or absolute is the same
            deflection_angle = inputstruct.length_whisker*curvature/2;
        end

        function deflection_angle = curvature_to_deflection_angle_length_whisker_fixed_grid(curvature, inputstruct)
            deflection_angle = Whisker_Recording.curvature_to_deflection_angle_length_whisker(curvature, inputstruct);
        end
        
        % from grid displacement and length_whisker (numerical approximation)
        function deflection_angle = grid_displacement_to_deflection_angle_length_whisker_fixed_grid(grid_displacement, inputstruct)
%             disp('Estimating deflection angle from grid displacement and whisker length')
            % x = sqrt((lw sin(D)/D)^2-s^2)
            disp('The deflection angle cannot be solved analytically (sin(D)/D = const), doing a numerical approximation')
            disp('NB: Assuming deflection angle values between 0 and 4.5 rad (~1.4 pi)')
            disp('Using a precision of 0.0001 rad')
            
            deflection_angle = nan*ones(size(grid_displacement));
            dd = 0.0001; % D precision
            
            % Calculate values sin(D)/D
            dn = dd:dd:5;
            sindd = sin(dn)./dn;
            % Use only from max (at 0) to min, to prevent double solutions
            [~, minsinddi] = min(sindd);
            sindd = sindd(1:minsinddi);
            dn = dn(1:minsinddi);

            for n = 1:length(grid_displacement)
                K = sqrt(inputstruct.distance_grid^2 + grid_displacement(n)^2)./inputstruct.length_whisker;
                % Find D so that sin(D)/D = K
                [~, mindi] = min(abs(sindd-K));
                
                deflection_angle(n) = dn(mindi);
            end
        end
        
        % from whisker tip position and whisker length (numerical approximation)
        function deflection_angle = whisker_tip_position_to_deflection_angle_length_whisker(whisker_tip_position, inputstruct)
            % D/sinD = b/lw; alpha+D = atan(x/s)
            % s = wp(1), x = wp(2)
            disp('The deflection angle cannot be solved analytically, doing a numerical approximation')
            disp('NB: Assuming positive deflection angle values')
            disp('Using a precision of 0.0001 rad')
            
            
            [Ntime, Ndim] = size(whisker_tip_position);
            if ~(Ndim==2)
                error('Whisker tip position should be given as an Ntime x 2 array (1 = perp. to midline (s), 2 = parallel midline (x))')
            end
            deflection_angle = nan*ones(Ntime,1);
            
            dd = 0.0001; % D precision
            % Make D values
            dn = 0:dd:2*pi;
            % Calculate values D/sin(D)
            dsind = dn./sin(dn);
            % Use only to min at arount 3/2 pi, to prevent double solutions
            [~, dnpi] = min(abs(dn-pi));
            [~, dnmax2] = max(dsind(dnpi:end));
            dnmax2 = dnpi+dnmax2-1;
            dsind = dsind(1:dnmax2);
            dn = dn(1:dnmax2);
            
            b = sqrt(whisker_tip_position(:,1).^2+ whisker_tip_position(:,2).^2);
            
            for nt = 1:Ntime
                K = inputstruct.length_whisker./b(nt);
                [~, mindi] = min(abs(dsind-K));
                
                deflection_angle(nt) = dn(mindi);
            end
        end
        
        % from base angle and curvature (numerical approximation)
        function deflection_angle = base_angle_to_deflection_angle_curvature_fixed_grid(base_angle, curvature, inputstruct)
%             disp('Calculating deflection angle from base angle and curvature')
            % curvature = 2.*sin(deflection_angle).*cos(base_angle+deflection_angle)/distance_grid;
            % sin(D)cos(a+D) = (sin(a+2D)-sin(a))/2
            % So minima/maxima at D = -a/2-pi/4+pi n/2, 
            % where n=0 is a minimum and n=1 a maximum 
            disp('The deflection angle cannot be solved analytically, doing a numerical approximation')
%             disp('NB: Assuming positive deflection angle values')
            disp('Using a precision of 0.0001 rad')
            
            NC = length(curvature);
            NA = length(base_angle);

            if ~(NC == NA)
                error('Give same size arrays for base angles and curvatures')
            end
            
            deflection_angle = nan*ones(size(curvature));
            dd = 0.0001; % D precision
            
            % Make D values
            
            for n = 1:NA
                a = base_angle(n);
                dn = -a/2-pi/4:dd:-a/2-pi/4+pi/2; % Use only from min to max, to prevent double solutions
                
                sindcosad = sin(dn).*cos(a+dn);
                % Find D so that sindcosad = curvature(n)
                [~, mindi] = min(abs(sindcosad-inputstruct.distance_grid.*curvature(n)./2));
                
                deflection_angle(n) = dn(mindi);
            end
            disp('Note that in a freely whisking case, a whisker can make movements that are impossible in a grid experiment')
            lw = inputstruct.distance_grid*deflection_angle./(sin(deflection_angle).*cos(deflection_angle + base_angle));
            f = figure;
            subplot(2,1,1)
            plot(lw , '.')
            title('Whisker length (mm)')
            subplot(2,1,2)
            plot(base_angle)
            hold all
            plot(curvature)
            legend('base angle','curvature')
            ip = input('Constant whisker length ok, not systematically violated? (y/n)', 's');
            if strcmp(ip, 'y')
            else
                error('Constant whisker length violated. Try transforming to whisker tip position ')
            end
            close(f)
            
        end
        
        % from base angle and grid displacement  
        function deflection_angle = base_angle_to_deflection_angle_grid_displacement_fixed_grid(base_angle, grid_displacement, inputstruct)
%             disp('Calculating deflection angle from base angle and grid displacement length')
            deflection_angle = atan(grid_displacement/inputstruct.distance_grid)-base_angle;
        end
        
        % from curvature and grid displacement 
        function deflection_angle = curvature_to_deflection_angle_grid_displacement_fixed_grid(curvature, grid_displacement, inputstruct)
%             disp('Calculating deflection angle from curvature and grid displacement')
            deflection_angle = asin(curvature.*sqrt(grid_displacement.^2 + inputstruct.distance_grid^2)./2);
        end
        
                   
        %% To curvature
        % from base angle and whisker length (numerical approximation)
        function curvature = base_angle_to_curvature_length_whisker_fixed_grid(base_angle, inputstruct)
%             disp('Calculating curvature from base angle and whisker length')
            deflection_angle = Whisker_Recording.base_angle_to_deflection_angle_length_whisker_fixed_grid(base_angle, inputstruct);
            curvature = Whisker_Recording.deflection_angle_to_curvature_length_whisker(deflection_angle, inputstruct);
        end
        
        % from deflection angle and whisker length
        function curvature = deflection_angle_to_curvature_length_whisker(deflection_angle, inputstruct)
%             disp('Calculating curvature from deflection angle and whisker length')
            curvature = 2*deflection_angle/inputstruct.length_whisker;
        end
        
        function curvature = deflection_angle_to_curvature_length_whisker_fixed_grid(deflection_angle, inputstruct)
%             disp('Calculating curvature from deflection angle and whisker length')
            curvature = 2*deflection_angle/inputstruct.length_whisker;
        end
        
        % from displacement and whisker length (numerical approximation)
        function curvature = grid_displacement_to_curvature_length_whisker_fixed_grid(grid_displacement, inputstruct)
%             disp('Calculating curvature from grid displacement and whisker length')
            deflection_angle = Whisker_Recording.grid_displacement_to_deflection_angle_length_whisker_fixed_grid(grid_displacement, inputstruct);
            curvature = Whisker_Recording.deflection_angle_to_curvature_length_whisker(deflection_angle, inputstruct);
        end
        
        % from whisker tip position and whisker length (numerical approximation)
        function curvature = whisker_tip_position_to_curvature_length_whisker(whisker_tip_position, inputstruct)
            deflection_angle = Whisker_Recording.whisker_tip_position_to_deflection_angle_length_whisker(whisker_tip_position, inputstruct);
            curvature = Whisker_Recording.deflection_angle_to_curvature_length_whisker(deflection_angle, inputstruct);
        end
        
        % from base angle and deflection angle
        function curvature = base_angle_to_curvature_deflection_angle_fixed_grid(base_angle, deflection_angle, inputstruct)
%             disp('Calculating curvature from base angle and deflection angle')
            curvature = 2.*sin(deflection_angle).*cos(base_angle+deflection_angle)/inputstruct.distance_grid;
            disp('Note that in a freely whisking case, a whisker can make movements that are impossible in a grid experiment')
            lw = 2*deflection_angle./curvature;
            f = figure;
            subplot(2,1,1)
            plot(lw , '.')
            title('Whisker length (mm)')
            subplot(2,1,2)
            plot(base_angle)
            hold all
            plot(deflection_angle)
            legend('base angle','deflection angle')
            ip = input('Constant whisker length ok, not violated? (y/n)', 's');
            if strcmp(ip, 'y')
            else
                error('Constant whisker length violated. Try transforming to whisker tip position ')
            end
            close(f)
        end
        
        % from base angle and grid displacement 
        function curvature = base_angle_to_curvature_grid_displacement_fixed_grid(base_angle, grid_displacement, inputstruct)
%             disp('Calculating curvature from base angle and grid displacement')
            deflection_angle = Whisker_Recording.base_angle_to_deflection_angle_grid_displacement_fixed_grid(base_angle, grid_displacement, inputstruct);
            curvature = Whisker_Recording.base_angle_to_curvature_deflection_angle_fixed_grid(base_angle, deflection_angle, inputstruct);
        end
        
        % from deflection angle and grid displacement
        function curvature = deflection_angle_to_curvature_grid_displacement_fixed_grid(deflection_angle, grid_displacement, inputstruct)
%             disp('Calculating curvature from deflection angle and grid displacement')
            curvature = 2.*sin(deflection_angle)./sqrt(grid_displacement.^2+inputstruct.distance_grid.^2);
        end        
        
        %% To base angle
        % deflection angle and whisker length
        function base_angle = deflection_angle_to_base_angle_length_whisker_fixed_grid(deflection_angle, inputstruct)
%             disp('Calculating base angle from deflection angle and whisker length')
            grid_displacement = Whisker_Recording.deflection_angle_to_grid_displacement_length_whisker_fixed_grid(deflection_angle, inputstruct);
            base_angle = Whisker_Recording.deflection_angle_to_base_angle_grid_displacement_fixed_grid(deflection_angle, grid_displacement, inputstruct);
        end
        
        % curvature and whisker length
        function base_angle = curvature_to_base_angle_length_whisker_fixed_grid(curvature, inputstruct)
%             disp('Calculating base angle from curvature and whisker length')
            deflection_angle = Whisker_Recording.curvature_to_deflection_angle_length_whisker(curvature, inputstruct);
            grid_displacement = Whisker_Recording.curvature_to_grid_displacement_length_whisker_fixed_grid(curvature, inputstruct);
            base_angle = Whisker_Recording.deflection_angle_to_base_angle_grid_displacement_fixed_grid(deflection_angle, grid_displacement, inputstruct);
        end
        
        % grid displacement and whisker length (numerical approximation)
        function base_angle = grid_displacement_to_base_angle_length_whisker_fixed_grid(grid_displacement, inputstruct)
%             disp('Estimating base angle from grid displacement and whisker length')
            deflection_angle = Whisker_Recording.grid_displacement_to_deflection_angle_length_whisker_fixed_grid(grid_displacement, inputstruct); % numerical approximation
            base_angle = Whisker_Recording.deflection_angle_to_base_angle_grid_displacement_fixed_grid(deflection_angle, grid_displacement, inputstruct);
        end
         
        % from whisker tip position and whisker length (numerical approximation)
        function base_angle = whisker_tip_position_to_base_angle_length_whisker(whisker_tip_position, inputstruct)
            deflection_angle = Whisker_Recording.whisker_tip_position_to_deflection_angle_length_whisker(whisker_tip_position, inputstruct);
            apD = atan(whisker_tip_position(:,2)./whisker_tip_position(:,1));
            base_angle = apD - deflection_angle;
        end
        
        % from whisker tip position and deflection angle (numerical approximation)
        function base_angle = deflection_angle_to_base_angle_whisker_tip_position(deflection_angle, whisker_tip_position, ~)
            apD = atan(whisker_tip_position(:,2)./whisker_tip_position(:,1));
            base_angle = apD - deflection_angle;
        end
        
        % deflection angle and curvature
        function base_angle = deflection_angle_to_base_angle_curvature_fixed_grid(deflection_angle, curvature, inputstruct)
%             disp('Calculating base angle from deflection angle and curvature')
            NC = length(curvature);
            ND = length(deflection_angle);

            if ~(NC == ND)
                error('Give same size arrays for deflection angles and curvatures')
            end
            
            if ~isfield(inputstruct, 'length_whisker')
                inputstruct.length_whisker = Whisker_Recording.deflection_angle_to_length_whisker_curvature_fixed_grid(deflection_angle, curvature, inputstruct);
            else
                aw = input('This calculation is not needed when the whisker length is known. Are you sure this is correct? (y/n)', 's');
                if strcmp(aw, 'n')
                    error('Ending calculation')
                elseif strcmp(aw, 'y')
                    disp('Continuing calculation')
                else
                    disp('Please answer y or n')
                end
            end
            grid_displacement = Whisker_Recording.deflection_angle_to_grid_displacement_length_whisker_fixed_grid(deflection_angle, inputstruct);
            base_angle = Whisker_Recording.deflection_angle_to_base_angle_grid_displacement_fixed_grid(deflection_angle, grid_displacement, inputstruct);
            
        end 
         
        % deflection angle and grid displacement
        function base_angle = deflection_angle_to_base_angle_grid_displacement_fixed_grid(deflection_angle, grid_displacement, inputstruct)
%             disp('Calculating base angle from deflection angle and grid displacement')
            base_angle = atan(grid_displacement/inputstruct.distance_grid)-deflection_angle;
        end
        
        % curvature and grid displacment
        function base_angle = curvature_to_base_angle_grid_displacement_fixed_grid(curvature, grid_displacement, inputstruct)
%             disp('Calculating base angle from curvature and grid displacement')
            deflection_angle = Whisker_Recording.curvature_to_deflection_angle_grid_displacement_fixed_grid(curvature, grid_displacement, inputstruct);
            base_angle = Whisker_Recording.deflection_angle_to_base_angle_grid_displacement_fixed_grid(deflection_angle, grid_displacement, inputstruct);
        end
        
        % from grid displacement and grid distance (approximation)
        function base_angle = grid_displacement_to_base_angle_fixed_grid(grid_displacement, inputstruct)
%             disp('Approximating base angle from grid displacement')
            disp('Ignoring curvature')
            % Absolute base angles
            base_angle = atan(grid_displacement/inputstruct.distance_grid);
        end
  
        function base_angle = grid_displacement_to_base_angle_rel_nobaseline(grid_displacement, inputstruct)
%             disp('Approximating base angle from grid displacement')
            disp('Ignoring curvature, approximation without baseline (small base angle): Delta alpha = Delta x/s ')
            base_angle = grid_displacement./inputstruct.distance_grid;
        end     
        
        %% To whisker length
        
        % from base angle and deflection angle
        function length_whisker = base_angle_to_length_whisker_deflection_angle_fixed_grid(base_angle, deflection_angle, inputstruct)
%             disp('Calculating whisker length from base angle and deflection angle')
            ND = length(deflection_angle);
            Na = length(base_angle);
            
            if ~(ND == Na)
                error('Give same size arrays for base and deflection angles')
            end
            if ND == 1
                length_whisker = inputstruct.distance_grid*deflection_angle./(sin(deflection_angle).*cos(deflection_angle + base_angle));
            elseif ND>1
%                 disp('Making a fit for the whisker length')
                disp('Calculating the mean of the measurements')
                % Transpose if needed
                [NDx, ~] = size(deflection_angle);
                [Nax, ~] = size(base_angle);
                if NDx == 1
                    deflection_angle = deflection_angle';
                end
                if Nax == 1
                    base_angle = base_angle';
                end
                disp('Note that in a freely whisking case, a whisker can make movements that are impossible in a grid experiment')
                lw = inputstruct.distance_grid*deflection_angle./(sin(deflection_angle).*cos(deflection_angle + base_angle));
                f = figure;
                subplot(2,1,1)
                plot(lw , '.')
                title('Whisker length (mm)')
                subplot(2,1,2)
                plot(base_angle)
                hold all
                plot(deflection_angle)
                legend('base angle','deflection_angle')
                ip = input('Constant whisker length ok, not violated? (y/n)', 's');
                if strcmp(ip, 'y')
                else
                    error('Constant whisker length violated. Try transforming to whisker tip position ')
                end
                close(f)
                
                length_whisker = nanmean(lw);
            end
        end
        
        % from base angle and curvature
        function length_whisker = base_angle_to_length_whisker_curvature_fixed_grid(base_angle, curvature, inputstruct)
%             disp('Calculating whisker length from base angle and curvature')
            deflection_angle = Whisker_Recording.base_angle_to_deflection_angle_curvature_fixed_grid(base_angle, curvature, inputstruct);
            length_whisker = Whisker_Recording.deflection_angle_to_length_whisker_curvature(deflection_angle, curvature, inputstruct);
        end
        
        % from grid displacement and base angle
        function length_whisker = base_angle_to_length_whisker_grid_displacement_fixed_grid(base_angle, grid_displacement, inputstruct)
%             disp('Calculating whisker length from base angle and grid displacement')
            deflection_angle = Whisker_Recording.base_angle_to_deflection_angle_grid_displacement_fixed_grid(base_angle, grid_displacement, inputstruct);
            length_whisker = Whisker_Recording.deflection_angle_to_length_whisker_grid_displacement_fixed_grid(deflection_angle, grid_displacement, inputstruct);
        end
        
        % from deflection angle and curvature (mean)
        function length_whisker = deflection_angle_to_length_whisker_curvature_fixed_grid(deflection_angle, curvature, ~)
            disp('NB: Fixed grid assumption not needed.')
            length_whisker = Whisker_Recording.deflection_angle_to_length_whisker_curvature(deflection_angle, curvature);
        end
        
        % from deflection angle and curvature (mean)
        function length_whisker = deflection_angle_to_length_whisker_curvature(deflection_angle, curvature, ~)
%             disp('Calculating whisker length from deflection angle and curvature')
            % c = D*(lw/2)
            % NB Note that due to deformations that are stronger for large
            % |D|, the deflection angle is 'too big' at large |D|. Taking
            % the mean is better at compensating for this than fitting,
            % which is more sensitive to the large |D| values.
            ND = length(deflection_angle);
            NC = length(curvature);

            if ~(ND == NC)
                error('Give same size arrays for deflection angles and curvatures')
            end
            if ND == 1
                length_whisker = 2*deflection_angle./curvature;
            elseif ND>1
%                 disp('Making a fit for the whisker length')
                disp('Calculating the mean of the measurements')
                % Transpose if needed
                [NDx, ~] = size(deflection_angle);
                [NCx, ~] = size(curvature);
                if NDx == 1
                    deflection_angle = deflection_angle';
                end
                if NCx == 1
                    curvature = curvature';
                end
                % Make a linear fit
%                 f = fit(curvature, deflection_angle, 'poly1');
%                 length_whisker1 = 2*f.p1;
                
                % Take the mean
                div = deflection_angle./curvature;
                length_whisker = 2*nanmean(div(isfinite(div)));
                
%                 length_whisker2 = 2*nanmean(div(isfinite(div)));
                
%                 figure
%                 subplot(1,2,1)
%                 hold all
% %                 plot(curvature)
%                 plot(deflection_angle)
%                 plot(length_whisker1*curvature/2)
%                 plot(length_whisker2*curvature/2)
%                 subplot(1,2,2)
%                 hold all
%                 plot(curvature, deflection_angle, '.')
%                 plot(curvature, length_whisker1*curvature/2, '.')
%                 plot(curvature, length_whisker2*curvature/2, '.')
%                 keyboard

            else
                error('Recordings are empty')
            end
        end
                                        
        % from grid displacement and deflection angle (mean)
        function length_whisker = deflection_angle_to_length_whisker_grid_displacement_fixed_grid(deflection_angle, grid_displacement, inputstruct)
%             disp('Calculating whisker length from deflection angle and grid displacement')
            % Using x = sqrt((lw.*sin(D)./D)^2 - s^2);
            Nx = length(grid_displacement);
            ND = length(deflection_angle);

            if ~(Nx == ND)
                error('Give same size arrays for deflection angles and grid displacement')
            end
            lw = deflection_angle.*sqrt(grid_displacement.^2 + inputstruct.distance_grid^2)./ sin(deflection_angle);
            if ND == 1
                length_whisker = lw;
            elseif ND>1
                disp('Using the mean over the recording values')
                length_whisker = nanmean(lw);
            else
                error('Recordings are empty')
            end
        end
        
        % from grid displacement and curvature
        function length_whisker = curvature_to_length_whisker_grid_displacement_fixed_grid(curvature, grid_displacement, inputstruct)
%             disp('Calculating whisker length from curvature and grid displacement')
            deflection_angle = Whisker_Recording.curvature_to_deflection_angle_grid_displacement_fixed_grid(curvature, grid_displacement, inputstruct);
            length_whisker = Whisker_Recording.deflection_angle_to_length_whisker_grid_displacement_fixed_grid(deflection_angle, grid_displacement, inputstruct);
        end
        
        function length_whisker = whisker_tip_position_to_length_whisker(whisker_tip_position, ~)
            length_whisker = NaN;
            error('This is not possible!')
            % Not possible!
%             s = whisker_tip_position(:,1);
%             x = whisker_tip_position(:,2);
%             Ns = length(s);
%             
%             b = sqrt(x.^2 + s.^2);
%             apD = atan(x./s);
%             
%                      
%             if Ns == 1
%                 length_whisker = 
%             elseif ND>1
%                 len
%             end
        end
        
        %% To whisker tip position
        function whisker_tip_position = base_angle_to_whisker_tip_position_deflection_angle_length_whisker(base_angle, deflection_angle, inputstruct)
            Ntime_a = length(base_angle);
            Ntime_D = length(deflection_angle);
            if ~(Ntime_a == Ntime_D)
                error('Recordings for base angle and deflection angle should have the same length')
            end
            whisker_tip_position = nan*ones(Ntime_a, 2);
            b  = inputstruct.length_whisker.*sin(deflection_angle)./deflection_angle;
            b(find(deflection_angle==0)) = inputstruct.length_whisker;

            whisker_tip_position(:,1) = b.*cos(base_angle + deflection_angle); % s
            whisker_tip_position(:,2) = b.*sin(base_angle + deflection_angle); % x
            
        end
        
        function whisker_tip_position = base_angle_to_whisker_tip_position_curvature_length_whisker(base_angle, curvature, inputstruct)
            Ntime_a = length(base_angle);
            Ntime_c = length(curvature);
            if ~(Ntime_a == Ntime_c)
                error('Recordings for base angle and curvature should have the same length')
            end
            whisker_tip_position = nan*ones(Ntime_a, 2);
            deflection_angle = inputstruct.length_whisker.*curvature./2;
            b  = 2.*sin(deflection_angle)./curvature;
            b(find(deflection_angle==0)) = inputstruct.length_whisker;
            whisker_tip_position(:,1) = b.*cos(base_angle + deflection_angle); % s
            whisker_tip_position(:,2) = b.*sin(base_angle + deflection_angle); % x
        end
        
        
        %% Degrees and radians
        function radian = degree_to_radian(degree)
            radian = pi.*degree./180;
        end

        function degree = radian_to_degree(radian)
            degree = 180.*radian./pi;
        end
        
        %% Unit names
        function unitnamevec = set_unitname(fieldnamevec)
           fieldnamevec = cellstr(fieldnamevec);
           Nname = length(fieldnamevec);
           
           unitnamevec = cell(size(fieldnamevec));
           nout = 0;
           for n=1:Nname
               nout = nout+1;
               fieldname = fieldnamevec{n};
               if strcmp(fieldname, 'curvature')
                   unitnamevec{nout} = 'mm-1';
               elseif strcmp(fieldname, 'grid_displacement')
                   unitnamevec{nout} = 'mm';
               elseif strcmp(fieldname, 'deflection_angle')
                   unitnamevec{nout} = 'radian';
               elseif strcmp(fieldname, 'base_angle')
                   unitnamevec{nout} = 'radian';
               elseif strcmp(fieldname, 'whisker_tip_position')
%                    unitnamevec{nout} = 'mm';
%                    unitnamevec{nout+1} = 'mm';
%                    nout=nout+1;
                    unitnamevec{nout} = 'mm';
               end
           end
           if nout == 1  
               unitnamevec = unitnamevec{1};
           end
           unitnamevec = cellstr(unitnamevec);
        end
        
    end
    
end