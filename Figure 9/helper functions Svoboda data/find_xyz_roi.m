function M = find_xyz_roi( s, volume )
% find_xyz_roi finds the ids and x, y and z-coordinates of each roi in a
% volume
% INPUT: Svoboda data structure s and a volume number (>=2)
% OUTPUT: matrix M (4, # cells in volume), where first row is the cell id,
% and the second, third and fourth the x, y and z-coordinate 

[ncellsvolume, ntime] = size(s.timeSeriesArrayHash.value{volume}.valueMatrix);
%% check number of cells
ncellsok = 1;
nctot = 0;
nplanes = length(s.timeSeriesArrayHash.value{volume}.imagingPlane);
for np = 1:nplanes
    nctot = nctot+length(s.timeSeriesArrayHash.value{volume}.imagingPlane{np}.ids);
end
if ~(nctot == ncellsvolume)
    disp('Warning: Number of recorded cells not equal to sum of cells in planes')
    ncellsok = 0;
end

M = nan*ones(4,ncellsvolume);

nn=0;
for np = 1:nplanes
    %% find cell ids for this plane
    cellidsthisplane = s.timeSeriesArrayHash.value{volume}.imagingPlane{np}.ids;
    ncellsthisplane = length(cellidsthisplane);
    if ncellsok == 1
        % Extra check: are cell ids from different locations (value matrix and imaging planes) the same?
        if ~ (sum(s.timeSeriesArrayHash.descrHash{volume}.value{1}(np).roiIds == cellidsthisplane) == ncellsthisplane)
            disp('Something went wrong, cell ids not correct')
            keyboard
        end
    else
        % Since # cells in value matrix is not the same as # cells in
        % imaging planes, check for each cell id in imaging planes whether it was recorded
        cellokvec = zeros(1,ncellsthisplane);
        for nc = 1:ncellsthisplane
            if ~isempty(find(s.timeSeriesArrayHash.value{volume}.ids == cellidsthisplane(nc)))
                cellokvec(nc) = 1;
            end
        end
        ncellsthisplane = sum(cellokvec);
        if ncellsthisplane>0
            cellidsthisplane = cellidsthisplane(cellokvec);
        else
            disp(['Plane ' num2str(np) ' has no recorded cells'])
            cellidsthisplane = [];
        end
        
    end
    M(1, nn+1:nn+ncellsthisplane) = cellidsthisplane;
    
    
    %% z-coordinate (depth): 
    % planes are 15 micrometers separated
    % each volume contains 3 planes
    % so volume 2, plane 1 correspondse to '0', plane 2 to '15', etc
    M(4, nn+1:nn+ncellsthisplane) = ((np-1)+(volume-2)*3)*15;
    
    %% find x and y coordinate: take average over pixel locations
    for nc = 1:ncellsthisplane
        M(2:3, nn+nc) =  mean(s.timeSeriesArrayHash.descrHash{volume}.value{1}(np).rois(nc).cornersXY,2);
    end

    nn = nn+ncellsthisplane;
end

end

