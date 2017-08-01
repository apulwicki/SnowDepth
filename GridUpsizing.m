function [ SWEoutput ] = GridUpsizing( gridSize )
%Combines all measurement locations within 'gridSize' combo of cells for
%each glacier and outputs the mean and std

global options
load TopoSWE.mat topo_sampled
load Full.mat fullSWE

for g = 1:3; glacier = options.glacier{g};
    SWEinput.(glacier) = fullSWE.S2.input.(glacier);
end
SWEinput = ObsInCell(SWEinput,topo_sampled);

 % Set up matrix of cell number
 for g = 1:3;   glacier = options.glacier{g};
t = options.mapsize(g,1)*options.mapsize(g,2);
num = 1:t; num = reshape(num,options.mapsize(g,2),options.mapsize(g,1))';    numO = num;

sizeG = gridSize;
    C = sqrt(sizeG);
    
 % Add nan rows to make the size a multiple of the new grid size
    if floor(size(num,1)/C) ~= size(num,1)/C
        numrowtoadd = ceil(size(num,1)/C)*C-size(num,1);
        num = [num; nan(numrowtoadd,size(num,2))];
    end
    if floor(size(num,2)/C) ~= size(num,2)/C
        numrowtoadd = ceil(size(num,2)/C)*C-size(num,2);
        num = [num nan(size(num,1),numrowtoadd)];
    end

% Change all grid cell numbers to the middle or one over within the size window
for i = 1:C:size(num,1)
   for j =  1:C:size(num,2)
       if mod(C,2) == 0    
                T = num(i+C/2,j+C/2);
       else;    T = num(i+(C-1)/2, j+(C-1)/2);
       end
       num(i:i+C-1,j:j+C-1) = T;
   end
end
    %Change back to original size
num = num(1:size(numO,1), 1:size(numO,2));
AllCellN.(glacier) = num;
 end
 
% Find sampled cells
sampledP = SampledCell(AllCellN);
for g = 1:3; glacier = options.glacier{g};
    SWEinput.(glacier)(:,4) = sampledP.(glacier);
end

[SWEoutput, ~, stdGrid] = ObsInCell(SWEinput, topo_sampled);
for g = 1:3; glacier = options.glacier{g};
    SWEoutput.(glacier)(:,4) = stdGrid.(glacier);
end

end

