function [cellstruct] = Neighborhood_SYC(Cxloc, Cyloc, tim)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-Cyloc and Cxloc are information about the xy locations
%-tim are the frame indices
%-cellx and celly are mxn matrices where m = number of frames and n = cell number
%-produces a cellstruct = indexed the same way as cellx and celly but
%contains superfield:
    %cells = each individual cell in a particular frame
%contains fields: 
    %x = vector of x coords of neighboring cells 
    %y = vector of y coords of neighboring cells
    %len = number of neighboring cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize the megastruct
cellstruct = struct();
%initialize the smaller struct
smallstruct = struct('x',[], 'y', [], 'len', []);

frames = max(tim);
%loop through each frame
for j = 1:frames
    cellx = Cxloc(find(tim==j));
    celly = Cyloc(find(tim==j));
%loop through each cell
for i = 1:length(cellx) 
    diffx = cellx- cellx(i);
    diffy = celly -celly(i);
    threshold = 25; %20pixels is upper limit radius within which cells have neighbor
    indx = find(abs(diffx) <= threshold); %indices x %%so instead of doing this, threshold on x and then threshold on y, then apply distance formula to remaining
    indy = find(abs(diffy) <= threshold); %indices y %%also test how the indices of x and y line up??, display on an image
    %find the intersection between indx and indy
    indxy = intersect(indx, indy);
    near_x = cellx(indxy);
    near_y = celly(indxy);
    
    %test code
    num_neighx = length(near_x);
    num_neighy = length(near_y);
    %if num_neighx ~= num_neighy
    %   display('what what??')
    %end
    
    smallstruct(i).x = near_x;
    smallstruct(i).y =  near_y;
    smallstruct(i).len = (num_neighx-1);
end
    cellstruct(j).cells = smallstruct;

end

end