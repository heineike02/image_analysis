function images = get_image_fnames(fname_conv,imdir,channel_to_image,pos_fnames_nn,Nchan)

if strcmp(fname_conv,'JSO')
    images = struct();
    files=dir(char(strcat(imdir,channel_to_image,'-wheel_',pos_fnames_nn,'_t*.tiff'))); 
    nTimes = size(files,1);

    %Get initial number
    init_fname = files(1).name;
    init_ind = strfind(init_fname,'_t');
    dot_ind = strfind(init_fname,'.');
    init_num = str2double(init_fname(init_ind+2:(dot_ind-1)));
    %picks out the first number and checks what it is modulo the number of channels 
    channel_mod = mod(init_num,Nchan)+Nchan*(mod(init_num,Nchan)==0) ;
    %This is a pain because the numbering is not 01,02,03 and
    %messes up the ordering of the files. 
    images_cell = cell(1,nTimes);
    for mm = 1:nTimes  
       %Nchan should be an input automatically calculated -
       %analysis_parameters.channels should be a cell with all the channels
       %you want. 
       %Nchan = 2; % if only want the results of 1 color, but imaged in 2 colors
       fname = char(strcat(channel_to_image,'_',pos_fnames_nn,'_t',int2str(channel_mod+(mm-1)*Nchan),'.tiff')); 
       images_cell{mm} = fname;
    end
    images.(channel_to_image) = images_cell;
elseif strcmp(fname_conv,'Micromanager')
    images = struct();
    files=dir(strcat(imdir,pos_fnames_nn,filesep,'*',channel_to_image,'*'));
    %images should be in the correct time order already
    nTimes = size(files,1);
    images_cell = cell(1,nTimes);
    for mm = 1:nTimes
        images_cell{mm} = strcat(pos_fnames_nn,filesep,files(mm).name);
    end
    images.(channel_to_image) = images_cell;
elseif strcmp(fname_conv,'HCS_Nikon')
elseif strcmp(fname_conv,'Metamorph')
    images = struct();
    files = dir(strcat(imdir,pos_fnames_nn,'_*.TIF'));
    nTimes = size(files,1);
    
    %Get initial number
    init_fname = files(1).name;
    init_ind = strfind(init_fname,'_t');
    dot_ind = strfind(init_fname,'.');
    init_num = str2double(init_fname(init_ind+2:(dot_ind-1)));
    %picks out the first number and checks what it is modulo the number of channels 
    channel_mod = mod(init_num,Nchan)+Nchan*(mod(init_num,Nchan)==0) ;
    %This is a pain because the numbering is not 01,02,03 and
    %messes up the ordering of the files. 
    images_cell = cell(1,nTimes);
    for mm = 1:nTimes 
       fname = char(strcat(pos_fnames_nn,'_t',int2str(channel_mod+(mm-1)*Nchan),'.TIF'));
       images_cell{mm} = fname;
    end
    images.(channel_to_image) = images_cell;
    
end


