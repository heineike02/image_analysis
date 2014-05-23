function images = get_image_fnames(fname_conv,imdir,channel_to_image,pos_fnames_nn,Nchan)

if strcmp(fname_conv,'JSO')
    images = struct();
    files=dir(char(strcat(imdir,channel_to_image,'_',pos_fnames_nn,'_t*.tiff'))); 
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
    for mm = 1:nTimes 
       fname = char(strcat(channel_to_image,'_',pos_fnames_nn,'_t',int2str(channel_mod+(mm-1)*Nchan),'.tiff'));
       images(mm).(channel_to_image) = fname;
    end
elseif strcmp(fname_conv,'Micromanager')
    images = struct();
    files=dir(strcat(imdir,pos_fnames_nn,'\*',channel_to_image,'*'));
    %images should be in the correct time order already
    nTimes = size(files,1);
    for mm = 1:nTimes
        images(mm).(channel_to_image) = strcat(pos_fnames_nn,'\',files(mm).name);
    end
elseif strcmp(fname_conv,'HCS_Nikon')
end


