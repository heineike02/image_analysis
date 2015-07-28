function varargout = image_processing_GUI(varargin)
% RUN_IMAGE_PROCESSING MATLAB code for run_image_processing.fig
%      RUN_IMAGE_PROCESSING, by itself, creates a new RUN_IMAGE_PROCESSING or raises the existing
%      singleton*.
%
%      H = RUN_IMAGE_PROCESSING returns the handle to a new RUN_IMAGE_PROCESSING or the handle to
%      the existing singleton*.
%
%      RUN_IMAGE_PROCESSING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RUN_IMAGE_PROCESSING.M with the given input arguments.
%
%      RUN_IMAGE_PROCESSING('Property','Value',...) creates a new RUN_IMAGE_PROCESSING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before image_processing_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to image_processing_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help run_image_processing

% Last Modified by GUIDE v2.5 26-Jul-2015 22:48:59

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @image_processing_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @image_processing_GUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


%To Do: 
%
% Display species magnification params
% alter params

% --- Executes just before run_image_processing is made visible.
function image_processing_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to run_image_processing (see VARARGIN)

% Choose default command line output for run_image_processing
handles.output = hObject;

% Initialize various parameters as the listed defaults
handles.ipdir = get(handles.ipdir,'String');
handles.base_dir = get(handles.base_dir,'String');
fname_conv_options = get(handles.fname_conv,'String');
handles.fname_conv = fname_conv_options{1};

%handles.phases = get(handles.phases,'String');
%handles.shift_spacing = str2num(get(handles.shift_spacing,'String'));

% handles.thresh = str2double(get(handles.thresh,'String'));
% handles.storeim = get(handles.storeim,'Value');
% 
% handles.coarse_smooth = str2double(get(handles.coarse_smooth,'String'));
% handles.local_smooth = str2double(get(handles.local_smooth,'String'));
% handles.maxcells = str2double(get(handles.maxcells,'String'));
% handles.deconvlucy_iterations = str2double(get(handles.deconvlucy_iterations,'String'));
% handles.close_max = str2double(get(handles.close_max,'String'));
% handles.far_max = str2double(get(handles.far_max,'String'));
% handles.ne_pixels = str2double(get(handles.ne_pixels,'String'));
% handles.edge_margin = str2double(get(handles.edge_margin,'String'));
% handles.track_memory = str2double(get(handles.track_memory,'String'));
% handles.min_points_for_traj = str2double(get(handles.min_points_for_traj,'String'));


%For species_magnification_params function
op_amp_options = get(handles.op_amp,'String');
handles.op_amp = op_amp_options{1};
species_options = get(handles.species,'String');
handles.species = species_options{1};
handles.maxdisp_1x = str2double(get(handles.maxdisp_1x,'String'));


% Update handles structure
guidata(hObject, handles);



% UIWAIT makes run_image_processing wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = image_processing_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function ipdir_Callback(hObject, eventdata, handles)
% hObject    handle to ipdir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.ipdir = get(hObject,'String');

% Hints: get(hObject,'String') returns contents of ipdir as text
%        str2double(get(hObject,'String')) returns contents of ipdir as a double
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function ipdir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ipdir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function thresh_Callback(hObject, eventdata, handles)
% hObject    handle to thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.thresh = str2double(get(hObject,'String'));

% Hints: get(hObject,'String') returns contents of thresh as text
%        str2double(get(hObject,'String')) returns contents of thresh as a double
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function thresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function phases_Callback(hObject, eventdata, handles)
% hObject    handle to phases (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of phases as text
%        str2double(get(hObject,'String')) returns contents of phases as a double
% --- Executes during object creation, after setting all properties.

handles.phases = regexp(get(hObject,'String'), ',', 'split');

guidata(hObject, handles)

function phases_CreateFcn(hObject, eventdata, handles)
% hObject    handle to phases (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in fname_conv.
function fname_conv_Callback(hObject, eventdata, handles)
% hObject    handle to fname_conv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns fname_conv contents as cell array
%        contents{get (hObject,'Value')} returns selected item from fname_conv
items = get(hObject,'String');
index_selected = get(hObject,'Value');
item_selected = items{index_selected};
handles.fname_conv = item_selected;
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function fname_conv_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fname_conv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'String',{'Micromanager';'JSO';'HCS_Nikon';'Metamorph'});


% --- Executes on button press in storeim.
function storeim_Callback(hObject, eventdata, handles)
% hObject    handle to storeim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of storeim

if (get(hObject,'Value') == get(hObject,'Max'))
	handles.storeim = 1;
else
	handles.storeim = 0;
end
guidata(hObject, handles)


% --- Executes on button press in species_magnification_params.
function species_magnification_params_Callback(hObject, eventdata, handles)
% hObject    handle to species_magnification_params (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[circ, siz, rad, maxdisp, std_thresh] = species_magnification_params(handles.species, handles.op_amp, handles.ipdir, handles.maxdisp_1x);

%handles.circ_file = circ_file;
set(handles.circ_file,'String',['circ',handles.species,'_',handles.op_amp,'.mat'])
set(handles.sizx,'String',siz(1))
set(handles.sizy,'String',siz(2))
set(handles.rad,'String',rad)
set(handles.maxdisp,'String',maxdisp)
set(handles.std_thresh,'String',std_thresh)

set(handles.load_species_magnification_params,'Visible','on')

guidata(hObject,handles)

% --- Executes on selection change in op_amp.
function op_amp_Callback(hObject, eventdata, handles)
% hObject    handle to op_amp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns op_amp contents as cell array
%        contents{get(hObject,'Value')} returns selected item from op_amp
items = get(hObject,'String');
index_selected = get(hObject,'Value');
item_selected = items{index_selected};
handles.op_amp = item_selected;
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function op_amp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to op_amp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'String',{'1x';'1p5x'});


% --- Executes on selection change in species.
function species_Callback(hObject, eventdata, handles)
% hObject    handle to species (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns species contents as cell array
%        contents{get(hObject,'Value')} returns selected item from species
items = get(hObject,'String');
index_selected = get(hObject,'Value');
item_selected = items{index_selected};
handles.species = item_selected;
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function species_CreateFcn(hObject, eventdata, handles)
% hObject    handle to species (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'String',{'SC';'KL'});



function maxdisp_1x_Callback(hObject, eventdata, handles)
% hObject    handle to maxdisp_1x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxdisp_1x as text
%        str2double(get(hObject,'String')) returns contents of maxdisp_1x as a double


% --- Executes during object creation, after setting all properties.
function maxdisp_1x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxdisp_1x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function coarse_smooth_Callback(hObject, eventdata, handles)
% hObject    handle to coarse_smooth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of coarse_smooth as text
%        str2double(get(hObject,'String')) returns contents of coarse_smooth as a double
handles.coarse_smooth = str2double(get(hObject,'String'));
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function coarse_smooth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to coarse_smooth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function local_smooth_Callback(hObject, eventdata, handles)
% hObject    handle to local_smooth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of local_smooth as text
%        str2double(get(hObject,'String')) returns contents of local_smooth as a double
handles.local_smooth = str2double(get(hObject,'String'));
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function local_smooth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to local_smooth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxcells_Callback(hObject, eventdata, handles)
% hObject    handle to maxcells (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxcells as text
%        str2double(get(hObject,'String')) returns contents of maxcells as a double
handles.maxcells = str2double(get(hObject,'String'));
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function maxcells_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxcells (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function deconvlucy_iterations_Callback(hObject, eventdata, handles)
% hObject    handle to deconvlucy_iterations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of deconvlucy_iterations as text
%        str2double(get(hObject,'String')) returns contents of deconvlucy_iterations as a double
handles.deconvlucy_iterations = str2double(get(hObject,'String'));
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function deconvlucy_iterations_CreateFcn(hObject, eventdata, handles)
% hObject    handle to deconvlucy_iterations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function far_max_Callback(hObject, eventdata, handles)
% hObject    handle to far_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of far_max as text
%        str2double(get(hObject,'String')) returns contents of far_max as a double
handles.far_max = str2double(get(hObject,'String'));
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function far_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to far_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function close_max_Callback(hObject, eventdata, handles)
% hObject    handle to close_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of close_max as text
%        str2double(get(hObject,'String')) returns contents of close_max as a double
handles.close_max = str2double(get(hObject,'String'));
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function close_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to close_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ne_pixels_Callback(hObject, eventdata, handles)
% hObject    handle to ne_pixels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ne_pixels as text
%        str2double(get(hObject,'String')) returns contents of ne_pixels as a double
handles.ne_pixels = str2double(get(hObject,'String'));
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function ne_pixels_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ne_pixels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function edge_margin_Callback(hObject, eventdata, handles)
% hObject    handle to edge_margin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edge_margin as text
%        str2double(get(hObject,'String')) returns contents of edge_margin as a double
handles.edge_margin = str2double(get(hObject,'String'));
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function edge_margin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edge_margin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function track_memory_Callback(hObject, eventdata, handles)
% hObject    handle to track_memory (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of track_memory as text
%        str2double(get(hObject,'String')) returns contents of track_memory as a double
handles.track_memory = str2double(get(hObject,'String'));
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function track_memory_CreateFcn(hObject, eventdata, handles)
% hObject    handle to track_memory (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function min_points_for_traj_Callback(hObject, eventdata, handles)
% hObject    handle to min_points_for_traj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of min_points_for_traj as text
%        str2double(get(hObject,'String')) returns contents of min_points_for_traj as a double
handles.min_points_for_traj = str2double(get(hObject,'String'));
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function min_points_for_traj_CreateFcn(hObject, eventdata, handles)
% hObject    handle to min_points_for_traj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fname_save_Callback(hObject, eventdata, handles)
% hObject    handle to fname_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fname_save as text
%        str2double(get(hObject,'String')) returns contents of fname_save as a double
handles.fname_save = get(hObject,'String');
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function fname_save_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fname_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function base_dir_Callback(hObject, eventdata, handles)
% hObject    handle to base_dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of base_dir as text
%        str2double(get(hObject,'String')) returns contents of base_dir as a double
handles.base_dir = get(hObject,'String');
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function base_dir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to base_dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in generate_metadata_parsed.
function generate_metadata_parsed_Callback(hObject, eventdata, handles)
% hObject    handle to generate_metadata_parsed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
metadata_conv_fname = [handles.ipdir,'times_from_umanager_metadata.py'];
phases = handles.phases;
imdir_phase = handles.imdir_phase;
for ph = [1:length(phases)]
      phase = phases{ph};
      imdir = imdir_phase.(phase);
      system(['python ', metadata_conv_fname, ' ', imdir]);
end

'parsed metadata files generated'


% --- Executes on button press in bgimg.
function bgimg_Callback(hObject, eventdata, handles)
% hObject    handle to bgimg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of bgimg


% --- Executes on button press in load_bg_image.
function load_bg_image_Callback(hObject, eventdata, handles)
% hObject    handle to load_bg_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%convert everything to handles
bgimg = get(handles.bgimg,'Value');
channel_to_image = handles.channel_to_image;
if bgimg == 0  %set default to 1
    for jj = 1:length(channel_to_image)
        %ch2i_txt = ch2i{jj}
        imbg.(channel_to_image{jj}) = 1;
    end
    
else
    for jj = 1:length(channel_to_image)
        imbg_jj = imread(handles.imbg_fname.(channel_to_image{jj}));
           
        if strcmp(handles.fname_conv,'JSO')
            imbg_jj = imbg_jj';  %micromanager images are the inverse of images collected by JSO image collection scripts.
        end
    
        %Convert bgimg to double and median filter
        imbg_jj = double(imbg_jj);
        coarse_smooth = 25;
        %smooths background image using course_smooth parameter.  Boundary
        %conditions are symmetric because default 0 bc's causes strange artifacts
        %on the edges.  For these background images symmetric BCs is a good
        %assumption
        imbg_jj = medfilt2(imbg_jj,[coarse_smooth,coarse_smooth],'symmetric');
        imbg.(channel_to_image{jj}) = imbg_jj;
    end
end

handles.imbg = imbg;

'bg image loaded'

set(handles.uipanel7,'Visible','on')

guidata(hObject, handles)


% --- Executes on button press in load_posvec_table.
function load_posvec_table_Callback(hObject, eventdata, handles)
% hObject    handle to load_posvec_table (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

posvec = get(handles.posvec_table,'Data')
handles.posvec = posvec;

'posvec table loaded'

set(handles.run_image_processing,'Visible','on')

guidata(hObject, handles)



function posvec_workspace_Callback(hObject, eventdata, handles)
% hObject    handle to load_posvec_workspace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of load_posvec_workspace as text
%        str2double(get(hObject,'String')) returns contents of load_posvec_workspace as a double


% --- Executes during object creation, after setting all properties.
function posvec_workspace_CreateFcn(hObject, eventdata, handles)
% hObject    handle to load_posvec_workspace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in load_posvec_workspace.
function load_posvec_workspace_Callback(hObject, eventdata, handles)
% hObject    handle to load_posvec_workspace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

posvec = evalin('base',get(handles.posvec_workspace,'String'))
handles.posvec = posvec;
set(handles.posvec_table,'Data',posvec)
'position data laoded from workspace'

set(handles.run_image_processing,'Visible','on')

guidata(hObject, handles)

% --- Executes on button press in run_image_processing.
function run_image_processing_Callback(hObject, eventdata, handles)
% hObject    handle to run_image_processing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

path(handles.ipdir,path)
analysis_params_names = {'ipdir', 
    'fname_conv',
    'storeim',
    'circ',
    'siz',
    'rad'
    'coarse_smooth'
    'local_smooth'
    'maxcells'
    'deconvlucy_iterations',
    'close_max',
    'far_max',
    'ne_pixels',
    'edge_margin',
    'maxdisp',
    'track_memory',
    'min_points_for_traj',
    'std_thresh',
    'thresh',
    'channels',
    'channel_to_image',
    'imbg'};

for jj = 1:length(analysis_params_names)
    analysis_params.(analysis_params_names{jj}) = handles.(analysis_params_names{jj});
end

%convert single channel to text from cell
channel_to_image = analysis_params.channel_to_image;
if  ( iscell(channel_to_image) & (length(channel_to_image)==1) )
    channel_to_image = channel_to_image{1};
    analysis_params.channel_to_image = channel_to_image;
end

all_tracks_vec = [];
all_times_vec = [];
posvec = handles.posvec;
Nwells = size(posvec,1);
phases = handles.phases;
time_calc_phase = handles.time_calc_phase;
imdir_phase = handles.imdir_phase;
shift_timing = handles.shift_timing;
base_dir = handles.base_dir;
fname_save = handles.fname_save;

for jj = 1:Nwells
   for ph = 1:length(phases);
        phase = phases{ph};
        analysis_params.imdir = imdir_phase.(phase);
        analysis_params.time_calc = time_calc_phase.(phase);
        pos_fnames = posvec(jj,:); 
        %remove bad positions (NAs)
        analysis_params.pos_fnames = pos_fnames(~strcmp(pos_fnames,'NA'));
        %main function for data processing
        [tracks,times] = time_series_analysis(analysis_params);
        all_tracks.(phase) = tracks;
        all_times.(phase) = times + shift_timing(ph);
   end
   all_tracks_vec{jj} = all_tracks;
   all_times_vec{jj} = all_times;
end
%store data

save([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')

%guidata(hObject, handles)


% --- Executes on button press in load_phase_table.
function load_phase_table_Callback(hObject, eventdata, handles)
% hObject    handle to load_phase_table (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
phase_table_data = get(handles.phase_table,'Data');
Nphases = find(1-strcmp('',{phase_table_data{1,:}}), 1, 'last' );
phases = {phase_table_data{1,1:Nphases}};
handles.phases = phases;

for jj = 1:length(phases);
    imdir_phase.(phases{jj})=[handles.base_dir,phase_table_data{2,jj},filesep];
end
handles.imdir_phase = imdir_phase;

handles.shift_timing = [phase_table_data{3,1:Nphases}];

for jj = 1:length(phases);
    time_calc_phase_jj = phase_table_data{4,jj};
    if regexpi(time_calc_phase_jj, '\d+');
        time_calc_phase_jj = str2num(time_calc_phase_jj);
    end
    time_calc_phase.(phases{jj})=time_calc_phase_jj;
end
handles.time_calc_phase = time_calc_phase;


handles.phases
handles.imdir_phase
handles.shift_timing
handles.time_calc_phase
'phase data loaded'

set(handles.uipanel3,'Visible','on')
set(handles.uipanel8,'Visible','on')

guidata(hObject, handles)


% --- Executes on button press in load_phase_table_workspace.
function load_phase_table_workspace_Callback(hObject, eventdata, handles)
% hObject    handle to load_phase_table_workspace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
phase_table = evalin('base',get(handles.phase_table_workspace,'String'))
set(handles.phase_table,'Data',phase_table);
'phase_table data laoded from workspace'

set(handles.uipanel3,'Visible','on')
set(handles.uipanel8,'Visible','on')

guidata(hObject, handles)

function phase_table_workspace_Callback(hObject, eventdata, handles)
% hObject    handle to phase_table_workspace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of phase_table_workspace as text
%        str2double(get(hObject,'String')) returns contents of phase_table_workspace as a double


% --- Executes during object creation, after setting all properties.
function phase_table_workspace_CreateFcn(hObject, eventdata, handles)
% hObject    handle to phase_table_workspace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in load_channel_table.
function load_channel_table_Callback(hObject, eventdata, handles)
% hObject    handle to load_channel_table (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
channel_table_data =  get(handles.channel_table,'Data');
Nchan = find(1-strcmp('',{channel_table_data{:,1}}), 1, 'last' );
channels = {channel_table_data{1:Nchan,1}};
handles.channels = channels;
channel_to_image_ind_vec = strcmp('True',{channel_table_data{1:Nchan,2}});
channel_to_image = channels(channel_to_image_ind_vec);

channel_to_image_ind = find(channel_to_image_ind_vec);
for jj = 1:length(channel_to_image)
    handles.imbg_fname.(channel_to_image{jj}) = [handles.base_dir,'BG',filesep,channel_table_data{channel_to_image_ind(jj),3}];
end

%if length(channel_to_image) == 1
%    channel_to_image = channel_to_image(1); %if only imaging one channel converts cell to text string 
%end

handles.channel_to_image = channel_to_image;

handles.channels
handles.channel_to_image
handles.imbg_fname
'Channel data loaded'

set(handles.bgimg,'Visible','on')
set(handles.load_bg_image,'Visible','on')

guidata(hObject, handles)


% --- Executes on button press in load_analysis_tracking_params.
function load_analysis_tracking_params_Callback(hObject, eventdata, handles)
% hObject    handle to load_analysis_tracking_params (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.coarse_smooth = str2double(get(handles.coarse_smooth,'String'));
handles.local_smooth = str2double(get(handles.local_smooth,'String'));
handles.maxcells = str2double(get(handles.maxcells,'String'));
handles.deconvlucy_iterations = str2double(get(handles.deconvlucy_iterations,'String'));
handles.close_max = str2double(get(handles.close_max,'String'));
handles.far_max = str2double(get(handles.far_max,'String'));
handles.ne_pixels = str2double(get(handles.ne_pixels,'String'));
handles.edge_margin = str2double(get(handles.edge_margin,'String'));
handles.thresh = str2double(get(handles.thresh,'String'));

handles.track_memory = str2double(get(handles.track_memory,'String'));
handles.min_points_for_traj = str2double(get(handles.min_points_for_traj,'String'));
handles.storeim = get(handles.storeim,'Value');

set(handles.uipanel5,'Visible','on');

guidata(hObject, handles)



function rad_Callback(hObject, eventdata, handles)
% hObject    handle to rad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rad as text
%        str2double(get(hObject,'String')) returns contents of rad as a double


% --- Executes during object creation, after setting all properties.
function rad_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function circ_file_Callback(hObject, eventdata, handles)
% hObject    handle to circ_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of circ_file as text
%        str2double(get(hObject,'String')) returns contents of circ_file as a double


% --- Executes during object creation, after setting all properties.
function circ_file_CreateFcn(hObject, eventdata, handles)
% hObject    handle to circ_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sizy_Callback(hObject, eventdata, handles)
% hObject    handle to sizy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sizy as text
%        str2double(get(hObject,'String')) returns contents of sizy as a double


% --- Executes during object creation, after setting all properties.
function sizy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sizy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxdisp_Callback(hObject, eventdata, handles)
% hObject    handle to maxdisp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxdisp as text
%        str2double(get(hObject,'String')) returns contents of maxdisp as a double


% --- Executes during object creation, after setting all properties.
function maxdisp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxdisp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in hes_presets.
function hes_presets_Callback(hObject, eventdata, handles)
% hObject    handle to hes_presets (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of hes_presets

if get(hObject,'Value') == 0
    %Enables preset parameters
    set(handles.circ_file,'Enable','on')
    set(handles.sizx,'Enable','on')
    set(handles.sizy,'Enable','on')
    set(handles.std_thresh,'Enable','on')
    set(handles.rad,'Enable','on')
    set(handles.maxdisp,'Enable','on')

    %shows load_species_magnification_params button
    set(handles.load_species_magnification_params,'Visible','on')

    %Hide parameters for preset program and program button
    set(handles.species_magnification_params, 'Visible','off')

elseif get(hObject,'Value') == 1
    %Disables preset parameters
    set(handles.circ_file,'Enable','off')
    set(handles.sizx,'Enable','off')
    set(handles.sizy,'Enable','off')
    set(handles.std_thresh,'Enable','off')
    set(handles.rad,'Enable','off')
    set(handles.maxdisp,'Enable','off')
    
    %hides load_species_magnification_params button
    set(handles.load_species_magnification_params,'Visible','off')
    
    %shows params for preset program and program button
    set(handles.species_magnification_params, 'Visible','on')    
end




function std_thresh_Callback(hObject, eventdata, handles)
% hObject    handle to std_thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of std_thresh as text
%        str2double(get(hObject,'String')) returns contents of std_thresh as a double


% --- Executes during object creation, after setting all properties.
function std_thresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to std_thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in load_species_magnification_params.
function load_species_magnification_params_Callback(hObject, eventdata, handles)
% hObject    handle to load_species_magnification_params (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

circ_loaded = load([handles.ipdir, get(handles.circ_file,'String')]);
handles.circ = circ_loaded.circ;
sizx = str2double(get(handles.sizx,'String'));
sizy = str2double(get(handles.sizy,'String'));
handles.siz = [sizx,sizy];
handles.rad = str2double(get(handles.rad,'String'));
handles.maxdisp = str2double(get(handles.maxdisp,'String'));
handles.std_thresh = str2double(get(handles.std_thresh,'String'));

handles.circ_file
handles.sizy
handles.rad
handles.maxdisp
handles.std_thresh
'Species / Magnification Parameters Set'

set(handles.uipanel4,'Visible','on')

guidata(hObject, handles)


function sizx_Callback(hObject, eventdata, handles)
% hObject    handle to sizx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sizx as text
%        str2double(get(hObject,'String')) returns contents of sizx as a double


% --- Executes during object creation, after setting all properties.
function sizx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sizx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Nsites_Callback(hObject, eventdata, handles)
% hObject    handle to Nsites (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Nsites as text
%        str2double(get(hObject,'String')) returns contents of Nsites as a double


% --- Executes during object creation, after setting all properties.
function Nsites_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Nsites (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in populate_pos_table.
function populate_pos_table_Callback(hObject, eventdata, handles)
% hObject    handle to populate_pos_table (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Nsites = str2double(get(handles.Nsites,'String'));
wellvec = get(handles.well_list,'Data');
wellvec = wellvec(not(strcmp('',wellvec)));
Nwells = length(wellvec);
for jj = 1:length(wellvec);
    for kk = 1:Nsites;
        posvec{jj,kk} = [wellvec{jj},'-Site_',num2str(kk-1)];
    end
end

set(handles.posvec_table,'Data',posvec);

guidata(hObject, handles)
