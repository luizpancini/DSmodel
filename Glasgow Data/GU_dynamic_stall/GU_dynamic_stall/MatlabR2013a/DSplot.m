function varargout = DSplot(varargin)
% DSPLOT M-datafile_select for DSplot.fig
%      DSPLOT, by itself, creates a new DSPLOT or raises the existing
%      singleton*.
%
%      H = DSPLOT returns the handle to a new DSPLOT or the handle to
%      the existing singleton*.
%
%      DSPLOT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DSPLOT.M with the given input arguments.
%
%      DSPLOT('Property','Value',...) creates a new DSPLOT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DSplot_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DSplot_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DSplot

% Last Modified by GUIDE v2.5 10-Sep-2013 13:19:37

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DSplot_OpeningFcn, ...
                   'gui_OutputFcn',  @DSplot_OutputFcn, ...
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



% --- Executes just before DSplot is made visible.
function DSplot_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DSplot (see VARARGIN)

% clc
handles.start_path = pwd;
handles.folder_name = 'datafolder_error';
handles.file_name = '';
if ispc
    handles.slash='\';
else
    handles.slash='/';
end

% Plot default
set(0,'DefaultAxesFontSize',15)
set(0,'DefaultAxesFontName','ItalicT')
set(0,'DefaultTextInterpreter','Latex')
sSize = get(0,'ScreenSize');
handles.sSize = sSize(3:4);

% Set default states
handles.loc_found = [1 0 0];
handles.coord_found = [1 0 0];
handles.RIBs_found = [1 0 0];
handles.enablingFolder = [handles.datafile_select handles.datafile_type...
                          handles.datafile_typenumber handles.slider_typenumber...
                          handles.clear_data...
                          handles.list_models handles.list_experiments handles.list_motions...
                          handles.selectedfiles_select handles.selectedfiles_info];
handles.enablingFile = [handles.export_coeffs...
                        handles.datafile_open handles.plot_calib handles.rib...
                        handles.plot_profile handles.plot_taps...
                        handles.time_scaled handles.time_seconds...
                        handles.plot_angle handles.plot_wind...
                        handles.plot_cp2D handles.side_u handles.side_ul...
                        handles.plot_cp3D handles.surf_cp...
                        handles.plot_coeffs handles.c_nmt handles.c_lmd...
                        handles.plot_close_all];
set(handles.enablingFolder, 'Enable', 'Off');
set(handles.enablingFile, 'Enable', 'Off');

set(handles.selected_freestream, 'Enable', 'Off', 'Value', 0);
set(handles.selected_motion, 'Enable', 'Off', 'Value', 1);

% Choose default command line output for DSplot
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes DSplot wait for user response (see UIRESUME)
% uiwait(handles.figure1);



% --- Outputs from this function are returned to the command line.
function varargout = DSplot_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% -------------------------------------------------------------------------
% FOLDER
% -------------------------------------------------------------------------

% --- Executes on button press in datafolder_select ('Select data folder').
function datafolder_select_Callback(hObject, eventdata, handles)
% hObject    handle to datafolder_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

folder_name = uigetdir(handles.start_path);

if folder_name ~= 0
    
    handles.folder_name = folder_name;
    handles.start_path = handles.folder_name;
    set(handles.datafolder_view,'String',handles.folder_name);
    set(handles.datafile_typeWARNING,'String',[]);

    % Number of cases (.dat files)
    d = dir(strcat(handles.folder_name,handles.slash,'*.dat'));
    N = size(d,1);
    
    if N < 4 % N == 0
        set(handles.datafolder_selectWARNING,'String','No data files found.');
        handles.folder_name = 'datafolder_error';
        set(handles.datafiles_no,'String','');
        set(handles.enablingFolder, 'Enable', 'Off');
        set(handles.enablingFile, 'Enable', 'Off');
        handles.loc_found = [1 0 0];
        handles.coord_found = [1 0 0];
        handles.RIBs_found = [1 0 0];
    else
        set(handles.datafolder_selectWARNING,'String',[]);
        set(handles.datafile_typeWARNING,'String',[]);
        handles.count_ab = ones(1,N);
        handles.count_c = ones(1,N);
        handles.count_d = ones(1,N);
        set(handles.datafiles_no,'String',num2str(N));
        set(handles.cases_no,'String',num2str(N));
        set(handles.enablingFolder, 'Enable', 'On');
        set(handles.slider_typenumber,'sliderstep',[1 N+1]/(N+1),'max',N+1,'min',0,'Value',0)

        % Checking the pressure transducer locations file
        trans_locs_filename = get(handles.trans_locs_filename,'String');
        FID = fopen(strcat(handles.folder_name,handles.slash,trans_locs_filename),'r');
        if FID > 0
            handles.loc_found = [0 1 0];
            handles.locfile = strcat(handles.folder_name,handles.slash,trans_locs_filename);
        else
            slash = findstr(handles.folder_name,handles.slash);
            FID = fopen(strcat(handles.folder_name(1:slash(end)),trans_locs_filename),'r');
            if FID > 0
                handles.loc_found = [0 1 0];
                handles.locfile = strcat(handles.folder_name(1:slash(end)),trans_locs_filename);
            else
                handles.loc_found = [1 0 0];
                handles.locfile = '';
            end
        end

        % Checking the aerofoil coordinate database file
        coord_filename = get(handles.coords_filename,'String');
        FID = fopen(strcat(handles.folder_name,handles.slash,coord_filename),'r');
        if FID > 0
            handles.coord_found = [0 1 0];
            handles.coordfile = strcat(handles.folder_name,handles.slash,coord_filename);
        else
            slash = findstr(handles.folder_name,handles.slash);
            FID = fopen(strcat(handles.folder_name(1:slash(end)),coord_filename),'r');
            if FID > 0
                handles.coord_found = [0 1 0];
                handles.coordfile = strcat(handles.folder_name(1:slash(end)),coord_filename);
            else
                handles.coord_found = [1 0 0];
                handles.coordfile = '';
            end
        end
        
        if all(handles.loc_found==[0 1 0]) * all(handles.coord_found==[0 1 0]) == 1
            set(handles.datafolder_selectWARNING,'String',[]);
        elseif all(handles.loc_found==[0 1 0]) + all(handles.coord_found==[0 1 0]) == 0
            set(handles.datafolder_selectWARNING,'String','No transducer locations and aerofoil coordinates files found.');
        elseif all(handles.loc_found==[0 1 0])==0
            set(handles.datafolder_selectWARNING,'String','No transducer locations file found.');
        else
            set(handles.datafolder_selectWARNING,'String','No aerofoil coordinates file found.');
        end
        
        % Checking the RIBs info file
        RIBs_filename = get(handles.RIBs_filename,'String');
        if exist(strcat(handles.folder_name,handles.slash,RIBs_filename),'file') == 2
            handles.RIBs_found = [0 1 0];
            set(handles.selected_freestream, 'Enable', 'On');
            set(handles.selected_motion, 'Enable', 'On');
        else
            slash = findstr(handles.folder_name,handles.slash);
            if exist(strcat(handles.folder_name(1:slash(end)),RIBs_filename),'file') == 2
                handles.RIBs_found = [0 1 0];
                set(handles.selected_freestream, 'Enable', 'On');
                set(handles.selected_motion, 'Enable', 'On');
            else
                handles.RIBs_found = [1 0 0];
                set(handles.selected_freestream, 'Enable', 'Off');
                set(handles.selected_motion, 'Enable', 'Off');
            end
        end
        
    end
    
    set(handles.loc,'BackgroundColor',handles.loc_found);
    set(handles.coords,'BackgroundColor',handles.coord_found);
    set(handles.RIBs,'BackgroundColor',handles.RIBs_found);
    
    % Set texts as default values
    set(handles.datafile_view,'String','Run number');
    set(handles.datafile_type,'String','Type run number');
    set(handles.datafile_typenumber,'String','Type list number');
    set(handles.datafile_typeWARNING,'String',[]);
    set(handles.INFO_model,'String',[]);
    set(handles.INFOname_model,'String',[]);
    set(handles.INFO_experiment,'String',[]);
    set(handles.INFOname_experiment,'String',[]);
    set(handles.INFO_motion,'String',[]);
    set(handles.INFOname_motion,'String',[]);
    set(handles.INFO_Re,'String',[]);
    set(handles.INFOname_Re,'String',[]);
    set(handles.INFO_ForK,'String',[]);
    set(handles.INFOname_ForK,'String',[]);

    % Delete data
    handles.data = [];
    handles.file_name = '';
    handles.chord = [];
    handles.location = [];
    handles.coordinate = [];
    
    fclose all;
% else
%     handles.folder_name = 'datafolder_error';
end

guidata(hObject, handles);



function trans_locs_filename_Callback(hObject, eventdata, handles)
% hObject    handle to trans_locs_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trans_locs_filename as text
%        str2double(get(hObject,'String')) returns contents of trans_locs_filename as a double

% Checking the pressure transducer locations file
trans_locs_filename = get(hObject,'String');
FID = fopen(strcat(handles.folder_name,handles.slash,trans_locs_filename),'r');
if FID > 0
    handles.loc_found = [0 1 0];
    handles.locfile = strcat(handles.folder_name,handles.slash,trans_locs_filename);
else
    slash = findstr(handles.folder_name,handles.slash);
    FID = fopen(strcat(handles.folder_name(1:slash(end)),trans_locs_filename),'r');
    if FID > 0
        handles.loc_found = [0 1 0];
        handles.locfile = strcat(handles.folder_name(1:slash(end)),trans_locs_filename);
    else
        handles.loc_found = [1 0 0];
        handles.locfile = '';
    end
end
set(handles.loc,'BackgroundColor',handles.loc_found);

if all(handles.loc_found==[0 1 0]) * all(handles.coord_found==[0 1 0]) == 1
    set(handles.datafolder_selectWARNING,'String',[]);
elseif all(handles.loc_found==[0 1 0]) + all(handles.coord_found==[0 1 0]) == 0
    set(handles.datafolder_selectWARNING,'String','No transducer locations and aerofoil coordinates files found.');
elseif all(handles.loc_found==[0 1 0])==0
    set(handles.datafolder_selectWARNING,'String','No transducer locations file found.');
else
    set(handles.datafolder_selectWARNING,'String','No aerofoil coordinates file found.');
end

fclose all;

guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function trans_locs_filename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trans_locs_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function coords_filename_Callback(hObject, eventdata, handles)
% hObject    handle to coords_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of coords_filename as text
%        str2double(get(hObject,'String')) returns contents of coords_filename as a double

% Checking the aerofoil coordinate database file
coord_filename = get(hObject,'String');
FID = fopen(strcat(handles.folder_name,handles.slash,coord_filename),'r');
if FID > 0
    handles.coord_found = [0 1 0];
    handles.coordfile = strcat(handles.folder_name,handles.slash,coord_filename);
else
    slash = findstr(handles.folder_name,handles.slash);
    FID = fopen(strcat(handles.folder_name(1:slash(end)),coord_filename),'r');
    if FID > 0
        handles.coord_found = [0 1 0];
        handles.coordfile = strcat(handles.folder_name(1:slash(end)),coord_filename);
    else
        handles.coord_found = [1 0 0];
        handles.coordfile = '';
    end
end
set(handles.coords,'BackgroundColor',handles.coord_found);

if all(handles.loc_found==[0 1 0]) * all(handles.coord_found==[0 1 0]) == 1
    set(handles.datafolder_selectWARNING,'String',[]);
elseif all(handles.loc_found==[0 1 0]) + all(handles.coord_found==[0 1 0]) == 0
    set(handles.datafolder_selectWARNING,'String','No transducer locations and aerofoil coordinates files found.');
elseif all(handles.loc_found==[0 1 0])==0
    set(handles.datafolder_selectWARNING,'String','No transducer locations file found.');
else
    set(handles.datafolder_selectWARNING,'String','No aerofoil coordinates file found.');
end

fclose all;

guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function coords_filename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to coords_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function RIBs_filename_Callback(hObject, eventdata, handles)
% hObject    handle to RIBs_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RIBs_filename as text
%        str2double(get(hObject,'String')) returns contents of RIBs_filename as a double

RIBs_filename = get(hObject,'String');
if exist(strcat(handles.folder_name,handles.slash,RIBs_filename),'file') == 2
    handles.RIBs_found = [0 1 0];
    set(handles.selected_freestream, 'Enable', 'On');
    set(handles.selected_motion, 'Enable', 'On');
else
    slash = findstr(handles.folder_name,handles.slash);
    if exist(strcat(handles.folder_name(1:slash(end)),RIBs_filename),'file') == 2
       handles.RIBs_found = [0 1 0];
       set(handles.selected_freestream, 'Enable', 'On');
       set(handles.selected_motion, 'Enable', 'On');
    else
       handles.RIBs_found = [1 0 0];
       set(handles.selected_freestream, 'Enable', 'Off');
       set(handles.selected_motion, 'Enable', 'Off');
    end
end
set(handles.RIBs,'BackgroundColor',handles.RIBs_found);

guidata(hObject, handles);
        
        
% --- Executes during object creation, after setting all properties.
function RIBs_filename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RIBs_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in selectedfiles_info.
function selectedfiles_info_Callback(hObject, eventdata, handles)
% hObject    handle to selectedfiles_info (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

d = dir(strcat(handles.folder_name,handles.slash,'*.dat'));
N = size(d,1);

RIBs_info = zeros(N,32);
for i = 1 : N
    disp([num2str(i), ' of ', num2str(N), ': ', d(i).name])
    a = importdata(strcat(handles.folder_name,handles.slash,d(i).name));
    RIBs_info(i,:) = a(1,1:32);
end

slash = findstr(handles.folder_name,handles.slash);
filename = strcat(handles.folder_name(1:slash(end)),'RIBs_info.mat');
save(filename, 'RIBs_info')
handles.RIBs_found = [0 1 0];
set(handles.RIBs,'BackgroundColor',handles.RIBs_found);
set(handles.RIBs_filename,'String','RIBs_info.mat');
set(handles.selected_freestream, 'Enable', 'On');
set(handles.selected_motion, 'Enable', 'On');

guidata(hObject, handles);



% -------------------------------------------------------------------------
% DATA
% -------------------------------------------------------------------------

% ------------------------------------
% Select
% ------------------------------------

% --- Executes on button press in datafile_select  ('Select data file').
function datafile_select_Callback(hObject, eventdata, handles)
% hObject    handle to datafile_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if strcmp(handles.folder_name,'datafolder_error') == 1
    set(handles.datafile_typeWARNING,'String','Select valid data folder first.');
else
    set(handles.datafile_typeWARNING,'String',[]);
    handles.file_name = uigetfile('.dat','',handles.folder_name);

    if handles.file_name ~= 0
        set(handles.datafile_view,'String',handles.file_name(1:end-4));
        handles.data = importdata(strcat(handles.folder_name,handles.slash,handles.file_name));
        set(handles.enablingFile, 'Enable', 'On');
        if size(handles.data,1)~=handles.data(1,11)+4
            set(handles.plot_calib, 'Enable', 'Off');
        end

        % Finding the list number
        d = dir(strcat(handles.folder_name,handles.slash,'*.dat'));
        j = 1;
        found = 0;
        while found == 0
            if strcmp(d(j).name, handles.file_name) == 1
                found = 1;
                set(handles.datafile_typenumber,'String',num2str(j));
                set(handles.slider_typenumber,'Value',j);
            end
            j = j+1;
        end
    
        % INFO boxes
        models = get(handles.list_models,'String');
        set(handles.INFO_model,'String',models{str2double(handles.file_name(1:2))+1}(4:end));
        set(handles.INFOname_model,'String','Model');
        experiment = get(handles.list_experiments,'String');
        set(handles.INFO_experiment,'String',experiment{str2double(handles.file_name(3))+2}(3:end));
        set(handles.INFOname_experiment,'String','Type of experiment');
        motion = get(handles.list_motions,'String');
        set(handles.INFO_motion,'String',motion{str2double(handles.file_name(4))+2}(3:end));
        set(handles.INFOname_motion,'String','Motion type');
        set(handles.INFO_Re,'String',handles.data(1,20));
        set(handles.INFOname_Re,'String','Reynolds number');
        if handles.data(1,7) == 1 || handles.data(1,7) == 5 || handles.data(1,7) == 7
            set(handles.INFO_ForK,'String',handles.data(1,22));
            set(handles.INFOname_ForK,'String','Reduced frequency');
        elseif handles.data(1,7) == 2 || handles.data(1,7) == 3
            set(handles.INFO_ForK,'String',handles.data(1,22));
            set(handles.INFOname_ForK,'String','Reduced pitch-rate');
        else
            set(handles.INFO_ForK,'String',[]);
            set(handles.INFOname_ForK,'String',[]);
        end

        % Definition of the chord of the model
        if handles.data(1,30) == 12
            handles.chord = 0.55/2;
        else
            handles.chord = 0.55;
        end

        % Loading pressure transducer locations
        FID = fopen(handles.locfile,'r');
        if FID > 0
            locations = textscan(FID, '%s', 'Delimiter', '\n');
            handles.location = str2num(locations{1}{handles.data(1,30)+1});
        end

        % Loading aerofoil coordinates
        coordinates = xlsread(handles.coordfile);
        handles.coordinate = coordinates(:,handles.data(1,30)*4-3:handles.data(1,30)*4);
        handles.coordinate = handles.coordinate(isfinite(handles.coordinate(:,1)),:);
    else
        handles.file_name = '';
    end
    
    fclose all;
end

guidata(hObject, handles);



function datafile_type_Callback(hObject, eventdata, handles)
% hObject    handle to datafile_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of datafile_type as text
%        str2double(get(hObject,'String')) returns contents of datafile_type as a double

FID = fopen(strcat(handles.folder_name,handles.slash,get(hObject,'String'),'.dat'));

if strcmp(handles.folder_name,'datafolder_error') == 1
    set(handles.datafile_typeWARNING,'String','Select valid data folder first.');
    set(hObject,'String','Type run number');
elseif FID < 0
    set(handles.datafile_typeWARNING,'String','Typed case not found.');
else
    set(handles.datafile_typeWARNING,'String',[]);
    handles.data = importdata(strcat(handles.folder_name,handles.slash,get(hObject,'String'),'.dat'));
    set(handles.enablingFile, 'Enable', 'On');
    if size(handles.data,1)~=handles.data(1,11)+4
        set(handles.plot_calib, 'Enable', 'Off');
    end

    set(handles.datafile_view,'String',get(hObject,'String'));
    handles.file_name = strcat(get(hObject,'String'), '.dat');
    
    % Finding the list number
    d = dir(strcat(handles.folder_name,handles.slash,'*.dat'));
    j = 1;
    found = 0;
    while found == 0
        if strcmp(d(j).name, handles.file_name) == 1
            found = 1;
            set(handles.datafile_typenumber,'String',num2str(j));
            set(handles.slider_typenumber,'Value',j);
        end
        j = j+1;
    end
    
    % INFO boxes
    models = get(handles.list_models,'String');
    set(handles.INFO_model,'String',models{str2double(handles.file_name(1:2))+1}(4:end));
    set(handles.INFOname_model,'String','Model');
    experiment = get(handles.list_experiments,'String');
    set(handles.INFO_experiment,'String',experiment{str2double(handles.file_name(3))+2}(3:end));
    set(handles.INFOname_experiment,'String','Type of experiment');
    motion = get(handles.list_motions,'String');
    set(handles.INFO_motion,'String',motion{str2double(handles.file_name(4))+2}(3:end));
    set(handles.INFOname_motion,'String','Motion type');
    set(handles.INFO_Re,'String',handles.data(1,20));
    set(handles.INFOname_Re,'String','Reynolds number');
    if handles.data(1,7) == 1 || handles.data(1,7) == 5 || handles.data(1,7) == 7
        set(handles.INFO_ForK,'String',handles.data(1,22));
        set(handles.INFOname_ForK,'String','Reduced frequency');
    elseif handles.data(1,7) == 2 || handles.data(1,7) == 3
        set(handles.INFO_ForK,'String',handles.data(1,22));
        set(handles.INFOname_ForK,'String','Reduced pitch-rate');
    else
        set(handles.INFO_ForK,'String',[]);
        set(handles.INFOname_ForK,'String',[]);
    end
        
    % Definition of the chord of the model
    if handles.data(1,30) == 12
        handles.chord = 0.55/2;
    else
        handles.chord = 0.55;
    end

    % Loading pressure transducer locations
    FID = fopen(handles.locfile,'r');
    if FID > 0
        locations = textscan(FID, '%s', 'Delimiter', '\n');
        handles.location = str2num(locations{1}{handles.data(1,30)+1});
    end

    % Loading aerofoil coordinates
    coordinates = xlsread(handles.coordfile);
    handles.coordinate = coordinates(:,handles.data(1,30)*4-3:handles.data(1,30)*4);
    handles.coordinate = handles.coordinate(isfinite(handles.coordinate(:,1)),:);
end

fclose all;

guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function datafile_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to datafile_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function datafile_typenumber_Callback(hObject, eventdata, handles)
% hObject    handle to datafile_typenumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of datafile_typenumber as text
%        str2double(get(hObject,'String')) returns contents of datafile_typenumber as a double

d = dir(strcat(handles.folder_name,handles.slash,'*.dat'));
N = size(d,1);

if strcmp(handles.folder_name,'datafolder_error') == 1
    set(handles.datafile_typeWARNING,'String','Select valid data folder first.');
    set(hObject,'String','Type list number');
elseif isnan(str2double(get(hObject,'String'))) == 1
    set(handles.datafile_typeWARNING,'String','Invalid list number.');
    set(hObject,'String','Type list number');
elseif str2double(get(hObject,'String')) > N
    set(handles.datafile_typeWARNING,'String','Listing number has to be equal or less than the total number of data files.');
    set(handles.slider_typenumber,'Value',N+1)
elseif str2double(get(hObject,'String')) < 1
    set(handles.datafile_typeWARNING,'String','Listing number has to be greater than zero.');
    set(handles.slider_typenumber,'Value',0)
else
    set(handles.datafile_typeWARNING,'String',[]);
    handles.data = importdata(strcat(handles.folder_name,handles.slash,d(str2double(get(hObject,'String'))).name));
    set(handles.enablingFile, 'Enable', 'On');
    if size(handles.data,1)~=handles.data(1,11)+4
        set(handles.plot_calib, 'Enable', 'Off');
    end

    set(handles.datafile_view,'String',d(str2double(get(hObject,'String'))).name(1:end-4));
    handles.file_name = d(str2double(get(hObject,'String'))).name;
    set(handles.slider_typenumber,'Value',str2double(get(hObject,'String')))
    
    % INFO boxes
    models = get(handles.list_models,'String');
    set(handles.INFO_model,'String',models{str2double(handles.file_name(1:2))+1}(4:end));
    set(handles.INFOname_model,'String','Model');
    experiment = get(handles.list_experiments,'String');
    set(handles.INFO_experiment,'String',experiment{str2double(handles.file_name(3))+2}(3:end));
    set(handles.INFOname_experiment,'String','Type of experiment');
    motion = get(handles.list_motions,'String');
    set(handles.INFO_motion,'String',motion{str2double(handles.file_name(4))+2}(3:end));
    set(handles.INFOname_motion,'String','Motion type');
    set(handles.INFO_Re,'String',handles.data(1,20));
    set(handles.INFOname_Re,'String','Reynolds number');
    if handles.data(1,7) == 1 || handles.data(1,7) == 5 || handles.data(1,7) == 7
        set(handles.INFO_ForK,'String',handles.data(1,22));
        set(handles.INFOname_ForK,'String','Reduced frequency');
    elseif handles.data(1,7) == 2 || handles.data(1,7) == 3
        set(handles.INFO_ForK,'String',handles.data(1,22));
        set(handles.INFOname_ForK,'String','Reduced pitch-rate');
    else
        set(handles.INFO_ForK,'String',[]);
        set(handles.INFOname_ForK,'String',[]);
    end
    
    % Definition of the chord of the model
    if handles.data(1,30) == 12
        handles.chord = 0.55/2;
    else
        handles.chord = 0.55;
    end

    % Loading pressure transducer locations
    FID = fopen(handles.locfile,'r');
    if FID > 0
        locations = textscan(FID, '%s', 'Delimiter', '\n');
        handles.location = str2num(locations{1}{handles.data(1,30)+1});
    end

    % Loading aerofoil coordinates
    coordinates = xlsread(handles.coordfile);
    handles.coordinate = coordinates(:,handles.data(1,30)*4-3:handles.data(1,30)*4);
    handles.coordinate = handles.coordinate(isfinite(handles.coordinate(:,1)),:);
    
    fclose all;
end

guidata(hObject, handles);
    

       
% --- Executes during object creation, after setting all properties.
function datafile_typenumber_CreateFcn(hObject, eventdata, handles)
% hObject    handle to datafile_typenumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on slider movement.
function slider_typenumber_Callback(hObject, eventdata, handles)
% hObject    handle to slider_typenumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

d = dir(strcat(handles.folder_name,handles.slash,'*.dat'));
N = size(d,1);

slider_number = round(get(hObject,'Value'));
set(hObject,'Value',slider_number);

if slider_number > 0 && slider_number < N+1
    
    set(handles.datafile_typeWARNING,'String',[]);
    handles.data = importdata(strcat(handles.folder_name,handles.slash,d(slider_number).name));
    set(handles.enablingFile, 'Enable', 'On');
    if size(handles.data,1)~=handles.data(1,11)+4
        set(handles.plot_calib, 'Enable', 'Off');
    end

    set(handles.datafile_view,'String',d(slider_number).name(1:end-4));
    set(handles.datafile_typenumber,'String',num2str(slider_number));
    handles.file_name = d(slider_number).name;
    
    
    % INFO boxes
    models = get(handles.list_models,'String');
    set(handles.INFO_model,'String',models{str2double(handles.file_name(1:2))+1}(4:end));
    set(handles.INFOname_model,'String','Model');
    experiment = get(handles.list_experiments,'String');
    set(handles.INFO_experiment,'String',experiment{str2double(handles.file_name(3))+2}(3:end));
    set(handles.INFOname_experiment,'String','Type of experiment');
    motion = get(handles.list_motions,'String');
    set(handles.INFO_motion,'String',motion{str2double(handles.file_name(4))+2}(3:end));
    set(handles.INFOname_motion,'String','Motion type');
    set(handles.INFO_Re,'String',handles.data(1,20));
    set(handles.INFOname_Re,'String','Reynolds number');
    if handles.data(1,7) == 1 || handles.data(1,7) == 5 || handles.data(1,7) == 7
        set(handles.INFO_ForK,'String',handles.data(1,22));
        set(handles.INFOname_ForK,'String','Reduced frequency');
    elseif handles.data(1,7) == 2 || handles.data(1,7) == 3
        set(handles.INFO_ForK,'String',handles.data(1,22));
        set(handles.INFOname_ForK,'String','Reduced pitch-rate');
    else
        set(handles.INFO_ForK,'String',[]);
        set(handles.INFOname_ForK,'String',[]);
    end
    
    % Definition of the chord of the model
    if handles.data(1,30) == 12
        handles.chord = 0.55/2;
    else
        handles.chord = 0.55;
    end

    % Loading pressure transducer locations
    FID = fopen(handles.locfile,'r');
    if FID > 0
        locations = textscan(FID, '%s', 'Delimiter', '\n');
        handles.location = str2num(locations{1}{handles.data(1,30)+1});
    end

    % Loading aerofoil coordinates
    coordinates = xlsread(handles.coordfile);
    handles.coordinate = coordinates(:,handles.data(1,30)*4-3:handles.data(1,30)*4);
    handles.coordinate = handles.coordinate(isfinite(handles.coordinate(:,1)),:);
    
    fclose all;
end

guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function slider_typenumber_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_typenumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



% ------------------------------------
% Export
% ------------------------------------

% --- Executes on button press in export_coeffs.
function export_coeffs_Callback(hObject, eventdata, handles)
% hObject    handle to export_coeffs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if size(handles.data,1)==handles.data(1,11)+4
    dat_line = 5;
else
    dat_line = 2;
end

alpha = handles.data(dat_line:end,end);

T = [0:1:size(handles.data,1)-dat_line]*(1/handles.data(1,18));
if handles.data(1,7) == 1 || handles.data(1,7) == 5
    nonT = T * handles.data(1,10) * 2 * pi;
else
    nonT = T * handles.data(1,23) / handles.chord;
end

xs = [handles.coordinate(:,1), handles.coordinate(:,3)]/100;
ys = [handles.coordinate(:,2), handles.coordinate(:,4)]/100;
LE_u = find((handles.location(2:end)-handles.location(1:end-1))<0,1,'last')+1;
xp = [handles.location', [ones(LE_u,1); -ones(length(handles.location)-LE_u,1)]];

[yp, dyp, lp] = aerofoil_taps(xs, ys, xp);

% Lift coefficient calculation
cy = zeros(size(handles.data,1)-dat_line+1,1);
cx = zeros(size(handles.data,1)-dat_line+1,1);
cm_x = zeros(size(handles.data,1)-dat_line+1,1);

for i = 1 : length(alpha)
    
    angle = atan(dyp(:,1));
    angle = angle+pi/2;
    find_ls = find(dyp(:,2)==-1,1,'first');
    angle(find_ls:end) = pi+angle(find_ls:end,:);
    
    % Data
    cp = (handles.data(i+dat_line-1,2:end-1))';
    
    cy_i = -cp.*lp.*sin(angle);
    cx_i = -cp.*lp.*cos(angle);
    
    % Forces normal and parallel to the chord
    cy(i) = sum(cy_i);
    cx(i) = -sum(cx_i);
    
    % Moment about the quarter-chord
    cm_x(i) = sum(-cy_i.*(xp(:,1)-0.25) + cx_i.*yp(:,1));
    
end

head = '% Time() Angle(deg) Cn Ct Cm' ;
data = [round(nonT'*100)/100 alpha cy cx cm_x];
fName = [handles.file_name(1:8), '_coeffs', '.dat'];
fid = fopen(fName,'w');
if fid >= 0
    fprintf(fid, '%s\r\n', head);
    fclose(fid);
end
dlmwrite(fName, data, '-append', 'newline', 'pc', 'delimiter','\t');
disp(['Data exported in ', fName])



% ------------------------------------
% Info
% ------------------------------------

% --- Executes on button press in datafile_open ('Open data file').
function datafile_open_Callback(hObject, eventdata, handles)
% hObject    handle to datafile_open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if strcmp(handles.folder_name,'datafolder_error') == 1
    set(handles.datafile_typeWARNING,'String','Select valid data folder first.');
elseif strcmp(handles.file_name,'')==1
    set(handles.datafile_typeWARNING,'String','Select valid data file.');
else
    winopen([handles.folder_name handles.slash handles.file_name])
end



% --- Executes on button press in plot_calib ('Plot calib. values').
function plot_calib_Callback(hObject, eventdata, handles)
% hObject    handle to plot_calib (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if strcmp(handles.folder_name,'datafolder_error') == 1
    set(handles.datafile_typeWARNING,'String','Select valid data folder first.');
elseif strcmp(handles.file_name,'')==1
    set(handles.datafile_typeWARNING,'String','Select valid data file.');
else

    nCh = length(handles.data(1,:));
    figure
    
    [ax,h1,h2] = plotyy([1:1:nCh], handles.data(2,:), [1:1:nCh], handles.data(3,:));
    set(ax,{'ycolor'},{'k';'k'})

    % Calibration Values
    set(h1,'LineStyle','none','Marker','o','LineWidth',2,'Color','b')
    set(get(ax(1),'Ylabel'),'String','Calibration Values')
    set(get(ax(1),'Ylabel'),'FontSize',15)
    set(get(ax(1),'Ylabel'),'Color','b')
%     ylim(ax(1), [0.65 0.7])
%     set(ax(1),'YTick',[0.65:0.01:0.7])
    xlim(ax(1),[0 nCh+1])
    set(ax(1),'XTick',[0:5:40])

    % Gain Values
    set(h2,'LineStyle','none','Marker','d','LineWidth',2,'Color','r')
    set(get(ax(2),'Ylabel'),'String','Gain Values')
    set(get(ax(2),'Ylabel'),'FontSize',15)
    set(get(ax(2),'Ylabel'),'Color','r')
%     ylim(ax(2), [0 10])
%     set(ax(2),'YTick',[0:2:10])
    xlim(ax(2),[0 nCh+1])
    set(ax(2),'XTick',[0:5:40])
    
    % Offset Values
    hold on
    plot([0 nCh+1], [1 1], '--k', 'LineWidth',1)
    plot([1:1:nCh], handles.data(4,:), 'sb', 'LineWidth', 2)

    xlabel('Channel number')
    set(gca,'OuterPosition',[-0.11 0.02 1.1 1])
    set(gcf,'Name',handles.file_name)
    
end



% --- Executes on button press in rib ('RIB report').
function rib_Callback(hObject, eventdata, handles)
% hObject    handle to rib (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if strcmp(handles.folder_name,'datafolder_error') == 1
    set(handles.datafile_typeWARNING,'String','Select valid data folder first.');
elseif strcmp(handles.file_name,'')==1
    set(handles.datafile_typeWARNING,'String','Select valid data file.');
else
    nCh = length(handles.data(1,:));
    dat(:,1) = handles.data(1,:)';

    cnames = {' ','Ch. Number', 'Cal. Value','Gain', 'Offset'};
    rnames = {'Run Number', 'Date of Test: Day', 'Date of Test: Month', 'Date of Test: Year',...
              'Temperature (Celsius)', 'Barometric Pressure (mmHg)', 'Motion Type', 'Empty', 'Empty', 'Empty',...
              'Number of Sweeps per Cycle', 'Number of Values per Cycle', 'Number of Cycles',...
              'Empty', 'Empty', 'Empty', 'Empty',...
              'Sampling Frequency (Hz)', 'Dynamic Pressure (Volts)', 'Reynolds Number', 'Mach Number',...
              'Empty', 'Wind Velocity (m/s)', 'Empty', 'Empty', 'Averaged (1) or Unaveraged (0)'...
              'Empty', 'Empty', 'Dynamic Pressure (N/m2)', 'Model Number', 'Empty', 'Empty'};
    
    if dat(7,1) == 0
        rnames(:,8:9) = {'Starting Incidence (deg)', 'Arc (deg)'};
    elseif dat(7,1) == 4
        rnames(:,8) = {'Nominal Angle of Attack (deg)'};
        rnames(:,9) = {'Averaged Angle of Attack (deg)'};
    elseif dat(7,1) == 1 || dat(7,1) == 5
        rnames(:,8:10) = {'Mean Incidence (deg)', 'Amplitude (deg)', 'Oscillation Frequency (Hz)'};
        rnames(:,22) = {'Reduced Frequency'};
    elseif dat(7,1) == 2 || dat(7,1) == 3
        rnames(:,8:10) = {'Starting Incidence (deg)', 'Ramp Arc (deg)', 'Linear Pitch-Rate (deg/s)'};
        rnames(:,22) = {'Reduced Pitch-Rate'};
    elseif dat(7,1) == 7
        rnames(:,8:10) = {'Mean Incidence (deg)', 'Tip Speed Ratio', 'Oscillation Frequency (Hz)'};
        rnames(:,22) = {'Reduced Frequency'};
    end
    
    if nCh == 37
        rnames = [rnames, 'Empty', 'Empty', 'Empty', 'Empty', 'Empty'];
        posF4 = 720;
        posT4 = 688;
    else
        posF4 = 630;
        posT4 = 598;
    end
    
    if size(handles.data,1)==handles.data(1,11)+4
        posF3 = 732;
        posT3 = 700;
        dat(:,2) = [1:1:nCh]';
        dat(:,3:5) = handles.data(2:4,:)';
    else
        posF3 = 452;
        posT3 = 420;
        cnames = cnames{1};
        rnames(:,19) = rnames(:,29);
        rnames(:,29) = {'Empty'};
    end
    
    posF = [100 100 posF3 posF4];
    posT = [20 20 posT3 posT4];
    
    format short g
    f = figure('Position',posF);
    t = uitable('Parent',f,'Data',dat,...
                'ColumnName',cnames, 'ColumnWidth',{70 70 70 70 70}, 'ColumnFormat', {'short g',[],[],[],[]},...
                'RowName',rnames,'Position',posT);

    set(gcf,'Name',handles.file_name)
end



% --- Executes on button press in clear_data ('Reset').
function clear_data_Callback(hObject, eventdata, handles)
% hObject    handle to clear_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Set texts as default values
set(handles.datafile_view,'String','Run number');
set(handles.datafile_type,'String','Type run number');
set(handles.datafile_typenumber,'String','Type list number');
set(handles.slider_typenumber,'Value',0)
set(handles.datafile_typeWARNING,'String',[]);
set(handles.INFO_model,'String',[]);
set(handles.INFOname_model,'String',[]);
set(handles.INFO_experiment,'String',[]);
set(handles.INFOname_experiment,'String',[]);
set(handles.INFO_motion,'String',[]);
set(handles.INFOname_motion,'String',[]);
set(handles.INFO_Re,'String',[]);
set(handles.INFOname_Re,'String',[]);
set(handles.INFO_ForK,'String',[]);
set(handles.INFOname_ForK,'String',[]);

set(handles.enablingFile, 'Enable', 'Off');

% Delete data
handles.data = [];
handles.file_name = '';
handles.chord = [];
handles.location = [];
handles.coordinate = [];

% Close all plots
gui_name = get(0, 'CurrentFigure');
set(gui_name, 'HandleVisibility', 'off');
close all;
set(gui_name, 'HandleVisibility', 'on');

set(handles.loc,'BackgroundColor',handles.loc_found);
set(handles.coords,'BackgroundColor',handles.coord_found);

clc

guidata(hObject, handles);



% ------------------------------------
% Plots
% ------------------------------------
        
% --- Executes on button press in plot_profile ('Wing profile').
function plot_profile_Callback(hObject, eventdata, handles)
% hObject    handle to plot_profile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

xs = [handles.coordinate(:,1), handles.coordinate(:,3)]/100;
ys = [handles.coordinate(:,2), handles.coordinate(:,4)]/100;

figure
plot(xs,ys,'-k','LineWidth',2)
axis equal
grid on
xlim([-0.01 1.01])
ylim([-0.15 0.15])
set(gca,'XTick',[0:0.2:1])
set(gca,'YTick',[-0.15:0.05:0.15])
set(gca, 'XTickLabel', [])
set(gca, 'YTickLabel', []);
set(gcf, 'Position', [10 round(handles.sSize(2)*0.75) 500 160])
set(gca,'OuterPosition',[-0.15 -0.15 1.26 1.26])
set(gcf,'Name',handles.file_name(1:end-4))

if get(handles.plot_taps,'Value') == 1
    LE_u = find((handles.location(2:end)-handles.location(1:end-1))<0,1,'last')+1;
    xp = [handles.location', [ones(LE_u,1); -ones(length(handles.location)-LE_u,1)]];
    yp = aerofoil_taps(xs, ys, xp);
    hold on
    plot(xp, yp, 'or')
end



% --- Executes on button press in plot_taps ('Pressure tappings').
function plot_taps_Callback(hObject, eventdata, handles)
% hObject    handle to plot_profile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes on button press in plot_angle ('Angle of incidence').
function plot_angle_Callback(hObject, eventdata, handles)
% hObject    handle to plot_angle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if size(handles.data,1)==handles.data(1,11)+4
    dat_line = 5;
else
    dat_line = 2;
end
        
angle = handles.data(dat_line:end,end);
T = [0:1:size(handles.data,1)-dat_line]*(1/handles.data(1,18));

if get(handles.time_seconds, 'Value')==1
    nonT = T;
    Tlab = 'Time (s)';
elseif get(handles.time_scaled, 'Value')==1
    if handles.data(1,7) == 1 || handles.data(1,7) == 5
        nonT = T * handles.data(1,10) * 2 * pi;
        Tlab = '$\omega t$ (rad)';
    else
        nonT = T * handles.data(1,23) / handles.chord;
        Tlab = 'Non-dimensional time, $tV/c$';
    end
end

hFig = figure;
plot(nonT,angle,'-k','LineWidth',2)
xlim([nonT(1) nonT(end)])
xlabel(Tlab)
ylabel('Angle of attack (deg)')
grid on
set(hFig, 'Position', [10 round(handles.sSize(2)*0.4) 600 300])
set(gca,'OuterPosition',[-0.09 0.06 1.14 0.96])
set(hFig,'Name',handles.file_name)



% --- Executes on button press in plot_wind ('Wind speed').
function plot_wind_Callback(hObject, eventdata, handles)
% hObject    handle to plot_wind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if size(handles.data,1)==handles.data(1,11)+4
    dat_line = 5;
else
    dat_line = 2;
end
        
pdyn = handles.data(dat_line:end,1);
rho = handles.data(1,6)*133.3/(287*(handles.data(1,5)+273.15));

V = sqrt(2*pdyn./rho);
T = [0:1:size(handles.data,1)-dat_line]*(1/handles.data(1,18));

if get(handles.time_seconds, 'Value')==1
    nonT = T;
    Tlab = 'Time (s)';
elseif get(handles.time_scaled, 'Value')==1
    if handles.data(1,7) == 1 || handles.data(1,7) == 5
        nonT = T * handles.data(1,10) * 2 * pi;
        Tlab = '$\omega t$ (rad)';
    else
        nonT = T * handles.data(1,23) / handles.chord;
        Tlab = 'Non-dimensional time, $tV/c$';
    end
end

hFig = figure;
plot(nonT,V,'-k','LineWidth',2)
xlim([nonT(1) nonT(end)])
xlabel(Tlab)
ylabel('Wind tunnel speed (m/s)')
grid on
set(hFig, 'Position', [10 handles.sSize(2)*0.05 600 300])
set(gca,'OuterPosition',[-0.113 0.06 1.165 0.96])
set(hFig,'Name',handles.file_name)



% --- Executes on button press in plot_cp2D ('Pressure coeff. 2D').
function plot_cp2D_Callback(hObject, eventdata, handles)
% hObject    handle to plot_cp2D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if size(handles.data,1)==handles.data(1,11)+4
    dat_line = 5;
else
    dat_line = 2;
end

T = [0:1:size(handles.data,1)-dat_line]*(1/handles.data(1,18));
if get(handles.time_seconds, 'Value')==1
    nonT = T;
    Tlab = 'Time (s)';
elseif get(handles.time_scaled, 'Value')==1
    if handles.data(1,7) == 1 || handles.data(1,7) == 5
        nonT = T * handles.data(1,10) * 2 * pi;
        Tlab = '$\omega t$ (rad)';
    else
        nonT = T * handles.data(1,23) / handles.chord;
        Tlab = 'Non-dimensional time, $tV/c$';
    end
end

if get(handles.side_u, 'Value')==1
    % Plot of the upper surface only
    max_x = find((handles.location(2:end)-handles.location(1:end-1))<0,1,'last')+2;
elseif get(handles.side_ul, 'Value')==1
    % Plot of the upper&lower surfaces
    max_x = size(handles.data,2)-1;
end

hFig = figure;
step = 2;
max_cp = 0;
for i = 2 : max_x
    if i==2
        figure(hFig)
        [ax,h1,h2] = plotyy(nonT, step*(i-1)-handles.data(dat_line:end,i),...
                            [nonT(1) nonT(end)], [1 1]*step*(i-1));
    else
        hold(ax(1),'on')
        plot(ax(1),nonT,step*(i-1)-handles.data(dat_line:end,i),'-k','LineWidth',1)
        hold(ax(2),'on')
        plot(ax(2),[nonT(1) nonT(end)], [1 1]*step*(i-1), '--k', 'LineWidth',1)
        max_cp = max(max_cp, max(step*(i-1)-handles.data(dat_line:end,i)));
    end
        
end

set(h1,'LineStyle','-','LineWidth',1,'Color','k')
set(get(ax(1),'Ylabel'),'String','Pressure coefficient')
set(get(ax(1),'Ylabel'),'FontSize',get(0,'DefaultAxesFontSize'))
set(get(ax(1),'Ylabel'),'FontName',get(0,'DefaultAxesFontName'))
set(ax(1),'ycolor','k') 
ylim(ax(1), [0 1.05*max_cp])
set(ax(1),'YTick',[0:step:6*step])

set(h2,'LineStyle','--','LineWidth',1,'Color','k')
set(ax(2),'ycolor','k') 
ylim(ax(2), [0 1.05*max_cp])
set(ax(2),'YTick',[step:step:step*(max_x-1)])
locationString = num2str(handles.location');
for i = 1 : max_x-1
    if str2double(locationString(i,:))==0
        i_string = '0.';
        while size(i_string,2)<size(locationString,2)
            i_string = strcat(i_string,'0');
        end
        locationString(i,:) = i_string;
    elseif findstr(locationString(i,:),'0.')~=1
        i_string = locationString(i,findstr(locationString(i,:),'0.'):end);
        while size(i_string,2)<size(locationString,2)
            i_string = strcat(i_string,'0');
        end
        locationString(i,:) = i_string;
    end
end
        
set(ax(2),'YTickLabel',strcat('x/c=',locationString))

xlim(ax(2), [nonT(1) nonT(end)])
xlim(ax(1), [nonT(1) nonT(end)])
xlabel(Tlab)
set(hFig, 'Position', [60 60 handles.sSize(1)/2 handles.sSize(2)*0.85])
set(gca,'OuterPosition',[-0.09 -0.04 0.92 1.1])

set(hFig,'Name',handles.file_name)



% --- Executes on button press in plot_cp3D ('Pressure coeff. 3D').
function plot_cp3D_Callback(hObject, eventdata, handles)
% hObject    handle to plot_cp3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if size(handles.data,1)==handles.data(1,11)+4
    dat_line = 5;
else
    dat_line = 2;
end

angle = handles.data(dat_line:end,end);

T = [0:1:size(handles.data,1)-dat_line]*(1/handles.data(1,18));
if get(handles.time_seconds, 'Value')==1
    nonT = T;
    Tlab = 'Time (s)';
elseif get(handles.time_scaled, 'Value')==1
    if handles.data(1,7) == 1 || handles.data(1,7) == 5
        nonT = T * handles.data(1,10) * 2 * pi;
        Tlab = '$\omega t$ (rad)';
    else
        nonT = T * handles.data(1,23) / handles.chord;
        Tlab = 'Non-dimensional time, $tV/c$';
    end
end

max_x = find((handles.location(2:end)-handles.location(1:end-1))<0,1,'last')+2;

figure
set(gcf,'Position',[1 1 handles.sSize(1) handles.sSize(2)])
for i = length(nonT):-1:1
    fill3(nonT(i)*ones(1,max_x-1+2),...
          [handles.location(1:max_x-1), handles.location(max_x-1), handles.location(1)],...
          [-handles.data(i+dat_line-1,2:max_x), 0, 0],'w')
    hold on
end
view(-53,22)

xlim([0 nonT(end)/2*1.2])
set(gca,'ydir','reverse')
h = ylabel('$x/c$');
ylim([0 4])
set(gca,'YTick',[0 1])
zlabel('$-c_p$')
zlim([-2 10])
set(gca,'ZTick',[-2:2:10])
pos = get(h,'pos');
set(h,'pos',pos+[0 -1.9 1.2])

for i = [-2:2:10]
    plot3([0 nonT(end)], [0 0], i*[1 1], ':k')
end
plot3([0 nonT(end)], [1 1], [-2 -2], ':k')
plot3([0 0], [1 1], [-2 0], ':k')
plot3(nonT(end)*[1 1], [1 1], [-2 0], ':k')
xt = get(gca,'XTick');
set(gca,'XTick',[])
for i = 0:(xt(2)-xt(1)):nonT(end)
    plot3(i*[1 1], [1 3.5], [-2 -2], ':k')
    text(i,1.35,-1.8, num2str(round(i*100)/100), 'FontAngle','italic', 'FontSize',get(h,'FontSize'))
    text(i,1.7,-1.8, num2str(round(interp1(nonT,angle,i)*100)/100), 'FontAngle','italic', 'FontSize',get(h,'FontSize'))
end
plot3(nonT(end)*[1 1], [0 3.5], [-2 -2], ':k')
text(nonT(end),1.35,-1.8, num2str(round(nonT(end)*100)/100), 'FontAngle','italic', 'FontSize',get(h,'FontSize'))
text(0,1.1,-2.8, Tlab, 'FontAngle','italic', 'FontSize',get(h,'FontSize'))
text(nonT(end),1.7,-1.8, num2str(round(angle(end)*100)/100), 'FontAngle','italic', 'FontSize',get(h,'FontSize'))
text(0,1.45,-3.1, 'Angle of attack (deg)', 'FontAngle','italic', 'FontSize',get(h,'FontSize'))

set(gcf,'Name',handles.file_name)



% --- Executes on button press in surf_cp ('Surface Cp contour').
function surf_cp_Callback(hObject, eventdata, handles)
% hObject    handle to surf_cp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if size(handles.data,1)==handles.data(1,11)+4
    dat_line = 5;
else
    dat_line = 2;
end

T = [0:1:size(handles.data,1)-dat_line]*(1/handles.data(1,18));
if get(handles.time_seconds, 'Value')==1
    nonT = T;
    Tlab = 'Time (s)';
    Opos = [-0.03 -0.06 1.095 1.12];
elseif get(handles.time_scaled, 'Value')==1
    if handles.data(1,7) == 1 || handles.data(1,7) == 5
        nonT = T * handles.data(1,10) * 2 * pi;
        Tlab = '$\omega t$ (rad)';
        Opos = [-0.07 -0.06 1.14 1.12];
    else
        nonT = T * handles.data(1,23) / handles.chord;
        Tlab = 'Non-dimensional time, $tV/c$';
        Opos = [-0.16 -0.06 1.24 1.12];
    end
end


% Create matrices
max_x = find((handles.location(2:end)-handles.location(1:end-1))<0,1,'last')+2;
x_surf = [1-handles.location(1:max_x-1), 1+handles.location(max_x:end)];
cp_surf = handles.data(dat_line:end, 2:end-1);

hFig = figure;
set(hFig, 'Position', [110 60 handles.sSize(1)/2 handles.sSize(2)*0.85])
contourf(x_surf, nonT, cp_surf, 100, 'LineStyle', 'none')
hold on
plot([1 1], [0 nonT(end)], '--k')
xlabel('Upper surface $\;\;\;\;\;$ - $\;\;\;\;\;$ Lower surface')
ylabel(Tlab)
xlim([0 2])
set(gca,'XTick',[0:0.25:2])
set(gca,'XTickLabel',{'TE', '.75', '.5', '.25', 'LE', '.25', '.5', '.75', 'TE'})
set(gca,'OuterPosition',Opos)
hCol = colorbar;
set(gca,'Clim', [-7 1])
set(hCol,'YTick', [-7:1:1])

set(hFig,'Name',handles.file_name)



% --- Executes on button press in plot_coeffs ('Force coefficients').
function plot_coeffs_Callback(hObject, eventdata, handles)
% hObject    handle to lift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if size(handles.data,1)==handles.data(1,11)+4
    dat_line = 5;
else
    dat_line = 2;
end

alpha = handles.data(dat_line:end,end);

T = [0:1:size(handles.data,1)-dat_line]*(1/handles.data(1,18));
if get(handles.time_seconds, 'Value')==1
    nonT = T;
    Tlab = 'Time (s)';
elseif get(handles.time_scaled, 'Value')==1
    if handles.data(1,7) == 1 || handles.data(1,7) == 5
        nonT = T * handles.data(1,10) * 2 * pi;
        Tlab = '$\omega t$ (rad)';
    else
        nonT = T * handles.data(1,23) / handles.chord;
        Tlab = 'Non-dimensional time, $tV/c$';
    end
end

xs = [handles.coordinate(:,1), handles.coordinate(:,3)]/100;
ys = [handles.coordinate(:,2), handles.coordinate(:,4)]/100;
LE_u = find((handles.location(2:end)-handles.location(1:end-1))<0,1,'last')+1;
xp = [handles.location', [ones(LE_u,1); -ones(length(handles.location)-LE_u,1)]];

[yp, dyp, lp] = aerofoil_taps(xs, ys, xp);

% Lift coefficient calculation
cl = zeros(size(handles.data,1)-dat_line+1,1);
cd = zeros(size(handles.data,1)-dat_line+1,1);
cy = zeros(size(handles.data,1)-dat_line+1,1);
cx = zeros(size(handles.data,1)-dat_line+1,1);
cm_x = zeros(size(handles.data,1)-dat_line+1,1);
% cm_nx = zeros(size(handles.data,1)-dat_line+1,1);

for i = 1 : length(alpha)
    
    angle = atan(dyp(:,1));
    angle = angle+pi/2;
    find_ls = find(dyp(:,2)==-1,1,'first');
    angle(find_ls:end) = pi+angle(find_ls:end,:);
    
    % Data as they are
    cp = (handles.data(i+dat_line-1,2:end-1))';
    % Data applying offset and gain values
%     cp = ((handles.data(i+4,2:end-1)-handles.data(4,2:end-1)) .* handles.data(2,2:end-1) .* handles.data(3,2:end-1))';
    
    cy_i = -cp.*lp.*sin(angle);
    cx_i = -cp.*lp.*cos(angle);
    
    % Forces normal and parallel to the chord
    cy(i) = sum(cy_i);
    cx(i) = -sum(cx_i);
    
    % Forces normal and parallel to the freestream velocity
    cl(i) = sum(cx_i)*sin(alpha(i)/180*pi) + sum(cy_i)*cos(alpha(i)/180*pi);
    cd(i) = -sum(cx_i)*cos(alpha(i)/180*pi) + sum(cy_i)*sin(alpha(i)/180*pi);
    
    % Moment about the quarter-chord with x contribution
    cm_x(i) = sum(-cy_i.*(xp(:,1)-0.25) + cx_i.*yp(:,1));
%     % Moment about the quarter-chord without x contribution
%     cm_nx(i) = sum(-cy_i.*(xp(:,1)-0.25));

%     % Arrows plot Movie
%     figure(1)
%     hold on
%     quiver([0 0], [0 0],...
%            [sum(-cx_i)*cos(-alpha(i)/180*pi) -sum(cy_i)*sin(-alpha(i)/180*pi)], [sum(-cx_i)*sin(-alpha(i)/180*pi) sum(cy_i)*cos(-alpha(i)/180*pi)], '-b', 'LineWidth', 2)
%     quiver([0 0], [0 0],...
%            [cd(i) 0], [0 cl(i)], '-g', 'LineWidth', 2)
%     quiver((xp(:,1)-0.25)*cos(alpha(i)/180*pi), (xp(:,1)-0.25)*sin(-alpha(i)/180*pi)+yp(:,1)*cos(alpha(i)/180*pi),...
%            cx_i*cos(-alpha(i)/180*pi)-cy_i*sin(-alpha(i)/180*pi), cx_i*sin(-alpha(i)/180*pi)+cy_i*cos(-alpha(i)/180*pi), 'k')
%     plot((xp(:,1)-0.25)*cos(alpha(i)/180*pi), (xp(:,1)-0.25)*sin(-alpha(i)/180*pi)+yp(:,1)*cos(alpha(i)/180*pi), '-or')
%     plot(([0 1]-0.25)*cos(alpha(i)/180*pi), ([0 1]-0.25)*sin(-alpha(i)/180*pi), '-k')
%     axis equal
%     xlim([-0.5 1])
%     legend('c_n, c_t', 'c_l, c_d', 'c_p', 'taps', 'chord')
%     pos = get(gcf,'Position');
%     set(gcf, 'Position', [handles.sSize(1)/2 handles.sSize(2)/2 pos(3:4)])
%     M(i) = getframe(1);
%     clf(1)
    
end


if get(handles.c_nmt, 'Value')==1
    % Plot coefficients in chord reference
    Cv = cy;
    Ch = cx;
    Cm = cm_x;
    ylab_v = '$c_n$';
    ylab_h = '$c_t$';
    ylab_m = '$c_{m_{1/4}}$';
elseif get(handles.c_lmd, 'Value')==1
    % Plot coefficients in wind reference
    Cv = cl;
    Ch = cd;
    Cm = cm_x;
    ylab_v = '$c_l$';
    ylab_h = '$c_d$';
    ylab_m = '$c_{m_{1/4}}$';
end
    

hFig1 = figure;
% Plotting the lift coefficient against the angle of attack
subplot(3,2,2)
plot(alpha, Cv, '-k','LineWidth',2)
xlabel('Angle of attack (deg)')
if strcmp(handles.file_name(3),'4') == 0
    alpha_lim = [min(-20,min(alpha)*1.05) max(40,max(alpha)*1.05)];
    set(gca,'XTick',[-40:5:40])
    set(gca,'XTickLabel',{'-40' '' '-30' '' '-20' '' '-10' '' '0' '' '10' '' '20' '' '30' '' '40'})
else
    alpha_lim = [min(alpha)*0.98 max(alpha)*1.02];
end
axis([alpha_lim min(min(Cv),-2) max(max(Cv),2)])
ylabel(ylab_v)
set(gca,'YTick',[-10:0.5:10])
hold on
plot(get(gca,'XLim'), [0 0], '-k')
plot([0 0], get(gca,'YLim'), '-k')

% Plotting the lift coefficient against the time
subplot(3,2,1)
plot(nonT, Cv, '-k','LineWidth',2)
axis([nonT(1) nonT(end) min(min(Cv),-2) max(max(Cv),2)])
xlabel(Tlab)
ylabel(ylab_v)
hold on
plot(get(gca,'XLim'), [0 0], '-k')
set(gca,'YTick',[-10:0.5:10])

% Plotting the drag coefficient against the angle of attack
subplot(3,2,6)
plot(alpha, Ch, '-k','LineWidth',2)
xlabel('Angle of attack (deg)')
if strcmp(handles.file_name(3),'4') == 0
    alpha_lim = [min(-20,min(alpha)*1.05) max(40,max(alpha)*1.05)];
    set(gca,'XTick',[-40:5:40])
    set(gca,'XTickLabel',{'-40' '' '-30' '' '-20' '' '-10' '' '0' '' '10' '' '20' '' '30' '' '40'})
else
    alpha_lim = [min(alpha)*0.98 max(alpha)*1.02];
end
axis([alpha_lim min(min(Ch),-0.3) max(max(Ch),0.2)])
ylabel(ylab_h)
set(gca,'YTick',[-10:0.1:10])
hold on
plot(get(gca,'XLim'), [0 0], '-k')
plot([0 0], get(gca,'YLim'), '-k')

% Plotting the drag coefficient against the time
subplot(3,2,5)
plot(nonT, Ch, '-k','LineWidth',2)
axis([nonT(1) nonT(end) min(min(Ch),-0.3) max(max(Ch),0.2)])
xlabel(Tlab)
ylabel(ylab_h)
hold on
plot(get(gca,'XLim'), [0 0], '-k')
set(gca,'YTick',[-10:0.1:10])

% Plotting the moment coefficient against the angle of attack
subplot(3,2,4)
plot(alpha, Cm, '-k','LineWidth',2)
xlabel('Angle of attack (deg)')
if strcmp(handles.file_name(3),'4') == 0
    alpha_lim = [min(-20,min(alpha)*1.05) max(40,max(alpha)*1.05)];
    set(gca,'XTick',[-40:5:40])
    set(gca,'XTickLabel',{'-40' '' '-30' '' '-20' '' '-10' '' '0' '' '10' '' '20' '' '30' '' '40'})
else
    alpha_lim = [min(alpha)*0.98 max(alpha)*1.02];
end
axis([alpha_lim min(min(Cm),-0.5) max(max(Cm),0.1)])
ylabel(ylab_m)
set(gca,'YTick',[-10:0.1:10])
hold on
plot(get(gca,'XLim'), [0 0], '-k')
plot([0 0], get(gca,'YLim'), '-k')

% Plotting the moment coefficient against the time
subplot(3,2,3)
plot(nonT, Cm, '-k','LineWidth',2)
axis([nonT(1) nonT(end) min(min(Cm),-0.5) max(max(Cm),0.1)])
xlabel(Tlab)
ylabel(ylab_m)
hold on
plot(get(gca,'XLim'), [0 0], '-k')
set(gca,'YTick',[-10:0.1:10])

set(hFig1,'Name',handles.file_name)
set(hFig1,'Position',[1 1 handles.sSize(1) handles.sSize(2)])



% --- Executes on button press in plot_close_all ('Close all plots').
function plot_close_all_Callback(hObject, eventdata, handles)
% hObject    handle to plot_close_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

gui_name = get(0, 'CurrentFigure');
set(gui_name, 'HandleVisibility', 'off');
close all;
% clc
set(gui_name, 'HandleVisibility', 'on');

set(handles.loc,'BackgroundColor',handles.loc_found);
set(handles.coords,'BackgroundColor',handles.coord_found);



% -------------------------------------------------------------------------
% CASE COUNTER
% -------------------------------------------------------------------------

% --- Executes on selection change in list_models.
function list_models_Callback(hObject, eventdata, handles)
% hObject    handle to list_models (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns list_models contents as cell array
%        contents{get(hObject,'Value')} returns selected item from list_models

model = get(hObject, 'Value')-1;

data_name = strcat(handles.folder_name,handles.slash,'*.dat');
d = dir(data_name);
N = size(d,1);

if sum(model==0)==1
    count = ones(1,N);
else
    ab_name = zeros(1,N);
    for i = 1 : N
        ab_name(i) = str2num(d(i).name(1:2));
    end
    
    count = zeros(1,N);
    for i = 1 : length(model)
        count = count + (ab_name==model(i));
    end
end

handles.count_ab = count;
set(handles.cases_no,'String',num2str(sum(handles.count_ab.*handles.count_c.*handles.count_d)));


guidata(hObject, handles);
    


% --- Executes during object creation, after setting all properties.
function list_models_CreateFcn(hObject, eventdata, handles)
% hObject    handle to list_models (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in list_experiments.
function list_experiments_Callback(hObject, eventdata, handles)
% hObject    handle to list_experiments (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns list_experiments contents as cell array
%        contents{get(hObject,'Value')} returns selected item from list_experiments

exp = get(hObject, 'Value')-2;

data_name = strcat(handles.folder_name,handles.slash,'*.dat');
d = dir(data_name);
N = size(d,1);

if sum(exp==-1)==1
    count = ones(1,N);
else
    c_name = zeros(1,N);
    for i = 1 : N
        c_name(i) = str2num(d(i).name(3));
    end
    
    count = zeros(1,N);
    for i = 1 : length(exp)
        count = count + (c_name==exp(i));
    end
end

handles.count_c = count;
set(handles.cases_no,'String',num2str(sum(handles.count_ab.*handles.count_c.*handles.count_d)));


guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function list_experiments_CreateFcn(hObject, eventdata, handles)
% hObject    handle to list_experiments (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in list_motions.
function list_motions_Callback(hObject, eventdata, handles)
% hObject    handle to list_motions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns list_motions contents as cell array
%        contents{get(hObject,'Value')} returns selected item from list_motions

motion = get(hObject, 'Value')-2;

data_name = strcat(handles.folder_name,handles.slash,'*.dat');
d = dir(data_name);
N = size(d,1);

if sum(motion==-1)==1
    count = ones(1,N);
else
    d_name = zeros(1,N);
    for i = 1 : N
        d_name(i) = str2num(d(i).name(4));
    end
    
    count = zeros(1,N);
    for i = 1 : length(motion)
        count = count + (d_name==motion(i));
    end
end

handles.count_d = count;
set(handles.cases_no,'String',num2str(sum(handles.count_ab.*handles.count_c.*handles.count_d)));


guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function list_motions_CreateFcn(hObject, eventdata, handles)
% hObject    handle to list_motions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in selectedfiles_select ('Select data file from selected cases').
function selectedfiles_select_Callback(hObject, eventdata, handles)
% hObject    handle to selectedfiles_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

selected_files = handles.count_ab.*handles.count_c.*handles.count_d;

if sum(selected_files) > 0
    
    % Preparing RIBs information
    infos = strcat( num2str(strcmp(get(handles.selected_freestream,'Enable'),'on')),...
                    num2str(get(handles.selected_freestream,'Value')),...
                    num2str(get(handles.selected_motion,'Value')) );
    if strcmp(get(handles.selected_freestream,'Enable'),'on')*(get(handles.selected_freestream,'Value')+get(handles.selected_motion,'Value')) >= 1
        slash = findstr(handles.folder_name,handles.slash);
        filename = strcat(handles.folder_name(1:slash(end)),get(handles.RIBs_filename,'String'));
        load(filename)
    else
        RIBs_info = [];
    end
    
    outputarg = SelectedFiles(handles.folder_name, selected_files, infos, RIBs_info);

    if strcmp(outputarg,'no_file') ~= 1
        handles.file_name = outputarg;

        set(handles.datafile_typeWARNING,'String',[]);
        set(handles.datafile_view,'String',handles.file_name(1:end-4));

        % Finding the list number
        d = dir(strcat(handles.folder_name,handles.slash,'*.dat'));
        j = 1;
        found = 0;
        while found == 0
            if strcmp(d(j).name, handles.file_name) == 1
                found = 1;
                set(handles.datafile_typenumber,'String',num2str(j));
                set(handles.slider_typenumber,'Value',j);
            end
            j = j+1;
        end

        handles.data = importdata(strcat(handles.folder_name,handles.slash,handles.file_name));
        set(handles.enablingFile, 'Enable', 'On');
        if size(handles.data,1)~=handles.data(1,11)+4
            set(handles.plot_calib, 'Enable', 'Off');
        end

        % INFO boxes
        models = get(handles.list_models,'String');
        set(handles.INFO_model,'String',models{str2double(handles.file_name(1:2))+1}(4:end));
        set(handles.INFOname_model,'String','Model');
        experiment = get(handles.list_experiments,'String');
        set(handles.INFO_experiment,'String',experiment{str2double(handles.file_name(3))+2}(3:end));
        set(handles.INFOname_experiment,'String','Type of experiment');
        motion = get(handles.list_motions,'String');
        set(handles.INFO_motion,'String',motion{str2double(handles.file_name(4))+2}(3:end));
        set(handles.INFOname_motion,'String','Motion type');
        set(handles.INFO_Re,'String',handles.data(1,20));
        set(handles.INFOname_Re,'String','Reynolds number');
        if handles.data(1,7) == 1 || handles.data(1,7) == 5 || handles.data(1,7) == 7
            set(handles.INFO_ForK,'String',handles.data(1,22));
            set(handles.INFOname_ForK,'String','Reduced frequency');
        elseif handles.data(1,7) == 2 || handles.data(1,7) == 3
            set(handles.INFO_ForK,'String',handles.data(1,22));
            set(handles.INFOname_ForK,'String','Reduced pitch-rate');
        else
            set(handles.INFO_ForK,'String',[]);
            set(handles.INFOname_ForK,'String',[]);
        end
    
        % Definition of the chord of the model
        if handles.data(1,30) == 12
            handles.chord = 0.55/2;
        else
            handles.chord = 0.55;
        end

        % Loading pressure transducer locations
        FID = fopen(handles.locfile,'r');
        if FID > 0
            locations = textscan(FID, '%s', 'Delimiter', '\n');
            handles.location = str2num(locations{1}{handles.data(1,30)+1});
        end

        % Loading aerofoil coordinates
        handles.coordfile
        coordinates = xlsread(handles.coordfile);
        handles.coordinate = coordinates(:,handles.data(1,30)*4-3:handles.data(1,30)*4);
        handles.coordinate = handles.coordinate(isfinite(handles.coordinate(:,1)),:);

        fclose all;
    end

end

guidata(hObject, handles);


% --- Executes on button press in selected_freestream.
function selected_freestream_Callback(hObject, eventdata, handles)
% hObject    handle to selected_freestream (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of selected_freestream


% --- Executes on button press in selected_motion.
function selected_motion_Callback(hObject, eventdata, handles)
% hObject    handle to selected_motion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of selected_motion
