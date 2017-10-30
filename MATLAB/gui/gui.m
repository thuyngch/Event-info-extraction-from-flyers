%------------------------------------------------------------------------------
% Date    : Oct 28, 2017
% Description :
%   This file is the GUI of system that uses the improved algorithm.
%------------------------------------------------------------------------------



function varargout = gui(varargin)
% GUI MATLAB code for gui.fig
%      GUI, by itself, creates a new GUI or raises the existing
%      singleton*.
%
%      H = GUI returns the handle to a new GUI or the handle to
%      the existing singleton*.
%
%      GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI.M with the given input arguments.
%
%      GUI('Property','Value',...) creates a new GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui

% Last Modified by GUIDE v2.5 28-Oct-2017 21:11:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_OutputFcn, ...
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


% --- Executes just before gui is made visible.
function gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui (see VARARGIN)

% Choose default command line output for gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Clean command window and axis
clc
fprintf('>>> Start GUI\n')
set(handles.txtTime, 'String', '');
warning off;

% UIWAIT makes gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in btnSelectFile.
function btnSelectFile_Callback(hObject, eventdata, handles)
% hObject    handle to btnSelectFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Open window explorer to select an image file
[handles.filename, handles.filepath] = uigetfile('*.jpg', ...
	'Select an image file...', './im/capture/');
handles.fileim = [handles.filepath, handles.filename];
guidata(hObject, handles);
fprintf('\n>>> Selected image: %s\n', handles.fileim)

% Close all of figures
h = findobj('type','figure');
for i = 1: length(h)
    if strcmp(h(i).Name, 'gui'); continue; end
    clf(h(i).Number)
end


% --- Executes on button press in cboxObj.
function cboxObj_Callback(hObject, eventdata, handles)
% hObject    handle to cboxObj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cboxObj


% --- Executes on button press in btnExtract.
function btnExtract_Callback(hObject, eventdata, handles)
% hObject    handle to btnExtract (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Remove previous sub-images and text file
fprintf('\n>>> Start extraction...\n')
cd result
delete(sprintf('%c.txt', handles.filename(1:end-4)));
delete *.jpg
cd ../

% Load database
load dt_set.mat
load ocr.mat

% Checkbox of plotting
bound = get(handles.cboxBound, 'Value');
homo = get(handles.cboxHomo, 'Value');
obj = get(handles.cboxObj, 'Value');
flgPlot = [bound, homo, obj];

% Extract
output_txt_path = ['result/', handles.filename(1:end-4), '.txt'];
[str, time] = recognizeText(handles.fileim, output_txt_path, ocr_name, dtset, flgPlot);
guidata(hObject, handles);
str_time = timeFilter(str, dtset.time);
str_loc = locFilter(str, dtset.location);

% Display result to the static box
set(handles.txtTime, 'String', str_time);
set(handles.txtLocation, 'String', str_loc);
fprintf('\nDone!\n')
fprintf('Homography interval     : %f [s]\n', time.homo)
fprintf('Text detection interval : %f [s]\n', time.region)
fprintf('Tesseract interval      : %f [s]\n', time.tesseract)
fprintf('Total interval          : %f [s]\n\n',...
    time.homo+time.region+time.tesseract)

% --- Executes on button press in cboxHomo.
function cboxHomo_Callback(hObject, eventdata, handles)
% hObject    handle to cboxHomo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cboxHomo


% --- Executes on button press in cboxBound.
function cboxBound_Callback(hObject, eventdata, handles)
% hObject    handle to cboxBound (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cboxBound
