function varargout = top3dGUI(varargin)
% TOP3DGUI MATLAB code for top3dGUI.fig
%      TOP3DGUI, by itself, creates a new TOP3DGUI or raises the existing
%      singleton*.
%
%      H = TOP3DGUI returns the handle to a new TOP3DGUI or the handle to
%      the existing singleton*.
%
%      TOP3DGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TOP3DGUI.M with the given input arguments.
%
%      TOP3DGUI('Property','Value',...) creates a new TOP3DGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before top3dGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to top3dGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help top3dGUI

% Last Modified by GUIDE v2.5 08-Apr-2015 14:59:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @top3dGUI_OpeningFcn, ...
    'gui_OutputFcn',  @top3dGUI_OutputFcn, ...
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


% --- Executes just before top3dGUI is made visible.
function top3dGUI_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to top3dGUI (see VARARGIN)

% Choose default command line output for top3dGUI
handles.output = hObject;

% Welcom message
fprintf(1, ' Top3d Graphical User Interface (GUI)\n');
fprintf(1, ' An efficient 3D Topology Optimization program\n');
fprintf(1, ' written by K. Liu and A. Tovar\n');
fprintf(1, ' Indiana Univ-Purdue Univ Indianapolis, U.S.A\n');
fprintf(1, ' visit us at <a href="http://www.top3dapp.com">Top3dAPP.com</a>\n');
% Set unit
% set(hObject,'Unit','normalized')
% Warning
warning('off', 'MATLAB:hg:patch:CannotUseFaceVertexCDataOfSize0');
% Apply default values
set_DesignPanel(handles);
get_DesignPanel(hObject, eventdata, handles, 'initial');
handles = guidata(hObject);
% - BC
get_BCPanel(hObject, eventdata, handles, 'initial');
handles = guidata(hObject);
% - LC
get_LCPanel(hObject, eventdata, handles, 'initial');
handles = guidata(hObject);
% Initial domains
x = 0.5*ones(5,10,2);
axes(handles.axes_design);
plot_3d(x, 1, 0.5);
% Plot arrow
sp = [0 0 5];
epx = sp + [2 0 0];
epy = sp - [0 0 2];
epz = sp + [0 2 0];
ar = plot_arrow(sp, epx);
plot_arrow(ar,'Width',2);
ar = plot_arrow(sp, epy);
plot_arrow(ar,'Width',2);
ar = plot_arrow(sp, epz);
plot_arrow(ar,'Width',2);
xlabel('x'), ylabel('z'), zlabel('y')

set(handles.axes_res,'Visible','off');
set(handles.uipanel15,'Visible','off');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes top3dGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = top3dGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in btn_outer.
function btn_outer_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
% hObject    handle to btn_outer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.axes_design);
if isfield(handles,'patch') && ~isempty(handles.patch)
    if isfield(handles.patch,'design')
        delete(handles.patch.design)
    end
end

if ~isempty(handles.bc.fixeddof)
    btn_bc_reset_Callback(hObject, eventdata, handles);
    handles = guidata(hObject);
end

if ~isempty(handles.lc.loaddof) || ~isempty(handles.lc.loadscale)
    btn_lc_reset_Callback(hObject, eventdata, handles);
    handles = guidata(hObject);
end

if handles.optpar.randinit
    x = 0.001+(1-0.001)*rand([handles.domain.nely,handles.domain.nelx,handles.domain.nelz]);
else
    x = repmat(handles.domain.volfrac,[handles.domain.nely,handles.domain.nelx,handles.domain.nelz]);
end

p_design = plot_3d(x,1e-3,0.5);

handles.patch.design = p_design;

guidata(hObject,handles)

function gui_volfrac_Callback(hObject, eventdata, handles)
% hObject    handle to gui_volfrac (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gui_volfrac as text
%        str2double(get(hObject,'String')) returns contents of gui_volfrac as a double
volfrac = str2double(get(hObject, 'String'));
if isnan(volfrac)
    set(hObject, 'String', 0.3);
    errordlg('Input must be a number','Error');
    handles.domain.volfrac = 0.3;
elseif volfrac < 0 || volfrac > 1
    set(hObject, 'String', 0.3);
    errordlg('Input must between 0 and 1','Error');
    handles.domain.volfrac = 0.3;
elseif volfrac == 0
    set(hObject, 'String', 1e-3);
    handles.domain.volfrac = 1e-3;
else
    handles.domain.volfrac = volfrac;
end

guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function gui_volfrac_CreateFcn(hObject, eventdata, handles) %#ok<*INUSD>
% hObject    handle to gui_volfrac (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gui_nelx_Callback(hObject, eventdata, handles)
% hObject    handle to gui_nelx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gui_nelx as text
%        str2double(get(hObject,'String')) returns contents of gui_nelx as a double
nelx = str2double(get(hObject, 'String'));
if isnan(nelx)
    set(hObject, 'String', 30);
    errordlg('Input must be a number','Error');
    handles.domain.nelx = 30;
elseif mod(nelx,1);
    set(hObject, 'String', floor(nelx));
    errordlg('Input must be a real integer','Error');
    handles.domain.nelx = floor(nelx);
else
    handles.domain.nelx = nelx;
end

% Save the new nelx value
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function gui_nelx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gui_nelx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gui_nely_Callback(hObject, eventdata, handles)
% hObject    handle to gui_nely (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gui_nely as text
%        str2double(get(hObject,'String')) returns contents of gui_nely as a double
nely = str2double(get(hObject, 'String'));
if isnan(nely)
    set(hObject, 'String', 10);
    errordlg('Input must be a number','Error');
    handles.domain.nely = 10;
elseif mod(nely,1);
    set(hObject, 'String', floor(nely));
    errordlg('Input must be a real integer','Error');
    handles.domain.nely = floor(nely);
else
    handles.domain.nely = nely;
end

% Save the new nelx value
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function gui_nely_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gui_nely (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gui_nelz_Callback(hObject, eventdata, handles)
% hObject    handle to gui_nelz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gui_nelz as text
%        str2double(get(hObject,'String')) returns contents of gui_nelz as a double
nelz = str2double(get(hObject, 'String'));
if isnan(nelz)
    set(hObject, 'String', 2);
    errordlg('Input must be a number','Error');
    handles.domain.nelz = 2;
elseif mod(nelz,1);
    set(hObject, 'String', floor(nelz));
    errordlg('Input must be a real integer','Error');
    handles.domain.nelz = floor(nelz);
else
    handles.domain.nelz = nelz;
end

% Save the new nelx value
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function gui_nelz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gui_nelz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in chk_ux.
function chk_ux_Callback(hObject, eventdata, handles)
% hObject    handle to chk_ux (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chk_ux
if (get(hObject,'Value') == get(hObject,'Max'))
    % Checkbox is checked-take appropriate action
    % Save the new nelx value
    handles.bc.ux = 1;
else
    % Checkbox is not checked-take appropriate action
    handles.bc.ux = 0;
end
guidata(hObject,handles)

% --- Executes on button press in chk_uy.
function chk_uy_Callback(hObject, eventdata, handles)
% hObject    handle to chk_uy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chk_uy
if (get(hObject,'Value') == get(hObject,'Max'))
    % Checkbox is checked-take appropriate action
    handles.bc.uy = 1;
else
    % Checkbox is not checked-take appropriate action
    handles.bc.uy = 0;
end
guidata(hObject,handles)

% --- Executes on button press in chk_uz.
function chk_uz_Callback(hObject, eventdata, handles)
% hObject    handle to chk_uz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chk_uz
if (get(hObject,'Value') == get(hObject,'Max'))
    % Checkbox is checked-take appropriate action
    % Save the new nelx value
    handles.bc.uz = 1;
else
    % Checkbox is not checked-take appropriate action
    handles.bc.uz = 0;
end
guidata(hObject,handles)


function bc_x1_Callback(hObject, eventdata, handles)
% hObject    handle to bc_x1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bc_x1 as text
%        str2double(get(hObject,'String')) returns contents of bc_x1 as a double
x1 = str2double(get(hObject, 'String'));
if isnan(x1)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
    handles.bc.x1 = 0;
elseif mod(x1,1);
    set(hObject, 'String', floor(x1));
    errordlg('Input must be a real integer','Error');
    handles.bc.x1 = floor(x1);
elseif x1 > handles.domain.nelx;
    set(hObject, 'String', handles.domain.nelx);
    errordlg('Input must less than nelx','Error');
    handles.bc.x1 = handles.domain.nelx;
else
    handles.bc.x1 = x1;
end
% Save the new nelx value
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function bc_x1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bc_x1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bc_y1_Callback(hObject, eventdata, handles)
% hObject    handle to bc_y1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bc_y1 as text
%        str2double(get(hObject,'String')) returns contents of bc_y1 as a double
y1 = str2double(get(hObject, 'String'));
if isnan(y1)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
    handles.bc.y1 = handles.domain.nely;
elseif mod(y1,1);
    set(hObject, 'String', floor(y1));
    errordlg('Input must be a real integer','Error');
    handles.bc.y1 = floor(y1);
elseif y1 > handles.domain.nely;
    set(hObject, 'String', handles.domain.nely);
    errordlg('Input must less than nely','Error');
    handles.bc.y1 = handles.domain.nely;
else
    handles.bc.y1 = y1;
end

% Save the new nelx value
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function bc_y1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bc_y1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bc_z1_Callback(hObject, eventdata, handles)
% hObject    handle to bc_z1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bc_z1 as text
%        str2double(get(hObject,'String')) returns contents of bc_z1 as a double
z1 = str2double(get(hObject, 'String'));
if isnan(z1)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
    handles.bc.z1 = 0;
elseif mod(z1,1);
    set(hObject, 'String', floor(z1));
    errordlg('Input must be a real integer','Error');
    handles.bc.z1 = floor(z1);
elseif z1 > handles.domain.nelz;
    set(hObject, 'String', handles.domain.nelz);
    errordlg('Input must less than nelz','Error');
    handles.bc.z1 = handles.domain.nelz;
else
    handles.bc.z1 = z1;
end
% Save the new nelx value

guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function bc_z1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bc_z1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bc_x2_Callback(hObject, eventdata, handles)
% hObject    handle to bc_x2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bc_x2 as text
%        str2double(get(hObject,'String')) returns contents of bc_x2 as a double
x2 = str2double(get(hObject, 'String'));
if isnan(x2)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
    handles.bc.x2 = 0;
elseif mod(x2,1);
    set(hObject, 'String', floor(x2));
    errordlg('Input must be a real integer','Error');
    handles.bc.x2 = floor(x2);
elseif x2 > handles.domain.nelx;
    set(hObject, 'String', handles.domain.nelx);
    errordlg('Input must less than nelx','Error');
    handles.bc.x2 = handles.domain.nelx;
else
    handles.bc.x2 = x2;
end
% Save the new nelx value
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function bc_x2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bc_x2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bc_y2_Callback(hObject, eventdata, handles)
% hObject    handle to bc_y2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bc_y2 as text
%        str2double(get(hObject,'String')) returns contents of bc_y2 as a double
y2 = str2double(get(hObject, 'String'));
if isnan(y2)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
    handles.bc.y2 = 0;
elseif mod(y2,1);
    set(hObject, 'String', floor(y2));
    errordlg('Input must be a real integer','Error');
    handles.bc.y2 = floor(y2);
elseif y2 > handles.domain.nely;
    set(hObject, 'String', handles.domain.nely);
    errordlg('Input must less than nely','Error');
    handles.bc.y2 = handles.domain.nely;
else
    handles.bc.y2 = y2;
end
% Save the new nelx value
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function bc_y2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bc_y2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bc_z2_Callback(hObject, eventdata, handles)
% hObject    handle to bc_z2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bc_z2 as text
%        str2double(get(hObject,'String')) returns contents of bc_z2 as a double
z2 = str2double(get(hObject, 'String'));
if isnan(z2)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
elseif mod(z2,1);
    set(hObject, 'String', floor(z2));
    errordlg('Input must be a real integer','Error');
    handles.bc.z2 = floor(z2);
elseif z2 > handles.domain.nelz;
    set(hObject, 'String', handles.domain.nelz);
    errordlg('Input must less than nelz','Error');
    handles.bc.z2 = handles.domain.nelz;
else
    handles.bc.z2 = z2;
end

% Save the new nelx value
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function bc_z2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bc_z2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_bc.
function btn_bc_Callback(hObject, eventdata, handles)
% hObject    handle to btn_bc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

nelx = handles.domain.nelx;
nely = handles.domain.nely;
% nelz = handles.domain.nelz;
axes(handles.axes_design);

% USER-DEFINED LOAD DOFs
switch handles.bc.type
    case 'individual'
        ib = handles.bc.x1; jb = handles.bc.y1; kb = handles.bc.z1;  % Coordinates
        fixednid = kb*(nelx+1)*(nely+1)+ib*(nely+1)+(nely+1-jb); % Node IDs
        
        fixeddof = [];
        htri_x = []; htri_y = []; htri_z = [];
        if handles.bc.ux
            fixeddof(:,end+1) = 3*fixednid(:) - 2;
            htri_x = plot_triangle(ib,kb,handles.domain.nely-jb,1,1);
        end
        
        if handles.bc.uy
            fixeddof(:,end+1) = 3*fixednid(:) - 1;
            htri_y = plot_triangle(ib,handles.domain.nely-jb,kb,2,2);
        end
        
        if handles.bc.uz
            fixeddof(:,end+1) = 3*fixednid(:);
            htri_z = plot_triangle(ib,handles.domain.nely-jb,kb,3,3);
        end
        
    case 'distribute'
        ib = handles.bc.x1:handles.bc.x2;
        jb = handles.bc.y1:handles.bc.y2;
        kb = handles.bc.z1:handles.bc.z2;  % Coordinates
        [ib,jb,kb] = meshgrid(ib,jb,kb);
        
        fixednid = kb(:)*(nelx+1)*(nely+1)+ib(:)*(nely+1)+(nely+1-jb(:));  % Node IDs
        fixeddof = [];
        htri_x = []; htri_y = []; htri_z = [];
        if handles.bc.ux
            fixeddof(:,end+1)   = 3*fixednid(:) - 2;
            htri_x = plot_triangle(ib,kb,handles.domain.nely-jb,1,1);
        end
        
        if handles.bc.uy
            fixeddof(:,end+1) = 3*fixednid(:) - 1;
            htri_y = plot_triangle(ib,handles.domain.nely-jb,kb,2,2);
        end
        
        if handles.bc.uz
            fixeddof(:,end+1) = 3*fixednid(:);
            htri_z = plot_triangle(ib,handles.domain.nely-jb,kb,3,3);
        end
end

fixeddof = union(handles.bc.fixeddof, fixeddof);
handles.bc.fixeddof = sort(fixeddof(:));
handles.bc.htri = [htri_x; htri_y; htri_z];
% Save the new nelx value
guidata(hObject,handles)

% il = handles.bc.x1; jl = handles.bc.y1; kl = handles.bc.z1;  % Coordinates
% loadnid = kl*(nelx+1)*(nely+1)+il*(nely+1)+(nely+1-jl); % Node IDs
% loaddof = 3*loadnid(:) - 1;                             % DOFs

% --- Executes when selected object is changed in support_type.
function support_type_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in support_type
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'btn_bc1'
        % Code for when btn_1 is selected.
        handles.bc.type = 'individual';
    case 'btn_bc2'
        % Code for when btn_2 is selected.
        handles.bc.type = 'distribute';
end
guidata(hObject,handles)


% --- Executes on button press in btn_lc.
function btn_lc_Callback(hObject, eventdata, handles)
% hObject    handle to btn_lc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

nelx = handles.domain.nelx;
nely = handles.domain.nely;
% nelz = handles.domain.nelz;
axes(handles.axes_design);
% USER-DEFINED LOAD DOFs

switch handles.lc.type
    case 'individual'
        il = handles.lc.x1; jl = handles.lc.y1; kl = handles.lc.z1;  % Coordinates
        loadnid = kl*(nelx+1)*(nely+1)+il*(nely+1)+(nely+1-jl);      % Node IDs
        % Arrow start point
        sp = [handles.lc.x1 handles.lc.z1 handles.domain.nely-handles.lc.y1];
        
        loaddof = [];
        loadscale = [];
        harr_x = []; harr_y = []; harr_z = [];
        if handles.lc.fx
            loaddof(:,end+1) = 3*loadnid(:) - 2;
            loadscale(:,end+1) = handles.lc.fx;
            % Plot arrow
            ep = sp + sign(handles.lc.fx)*[max(5,floor(handles.domain.nelx/4)) 0 0];
            ar = plot_arrow(sp,ep);
            harr_x = plot_arrow(ar,'Width',2);
        end
        
        if handles.lc.fy
            loaddof(:,end+1) = 3*loadnid(:) - 1;
            loadscale(:,end+1) = handles.lc.fy;
            % Plot arrow
            ep = sp - sign(handles.lc.fy)*[0 0 max(5,floor(handles.domain.nely/4))];
            ar = plot_arrow(sp,ep);
            harr_y = plot_arrow(ar,'Width',2);
        end
        
        if handles.lc.fz
            loaddof(:,end+1) = 3*loadnid(:);
            loadscale(:,end+1) = handles.lc.fz;
            ep = sp + sign(handles.lc.fz)*[0 max(5,floor(handles.domain.nelz/4)) 0];
            % Plot arrow
            ar = plot_arrow(sp,ep);
            harr_z = plot_arrow(ar,'Width',2);
        end
        
    case 'distribute'
        il = handles.lc.x1:handles.lc.x2;
        jl = handles.lc.y1:handles.lc.y2;
        kl = handles.lc.z1:handles.lc.z2;  % Coordinates
        [il,jl,kl] = meshgrid(il,jl,kl);
        
        loadnid = kl(:)*(nelx+1)*(nely+1)+il(:)*(nely+1)+(nely+1-jl(:));  % Node IDs
        loaddof = [];
        loadscale = [];
        harr_x = []; harr_y = []; harr_z = [];
        if handles.lc.fx
            loaddof(:,end+1)   = 3*loadnid(:) - 2;
            loadscale(:,end+1) = repmat(handles.lc.fx,size(loadnid(:)));
        end
        
        if handles.lc.fy
            loaddof(:,end+1) = 3*loadnid(:) - 1;
            loadscale(:,end+1) = repmat(handles.lc.fy,size(loadnid(:)));
        end
        
        if handles.lc.fz
            loaddof(:,end+1) = 3*loadnid(:);
            loadscale(:,end+1) = repmat(handles.lc.fz,size(loadnid(:)));
        end
        
        % Plot arrow
        %         [lds_x,lds_z,lds_y] = meshgrid(il,kl,handles.domain.nely-jl);
        sp = horzcat(il(:),kl(:),handles.domain.nely-jl(:));
        
        if handles.lc.fx
            ep = sp + sign(handles.lc.fx)*...
                repmat([max(5,floor(handles.domain.nelx/4)) 0 0],size(sp,1),1);
            ar = plot_arrow(sp,ep);
            harr_x = plot_arrow(ar,'Width',2);
        end
        
        if handles.lc.fy
            ep = sp - sign(handles.lc.fy)*...
                repmat([0 0 max(5,floor(handles.domain.nely/4))],size(sp,1),1);
            ar = plot_arrow(sp,ep);
            harr_y = plot_arrow(ar,'Width',2);
        end
        
        if handles.lc.fz
            ep = sp + sign(handles.lc.fz)*...
                repmat([0 max(5,floor(handles.domain.nelz/4)) 0],size(sp,1),1);
            ar = plot_arrow(sp,ep);
            harr_z = plot_arrow(ar,'Width',2);
        end
        
end

%
% loaddof = union(handles.lc.loaddof, loaddof);
handles.lc.loaddof = [handles.lc.loaddof; loaddof(:)];
handles.lc.loadscale = [handles.lc.loadscale; loadscale(:)];
handles.lc.harr = [handles.lc.harr; harr_x(:); harr_y(:); harr_z(:)];
% Save the new nelx value
guidata(hObject,handles)



% --- Executes on button press in gui_fx.
function gui_fx_Callback(hObject, eventdata, handles)
% hObject    handle to gui_fx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of gui_fx
fx = str2double(get(hObject, 'String'));
if isnan(fx)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save new value fx to handles
handles.lc.fx = fx;
guidata(hObject,handles)

% --- Executes on button press in gui_fy.
function gui_fy_Callback(hObject, eventdata, handles)
% hObject    handle to gui_fy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of gui_fy
fy = str2double(get(hObject, 'String'));
if isnan(fy)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save new value fy to handles
handles.lc.fy = fy;
guidata(hObject,handles)

% --- Executes on button press in gui_fz.
function gui_fz_Callback(hObject, eventdata, handles)
% hObject    handle to gui_fz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of gui_fz
fz = str2double(get(hObject, 'String'));
if isnan(fz)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save new value fx to handles
handles.lc.fz = fz;
guidata(hObject,handles)

function lc_x1_Callback(hObject, eventdata, handles)
% hObject    handle to lc_x1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lc_x1 as text
%        str2double(get(hObject,'String')) returns contents of lc_x1 as a double
x1 = str2double(get(hObject, 'String'));
if isnan(x1)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
    handles.lc.x1 = 0;
elseif mod(x1,1);
    set(hObject, 'String', floor(x1));
    errordlg('Input must be a real integer','Error');
    handles.lc.x1 = floor(x1);
elseif x1 > handles.domain.nelx;
    set(hObject, 'String', handles.domain.nelx);
    errordlg('Input must less than nelx','Error');
    handles.lc.x1 = handles.domain.nelx;
else
    handles.lc.x1 = x1;
end
% Save the new nelx value
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function lc_x1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lc_x1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lc_y1_Callback(hObject, eventdata, handles)
% hObject    handle to lc_y1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lc_y1 as text
%        str2double(get(hObject,'String')) returns contents of lc_y1 as a double
y1 = str2double(get(hObject, 'String'));
if isnan(y1)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
    handles.lc.y1 = handles.domain.nely;
elseif mod(y1,1);
    set(hObject, 'String', floor(y1));
    errordlg('Input must be a real integer','Error');
    handles.lc.y1 = floor(y1);
elseif y1 > handles.domain.nely;
    set(hObject, 'String', handles.domain.nely);
    errordlg('Input must less than nely','Error');
    handles.lc.y1 = handles.domain.nely;
else
    handles.lc.y1 = y1;
end

% Save the new nelx value
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function lc_y1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lc_y1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lc_z1_Callback(hObject, eventdata, handles)
% hObject    handle to lc_z1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lc_z1 as text
%        str2double(get(hObject,'String')) returns contents of lc_z1 as a double
z1 = str2double(get(hObject, 'String'));
if isnan(z1)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
    handles.lc.z1 = 0;
elseif mod(z1,1);
    set(hObject, 'String', floor(z1));
    errordlg('Input must be a real integer','Error');
    handles.lc.z1 = floor(z1);
elseif z1 > handles.domain.nelz;
    set(hObject, 'String', handles.domain.nelz);
    errordlg('Input must less than nelz','Error');
    handles.lc.z1 = handles.domain.nelz;
else
    handles.lc.z1 = z1;
end
% Save the new nelx value

guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function lc_z1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lc_z1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lc_x2_Callback(hObject, eventdata, handles)
% hObject    handle to lc_x2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lc_x2 as text
%        str2double(get(hObject,'String')) returns contents of lc_x2 as a double
x2 = str2double(get(hObject, 'String'));
if isnan(x2)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
    handles.lc.x2 = 0;
elseif mod(x2,1);
    set(hObject, 'String', floor(x2));
    errordlg('Input must be a real integer','Error');
    handles.lc.x2 = floor(x2);
elseif x2 > handles.domain.nelx;
    set(hObject, 'String', handles.domain.nelx);
    errordlg('Input must less than nelx','Error');
    handles.lc.x2 = handles.domain.nelx;
else
    handles.lc.x2 = x2;
end
% Save the new nelx value
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function lc_x2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lc_x2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lc_y2_Callback(hObject, eventdata, handles)
% hObject    handle to lc_y2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lc_y2 as text
%        str2double(get(hObject,'String')) returns contents of lc_y2 as a double
y2 = str2double(get(hObject, 'String'));
if isnan(y2)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
    handles.lc.y2 = 0;
elseif mod(y2,1);
    set(hObject, 'String', floor(y2));
    errordlg('Input must be a real integer','Error');
    handles.lc.y2 = floor(y2);
elseif y2 > handles.domain.nely;
    set(hObject, 'String', handles.domain.nely);
    errordlg('Input must less than nely','Error');
    handles.lc.y2 = handles.domain.nely;
else
    handles.lc.y2 = y2;
end
% Save the new nelx value
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function lc_y2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lc_y2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lc_z2_Callback(hObject, eventdata, handles)
% hObject    handle to lc_z2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lc_z2 as text
%        str2double(get(hObject,'String')) returns contents of lc_z2 as a double
z2 = str2double(get(hObject, 'String'));
if isnan(z2)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
elseif mod(z2,1);
    set(hObject, 'String', floor(z2));
    errordlg('Input must be a real integer','Error');
    handles.lc.z2 = floor(z2);
elseif z2 > handles.domain.nelz;
    set(hObject, 'String', handles.domain.nelz);
    errordlg('Input must less than nelz','Error');
    handles.lc.z2 = handles.domain.nelz;
else
    handles.lc.z2 = z2;
end

% Save the new nelx value
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function lc_z2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lc_z2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_run.
function btn_run_Callback(hObject, eventdata, handles)
% hObject    handle to btn_run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
eflag = 0;

if isempty(handles.bc.fixeddof)
    errordlg('Boundary Condition not defined','Error');
    eflag = eflag + 1;
end

if isempty(handles.lc.loaddof) && isempty(handles.lc.loadscale)
    errordlg('Loading Condition not defined','Error');
    eflag = eflag + 1;
end

if eflag
    return
end

% Freeze Colormap
axes(handles.axes_design)
colormap(flipud(gray));
freezeColors
    
axes(handles.axes_res);

if isfield(handles,'patch') && ~isempty(handles.patch)
    if isfield(handles.patch,'result') && ~isempty(handles.patch.result)
        delete(handles.patch.result);
    end
end

if handles.optpar.displayflag
    set(handles.axes_res,'Visible','on');
else
    set(handles.axes_res,'Visible','off');
end

[handles.patch.result, handles.xopt, handles.disp] = top3d(handles);

set_PlotCtrl(hObject, eventdata, handles);

rgb(1,1) = get(handles.edit_Red,'val');
rgb(1,2) = get(handles.edit_Green,'val');
rgb(1,3) = get(handles.edit_Blue,'val');
handles.optpar.fcolor = rgb./255;
handles.optpar.cutoff = get(handles.plot_cutoff,'Value');
handles.optpar.alpha  = get(handles.plot_alpha,'Value');
handles.optpar.plot_type = 'patch';
guidata(hObject, handles);



% --- Executes when selected object is changed in load_type.
function load_type_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in load_type
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'btn_lc1'
        % Code for when btn_1 is selected.
        handles.lc.type = 'individual';
    case 'btn_lc2'
        % Code for when btn_2 is selected.
        handles.lc.type = 'distribute';
end
guidata(hObject,handles)


% --- Executes on button press in btn_matlcon.
function btn_matlcon_Callback(hObject, eventdata, handles)
% hObject    handle to btn_matlcon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function gui_E0_Callback(hObject, eventdata, handles)
% hObject    handle to gui_E0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gui_E0 as text
%        str2double(get(hObject,'String')) returns contents of gui_E0 as a double
E0 = str2double(get(hObject, 'String'));
if isnan(E0)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new nelx value
handles.matlcon.E0 = E0;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function gui_E0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gui_E0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gui_nu_Callback(hObject, eventdata, handles)
% hObject    handle to gui_nu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gui_nu as text
%        str2double(get(hObject,'String')) returns contents of gui_nu as a double
nu = str2double(get(hObject, 'String'));
if isnan(nu)
    set(hObject, 'String', 0.3);
    errordlg('Input must be a number','Error');
end

% Save the new nelx value
handles.matlcon.nu = nu;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function gui_nu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gui_nu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gui_p_Callback(hObject, eventdata, handles)
% hObject    handle to gui_p (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gui_p as text
%        str2double(get(hObject,'String')) returns contents of gui_p as a double
p = str2double(get(hObject, 'String'));
if isnan(p)
    set(hObject, 'String', 0.3);
    errordlg('Input must be a number','Error');
end

% Save the new nelx value
handles.matlcon.p = p;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function gui_p_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gui_p (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit32_Callback(hObject, eventdata, handles)
% hObject    handle to gui_volfrac (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gui_volfrac as text
%        str2double(get(hObject,'String')) returns contents of gui_volfrac as a double


% --- Executes during object creation, after setting all properties.
function edit32_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gui_volfrac (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gui_rmin_Callback(hObject, eventdata, handles)
% hObject    handle to gui_rmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gui_rmin as text
%        str2double(get(hObject,'String')) returns contents of gui_rmin as a double
rmin = str2double(get(hObject, 'String'));
if isnan(rmin)
    set(hObject, 'String', 1.2);
    errordlg('Input must be a number','Error');
end

handles.optpar.rmin = rmin;
% Save the new nelx value
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function gui_rmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gui_rmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gui_maxloop_Callback(hObject, eventdata, handles)
% hObject    handle to gui_maxloop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gui_maxloop as text
%        str2double(get(hObject,'String')) returns contents of gui_maxloop as a double
maxloop = str2double(get(hObject, 'String'));
if isnan(maxloop)
    set(hObject, 'String', 200);
    errordlg('Input must be a number','Error');
end

if mod(maxloop,1);
    set(hObject, 'String', 200);
    errordlg('Input must be a real integer','Error');
end

handles.optpar.maxloop = maxloop;
% Save the new nelx value
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function gui_maxloop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gui_maxloop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gui_tolx_Callback(hObject, eventdata, handles)
% hObject    handle to gui_tolx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gui_tolx as text
%        str2double(get(hObject,'String')) returns contents of gui_tolx as a double
tolx = str2double(get(hObject, 'String'));
if isnan(tolx)
    set(hObject, 'String', 0.01);
    errordlg('Input must be a number','Error');
end

handles.optpar.tolx = tolx;
% Save the new nelx value
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function gui_tolx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gui_tolx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in chk_plotflag.
function chk_plotflag_Callback(hObject, eventdata, handles)
% hObject    handle to chk_plotflag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chk_plotflag
if (get(hObject,'Value') == get(hObject,'Max'))
    % Checkbox is checked-take appropriate action
    % Save the new nelx value
    handles.optpar.displayflag = 1;
else
    % Checkbox is not checked-take appropriate action
    handles.optpar.displayflag = 0;
end
guidata(hObject,handles)


% --- Executes on button press in chk_randinit.
function chk_randinit_Callback(hObject, eventdata, handles)
% hObject    handle to chk_randinit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chk_randinit
if (get(hObject,'Value') == get(hObject,'Max'))
    % Checkbox is checked-take appropriate action
    % Save the new nelx value
    handles.optpar.randinit = 1;
else
    % Checkbox is not checked-take appropriate action
    handles.optpar.randinit = 0;
end
guidata(hObject,handles)



function gui_iter_Callback(hObject, eventdata, handles)
% hObject    handle to gui_iter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gui_iter as text
%        str2double(get(hObject,'String')) returns contents of gui_iter as a double


% --- Executes during object creation, after setting all properties.
function gui_iter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gui_iter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function set_DesignPanel(handles, status)
if nargin == 1
    status = 'default';
end

switch status
    case 'default'
        set(handles.gui_nelx,'String',30);
        set(handles.gui_nely,'String',10);
        set(handles.gui_nelz,'String',2);
        set(handles.gui_volfrac,'String',0.3);
        set(handles.gui_E0,'String',1);
        set(handles.gui_nu,'String',0.3);
        set(handles.gui_p,'String',3);
        set(handles.gui_maxloop,'String',200);
        set(handles.gui_tolx,'String',0.01);
        set(handles.chk_plotflag,'Value',1);
        set(handles.gui_rmin,'String',1.2);
        set(handles.chk_randinit,'Value',0);
    case 'open'
        set(handles.gui_nelx,'String', handles.domain.nelx);
        set(handles.gui_nely,'String', handles.domain.nely);
        set(handles.gui_nelz,'String', handles.domain.nelz);
        set(handles.gui_volfrac,'String', handles.domain.volfrac);
        set(handles.gui_E0,'String', handles.matlcon.E0);
        set(handles.gui_nu,'String', handles.matlcon.nu);
        set(handles.gui_p,'String', handles.matlcon.p);
        set(handles.gui_maxloop,'String', handles.optpar.maxloop);
        set(handles.gui_tolx,'String', handles.optpar.tolx);
        set(handles.chk_plotflag,'Value', handles.optpar.displayflag);
        set(handles.gui_rmin,'String', handles.optpar.rmin);
        set(handles.chk_randinit,'Value', handles.optpar.randinit); 
end


function get_DesignPanel(hObject, eventdata, handles, status)
% - Domain
handles.domain.nelx = str2double(get(handles.gui_nelx, 'String'));
handles.domain.nely = str2double(get(handles.gui_nely, 'String'));
handles.domain.nelz = str2double(get(handles.gui_nelz, 'String'));
handles.domain.volfrac  = str2double(get(handles.gui_volfrac, 'String'));
% - Matlcon
handles.matlcon.E0 = str2double(get(handles.gui_E0, 'String'));
handles.matlcon.nu = str2double(get(handles.gui_nu, 'String'));
handles.matlcon.p  = str2double(get(handles.gui_p, 'String'));
% - OptPar
handles.optpar.maxloop = str2double(get(handles.gui_maxloop, 'String'));
handles.optpar.tolx    = str2double(get(handles.gui_tolx, 'String'));
handles.optpar.displayflag = get(handles.chk_plotflag,'Value');
handles.optpar.rmin    = str2double(get(handles.gui_rmin, 'String'));
handles.optpar.randinit = get(handles.chk_randinit,'Value');

if strcmp(status,'reset')
    axes(handles.axes_design);
    if isfield(handles,'patch') && ~isempty(handles.patch)
        if isfield(handles.patch,'design') && ~isempty(handles.patch.design)
            delete(handles.patch.design)
        end
    end
    
    if ~isempty(handles.bc.fixeddof)
        btn_bc_reset_Callback(hObject, eventdata, handles);
        handles = guidata(hObject);
    end
    
    if ~isempty(handles.lc.loaddof) || ~isempty(handles.lc.loadscale)
        btn_lc_reset_Callback(hObject, eventdata, handles);
        handles = guidata(hObject);
    end
    
    x = repmat(handles.domain.volfrac,[handles.domain.nely,handles.domain.nelx,handles.domain.nelz]);
    
    p_design = plot_3d(x, 1e-3, 0.5);
    
    handles.patch.design = p_design;
end

guidata(hObject,handles);


% --- Executes on button press in btn_outer_reset.
function btn_outer_reset_Callback(hObject, eventdata, handles)
% hObject    handle to btn_outer_reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set_DesignPanel(handles);
get_DesignPanel(hObject, eventdata, handles,'reset');

% --- Executes on button press in btn_bc_reset.
function btn_bc_reset_Callback(hObject, eventdata, handles)
% hObject    handle to btn_bc_reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set_BCPanel(handles, 'reset');
get_BCPanel(hObject, eventdata, handles,'initial')

function set_BCPanel(handles,status)
switch status
    case 'default'
        set(handles.chk_ux, 'Value', 1);
        set(handles.chk_uy, 'Value', 1);
        set(handles.chk_uz, 'Value', 1);
        set(handles.btn_bc1, 'Value', 0);
        set(handles.btn_bc2, 'Value', 1);
        set(handles.bc_x1, 'String', 0);
        set(handles.bc_y1, 'String', 0);
        set(handles.bc_z1, 'String', 0);
        set(handles.bc_x2, 'String', 0);
        set(handles.bc_y2, 'String', 10);
        set(handles.bc_z2, 'String', 2);
    case 'reset'
        set(handles.chk_ux,'Value',0);
        set(handles.chk_uy,'Value',0);
        set(handles.chk_uz,'Value',0);
        set(handles.btn_bc1,'Value',1);
        set(handles.btn_bc2,'Value',0);
        set(handles.bc_x1,'String','');
        set(handles.bc_y1,'String','');
        set(handles.bc_z1,'String','');
        set(handles.bc_x2,'String','');
        set(handles.bc_y2,'String','');
        set(handles.bc_z2,'String','');
        %         axes(handles.axes_design);
        if isfield(handles.bc,'htri') && ~isempty(handles.bc.htri)
            delete(handles.bc.htri);
        end
    case 'open'
        set(handles.chk_ux,'Value', handles.bc.ux);
        set(handles.chk_uy,'Value', handles.bc.uy);
        set(handles.chk_uz,'Value', handles.bc.uz);        
        switch handles.bc.type
            case 'individual'
                set(handles.btn_bc1,'Value',1);
                set(handles.btn_bc2,'Value',0);
            case 'distribute'
                set(handles.btn_bc1,'Value',0);
                set(handles.btn_bc2,'Value',1);
        end
        % x1
        if isfield(handles.bc, 'x1')
            set(handles.bc_x1,'String', handles.bc.x1);
        else
            set(handles.bc_x1,'String', '');
        end
        % y1
        if isfield(handles.bc, 'y1')
            set(handles.bc_y1,'String', handles.bc.y1);
        else
            set(handles.bc_y1,'String', '');
        end
        % z1
        if isfield(handles.bc, 'z1')
            set(handles.bc_z1,'String', handles.bc.z1);
        else
            set(handles.bc_z1,'String', '');
        end
        % x2
        if isfield(handles.bc, 'x2')
            set(handles.bc_x2,'String', handles.bc.x2);
        else
            set(handles.bc_x2,'String', '');
        end
        % y2
        if isfield(handles.bc, 'y2')
            set(handles.bc_y2,'String', handles.bc.y2);
        else
            set(handles.bc_y2,'String', '');
        end
        % z2
        if isfield(handles.bc, 'z2')
            set(handles.bc_z2,'String', handles.bc.z2);
        else
            set(handles.bc_z2,'String', '');
        end
end

function get_BCPanel(hObject, eventdata, handles, status)
handles.bc.ux       = get(handles.chk_ux,'Value');
handles.bc.uy       = get(handles.chk_uy,'Value');
handles.bc.uz       = get(handles.chk_uz,'Value');
if get(handles.btn_bc1,'Value') == get(handles.btn_bc1,'Max')
    handles.bc.type     = 'individual';
elseif get(handles.btn_bc2,'Value') == get(handles.btn_bc2,'Max')
    handles.bc.type     = 'distribute';
end

handles.bc.fixeddof = [];
handles.bc.htri     = [];


if strcmp(status, 'default')
    handles.bc.x1 = str2double(get(handles.bc_x1, 'String'));
    handles.bc.y1 = str2double(get(handles.bc_y1, 'String'));
    handles.bc.z1 = str2double(get(handles.bc_z1, 'String'));
    handles.bc.x2 = str2double(get(handles.bc_x2, 'String'));
    handles.bc.y2 = str2double(get(handles.bc_y2, 'String'));
    handles.bc.z2 = str2double(get(handles.bc_z2, 'String'));
end

% Save new data
guidata(hObject,handles)

% --- Executes on button press in btn_default.
function btn_default_Callback(hObject, eventdata, handles)
% hObject    handle to btn_default (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set_DesignPanel(handles);
get_DesignPanel(hObject, eventdata, handles,'initial');
handles = guidata(hObject);
set_BCPanel(handles, 'default');
get_BCPanel(hObject, eventdata, handles, 'default');
handles = guidata(hObject);
set_LCPanel(handles, 'default');
get_LCPanel(hObject, eventdata, handles, 'default');
handles = guidata(hObject);


btn_outer_Callback(hObject, eventdata, handles);
handles = guidata(hObject);
btn_bc_Callback(hObject, eventdata, handles);
handles = guidata(hObject);
btn_lc_Callback(hObject, eventdata, handles);
handles = guidata(hObject);

% Save new data
guidata(hObject, handles);


% --- Executes on button press in btn_lc_reset.
function btn_lc_reset_Callback(hObject, eventdata, handles)
% hObject    handle to btn_lc_reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set_LCPanel(handles, 'reset');
get_LCPanel(hObject, eventdata, handles,'initial')

function set_LCPanel(handles,status)
switch status
    case 'default'
        set(handles.gui_fx,'String',0);
        set(handles.gui_fy,'String',-1);
        set(handles.gui_fz,'String',0);
        set(handles.btn_lc1,'Value',1);
        set(handles.btn_lc2,'Value',0);
        set(handles.lc_x1,'String',30);
        set(handles.lc_y1,'String',5);
        set(handles.lc_z1,'String',1);
        set(handles.lc_x2,'String','');
        set(handles.lc_y2,'String','');
        set(handles.lc_z2,'String','');
    case 'reset'
        set(handles.gui_fx,'String',0);
        set(handles.gui_fy,'String',0);
        set(handles.gui_fz,'String',0);
        set(handles.btn_lc1,'Value',1);
        set(handles.btn_lc2,'Value',0);
        set(handles.lc_x1,'String','');
        set(handles.lc_y1,'String','');
        set(handles.lc_z1,'String','');
        set(handles.lc_x2,'String','');
        set(handles.lc_y2,'String','');
        set(handles.lc_z2,'String','');
        
        if isfield(handles.lc,'harr') && ~isempty(handles.lc.harr)
            delete(handles.lc.harr);
        end
    case 'open'
        set(handles.gui_fx,'String', handles.lc.fx);
        set(handles.gui_fy,'String', handles.lc.fy);
        set(handles.gui_fz,'String', handles.lc.fz);       
        switch handles.lc.type
            case 'individual'
                set(handles.btn_lc1,'Value',1);
                set(handles.btn_lc2,'Value',0);
            case 'distribute'
                set(handles.btn_lc1,'Value',0);
                set(handles.btn_lc2,'Value',1);
        end
        % x1
        if isfield(handles.lc, 'x1')
            set(handles.lc_x1,'String', handles.bc.x1);
        else
            set(handles.lc_x1,'String', '');
        end
        % y1
        if isfield(handles.lc, 'y1')
            set(handles.lc_y1,'String', handles.lc.y1);
        else
            set(handles.lc_y1,'String', '');
        end
        % z1
        if isfield(handles.lc, 'z1')
            set(handles.lc_z1,'String', handles.lc.z1);
        else
            set(handles.lc_z1,'String', '');
        end
        % x2
        if isfield(handles.lc, 'x2')
            set(handles.lc_x2,'String', handles.lc.x2);
        else
            set(handles.lc_x2,'String', '');
        end
        % y2
        if isfield(handles.lc, 'y2')
            set(handles.lc_y2,'String', handles.lc.y2);
        else
            set(handles.lc_y2,'String', '');
        end
        % z2
        if isfield(handles.lc, 'z2')
            set(handles.lc_z2,'String', handles.lc.z2);
        else
            set(handles.lc_z2,'String', '');
        end
end

function get_LCPanel(hObject, eventdata, handles, status)

handles.lc.fx       = str2double(get(handles.gui_fx, 'String'));
handles.lc.fy       = str2double(get(handles.gui_fy, 'String'));
handles.lc.fz       = str2double(get(handles.gui_fz, 'String'));
if get(handles.btn_lc1,'Value') == get(handles.btn_bc1,'Max')
    handles.lc.type     = 'individual';
elseif get(handles.btn_lc2,'Value') == get(handles.btn_bc2,'Max')
    handles.lc.type     = 'distribute';
end
handles.lc.loaddof  = [];
handles.lc.loadscale  = [];
handles.lc.harr     = [];

if strcmp(status, 'default')
    handles.lc.x1 = str2double(get(handles.lc_x1, 'String'));
    handles.lc.y1 = str2double(get(handles.lc_y1, 'String'));
    handles.lc.z1 = str2double(get(handles.lc_z1, 'String'));
    handles.lc.x2 = str2double(get(handles.lc_x2, 'String'));
    handles.lc.y2 = str2double(get(handles.lc_y2, 'String'));
    handles.lc.z2 = str2double(get(handles.lc_z2, 'String'));
end

% Save new data
guidata(hObject, handles);


% --- Executes on button press in btn_reset.
function btn_reset_Callback(hObject, eventdata, handles)
% hObject    handle to btn_reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(top3dGUI);
top3dGUI;

% --- Executes on button press in btn_close.
function btn_close_Callback(hObject, eventdata, handles)
% hObject    handle to btn_close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg('Are you sure you want to close the Top3d GUI?',...
    'Close Request','Yes','No','Yes');

switch selection
    case 'Yes',
        delete(handles.figure1);
    case 'No'
        return
end


% --- Executes on button press in plot_box.
function plot_box_Callback(hObject, eventdata, handles)
% hObject    handle to plot_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plot_box
button_state = get(hObject,'Value');
if button_state == get(hObject,'Max')
    % Toggle button is pressed-take appropriate action
    axes(handles.axes_res); box on;
elseif button_state == get(hObject,'Min')
    % Toggle button is not pressed-take appropriate action
    axes(handles.axes_res); box off;
end


% --- Executes on button press in plot_axis.
function plot_axis_Callback(hObject, eventdata, handles)
% hObject    handle to plot_axis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plot_axis
button_state = get(hObject,'Value');
if button_state == get(hObject,'Max')
    % Toggle button is pressed-take appropriate action
    axes(handles.axes_res); axis on;
elseif button_state == get(hObject,'Min')
    % Toggle button is not pressed-take appropriate action
    axes(handles.axes_res); axis off;
end


% --- Executes on slider movement.
function plot_alpha_Callback(hObject, eventdata, handles)
% hObject    handle to plot_alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
alpha = get(hObject,'Value');
set(handles.plot_alpha_gui,'String', num2str(alpha,'%1.2f'));

axes(handles.axes_res);

handles.optpar.alpha = alpha;
handles.patch.result  = plot_3d(handles.xopt,handles.optpar.cutoff,...
    handles.optpar.alpha,handles.optpar.plot_type,handles.optpar.fcolor,handles.disp);

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function plot_alpha_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plot_alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function plot_alpha_gui_Callback(hObject, eventdata, handles)
% hObject    handle to plot_alpha_gui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of plot_alpha_gui as text
%        str2double(get(hObject,'String')) returns contents of plot_alpha_gui as a double
alpha = str2double(get(hObject, 'String'));
if isnan(alpha)
    set(hObject, 'String', 0.5);
    errordlg('Input must be a number','Error');
    return
end

set(handles.plot_alpha, 'Value', alpha)

axes(handles.axes_res);

handles.optpar.alpha = alpha;
handles.patch.result  = plot_3d(handles.xopt,handles.optpar.cutoff,...
    handles.optpar.alpha,handles.optpar.plot_type,handles.optpar.fcolor,handles.disp);


guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function plot_alpha_gui_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plot_alpha_gui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function plot_cutoff_Callback(hObject, eventdata, handles)
% hObject    handle to plot_cutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
cutoff = get(hObject, 'Value');
set(handles.plot_cutoff_gui, 'String', num2str(cutoff,'%1.2f'));

axes(handles.axes_res);

handles.optpar.cutoff = cutoff;
handles.patch.result  = plot_3d(handles.xopt,handles.optpar.cutoff,...
    handles.optpar.alpha,handles.optpar.plot_type,handles.optpar.fcolor,handles.disp);

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function plot_cutoff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plot_cutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function plot_cutoff_gui_Callback(hObject, eventdata, handles)
% hObject    handle to plot_cutoff_gui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of plot_cutoff_gui as text
%        str2double(get(hObject,'String')) returns contents of plot_cutoff_gui as a double
cutoff = str2double(get(hObject, 'String'));
if isnan(cutoff)
    set(hObject, 'String', 0.5);
    errordlg('Input must be a number','Error');
    return
end

set(handles.plot_cutoff, 'Value', cutoff)

axes(handles.axes_res);

handles.optpar.cutoff = cutoff;
handles.patch.result  = plot_3d(handles.xopt,handles.optpar.cutoff,...
    handles.optpar.alpha,handles.optpar.plot_type,handles.optpar.fcolor,handles.disp);

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function plot_cutoff_gui_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plot_cutoff_gui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function set_PlotCtrl(hObject, eventdata, handles)

set(handles.plot_box,'Value',1);
set(handles.plot_axis,'Value',1);
set(handles.plot_cutoff,'Value',0.5);
set(handles.plot_cutoff_gui,'String', 0.50);
set(handles.plot_alpha,'Value',1);
set(handles.plot_alpha_gui,'String', 1.00);

set(handles.textR,'Visible','off');
set(handles.textG,'Visible','off');
set(handles.textB,'Visible','off');
rgb = randi(256,[1 3])-1;
set(handles.edit_Red,'str',num2str(rgb(1)),'val',rgb(1),'Visible','off');
set(handles.edit_Green,'str',num2str(rgb(2)),'val',rgb(2),'Visible','off');
set(handles.edit_Blue,'str',num2str(rgb(3)),'val',rgb(3),'Visible','off');
set_colors(handles);
set(handles.plot_setcolor,'Visible','off');

set(handles.uipanel15,'Visible','on');


% --- Executes on selection change in plot_type.
function plot_type_Callback(hObject, eventdata, handles)
% hObject    handle to plot_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns plot_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from plot_type
val = get(hObject,'Value');
switch val
    case 1   % User selected the first item
        type = 'patch';
        set(handles.textR,'Visible','off');
        set(handles.textG,'Visible','off');
        set(handles.textB,'Visible','off');
        set(handles.edit_Red,'Visible','off');
        set(handles.edit_Green,'Visible','off');
        set(handles.edit_Blue,'Visible','off');
        set(handles.plot_setcolor,'Visible','off');
    case 2   % User selected the second item
        type = 'isosurface';
        set(handles.textR,'Visible','on');
        set(handles.textG,'Visible','on');
        set(handles.textB,'Visible','on');
        rgb = randi(256,[1 3])-1;
        set(handles.edit_Red,'str',num2str(rgb(1)),'val',rgb(1),'Visible','on');
        set(handles.edit_Green,'str',num2str(rgb(2)),'val',rgb(2),'Visible','on');
        set(handles.edit_Blue,'str',num2str(rgb(3)),'val',rgb(3),'Visible','on');
        set_colors(handles);
        set(handles.plot_setcolor,'Visible','on');
    case 3
        type = 'isonormals';
        set(handles.textR,'Visible','on');
        set(handles.textG,'Visible','on');
        set(handles.textB,'Visible','on');
        rgb = randi(256,[1 3])-1;
        set(handles.edit_Red,'str',num2str(rgb(1)),'val',rgb(1),'Visible','on');
        set(handles.edit_Green,'str',num2str(rgb(2)),'val',rgb(2),'Visible','on');
        set(handles.edit_Blue,'str',num2str(rgb(3)),'val',rgb(3),'Visible','on');
        set_colors(handles);
        set(handles.plot_setcolor,'Visible','on');
    case 4
        type = 'disp';
        set(handles.textR,'Visible','off');
        set(handles.textG,'Visible','off');
        set(handles.textB,'Visible','off');
        set(handles.edit_Red,'Visible','off');
        set(handles.edit_Green,'Visible','off');
        set(handles.edit_Blue,'Visible','off');
        set(handles.plot_setcolor,'Visible','off');
end

rgb(1,1) = get(handles.edit_Red,'val');
rgb(1,2) = get(handles.edit_Green,'val');
rgb(1,3) = get(handles.edit_Blue,'val');
handles.optpar.fcolor = rgb./255;

axes(handles.axes_res);
handles.optpar.plot_type  = type;

handles.patch.result  = plot_3d(handles.xopt,handles.optpar.cutoff,...
    handles.optpar.alpha,handles.optpar.plot_type,handles.optpar.fcolor, handles.disp);

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function plot_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plot_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in plot_setcolor.
function plot_setcolor_Callback(hObject, eventdata, handles)
% hObject    handle to plot_setcolor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

clr = uisetcolor('Set Face Color');

rgb = round(clr.*255);
set(handles.edit_Red,'str',num2str(rgb(1)),'val',rgb(1));
set(handles.edit_Green,'str',num2str(rgb(2)),'val',rgb(2));
set(handles.edit_Blue,'str',num2str(rgb(3)),'val',rgb(3));

set_colors(handles);

handles.optpar.fcolor = clr;

handles.patch.result  = plot_3d(handles.xopt,handles.optpar.cutoff,...
    handles.optpar.alpha,handles.optpar.plot_type,handles.optpar.fcolor);

guidata(hObject, handles);


function edit_Red_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Red (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Red as text
%        str2double(get(hObject,'String')) returns contents of edit_Red as a double
str = get(hObject,'string');  num = str2double(str);
num0 = get(hObject,'value');  str0 = num2str(num0);

if isnan(num)
    errordlg('Input must be a number','Error');
end

if num>255 || num<0
    num = nan; 
    errordlg('Input must between 0 and 255','Error');
end

if num>0 && num<1, num = num.*255; end
num = round(num);  str = num2str(num);
if ~isnan(num)
    set(hObject,'value',num,'string',str);
    clr = [num get(handles.edit_Green,'val') get(handles.edit_Blue,'val')]./255;
    set_colors(handles);
    
    handles.optpar.fcolor = clr;
    
    handles.patch.result  = plot_3d(handles.xopt,handles.optpar.cutoff,...
        handles.optpar.alpha,handles.optpar.plot_type,handles.optpar.fcolor);
    
else
    set(hObject,'string',str0);
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_Red_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Red (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Green_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Green (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Green as text
%        str2double(get(hObject,'String')) returns contents of edit_Green as a double
str = get(hObject,'string');  num = str2double(str);
num0 = get(hObject,'value');  str0 = num2str(num0);

if isnan(num)
    errordlg('Input must be a number','Error');
end

if num>255 || num<0
    num = nan; 
    errordlg('Input must between 0 and 255','Error');
end

if num>0 && num<1, num = num.*255; end
num = round(num);  str = num2str(num);
if ~isnan(num)
    set(hObject,'value',num,'string',str);
    clr = [get(handles.edit_Red,'val') num get(handles.edit_Blue,'val')]./255;
    set_colors(handles);
    
    handles.optpar.fcolor = clr;
    
    handles.patch.result  = plot_3d(handles.xopt,handles.optpar.cutoff,...
        handles.optpar.alpha,handles.optpar.plot_type,handles.optpar.fcolor);
    
else
    set(hObject,'string',str0);
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_Green_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Green (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Blue_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Blue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Blue as text
%        str2double(get(hObject,'String')) returns contents of edit_Blue as a double
str = get(hObject,'string');  num = str2double(str);
num0 = get(hObject,'value');  str0 = num2str(num0);

if isnan(num)
    errordlg('Input must be a number','Error');
end

if num>255 || num<0
    num = nan; 
    errordlg('Input must between 0 and 255','Error');
end

if num>0 && num<1, num = num.*255; end
num = round(num);  str = num2str(num);
if ~isnan(num)
    set(hObject,'value',num,'string',str);
    clr = [get(handles.edit_Red,'val') get(handles.edit_Green,'val') num]./255;
    set_colors(handles);
    
    handles.optpar.fcolor = clr;
    
    handles.patch.result  = plot_3d(handles.xopt,handles.optpar.cutoff,...
        handles.optpar.alpha,handles.optpar.plot_type,handles.optpar.fcolor);
    
else
    set(hObject,'string',str0);
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_Blue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Blue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function set_colors(handles)
rgb(1,1) = get(handles.edit_Red,'val');
rgb(1,2) = get(handles.edit_Green,'val');
rgb(1,3) = get(handles.edit_Blue,'val');
clr = rgb./255;
set([handles.textR handles.textG handles.textB],'BackGroundColor',clr);
set([handles.textR handles.textG handles.textB],'ForeGroundColor',1-clr);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over info.
% hObject    handle to info (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Help_Callback(hObject, eventdata, handles)
% hObject    handle to Help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function About_Callback(hObject, eventdata, handles)
% hObject    handle to About (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
msgbox({'Top3d GUI beta' 'http://top3dapp.com'});

% --------------------------------------------------------------------
function Demo_menu_Callback(hObject, eventdata, handles)
% hObject    handle to Demo_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
btn_default_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function File_Callback(hObject, eventdata, handles)
% hObject    handle to File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function New_Callback(hObject, eventdata, handles)
% hObject    handle to New (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
btn_reset_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function Export_Callback(hObject, eventdata, handles)
% hObject    handle to Export (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function Quit_Callback(hObject, eventdata, handles)
% hObject    handle to Quit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
btn_close_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function Load_Callback(hObject, eventdata, handles)
% hObject    handle to Load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName, PathName] = uigetfile({'*.top3d', 'Top3d Design (*.top3d)'; '*.*', 'All Files (*.*)'});
if FileName~= 0
    load([PathName, FileName], '-mat')
    fprintf(1, 'Successfully loaded Top3d design file: %s\n', [PathName, FileName]);
else
    return
end

if ~exist('domain', 'var')
    errordlg('Invalid Top3d Design file','Import Error');
    return
end
% + Retrive data
handles.domain  = domain;
handles.matlcon = matlcon;
handles.optpar  = optpar;
handles.bc      = bc;
handles.lc      = lc;
% + Update GUI
set(handles.axes_res,'Visible','off');
handles.patch.result = [];
set(handles.uipanel15,'Visible','off');
axes(handles.axes_design), cla;
% + Design domain panel
set_DesignPanel(handles, 'open');
% --- Draw design domain
if handles.optpar.randinit
    x = 0.001+(1-0.001)*rand([handles.domain.nely,handles.domain.nelx,handles.domain.nelz]);
else
    x = repmat(handles.domain.volfrac,[handles.domain.nely,handles.domain.nelx,handles.domain.nelz]);
end
p_design = plot_3d(x, 1e-3, 0.5);
handles.patch.design = p_design;
% --- Update handles
guidata(hObject, handles)
% + Boundary Condtions panel
set_BCPanel(handles, 'open');
% --- Draw triangles
if isfield(handles.bc, 'x1')
    btn_bc_Callback(hObject, eventdata, handles);
    handles = guidata(hObject);
end
% + Load condition panel
set_LCPanel(handles, 'open');
if isfield(handles.lc, 'x1')
    btn_lc_Callback(hObject, eventdata, handles);
    handles = guidata(hObject);
end

% + Save again
guidata(hObject, handles);

% --------------------------------------------------------------------
function Save_Callback(hObject, eventdata, handles)
% hObject    handle to Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
domain  = handles.domain; %#ok<*NASGU>
matlcon = handles.matlcon;
optpar  = handles.optpar;
% Supprots Panel
bc      = handles.bc;
bc.fixeddof = [];
bc.loaddof  = [];
bc.htri     = [];
if isfield(bc, 'x1') && isnan(bc.x1), bc = rmfield(bc, 'x1'); end
if isfield(bc, 'y1') && isnan(bc.y1), bc = rmfield(bc, 'y1'); end
if isfield(bc, 'z1') && isnan(bc.z1), bc = rmfield(bc, 'z1'); end
if isfield(bc, 'x2') && isnan(bc.x2), bc = rmfield(bc, 'x2'); end
if isfield(bc, 'y2') && isnan(bc.y2), bc = rmfield(bc, 'y2'); end
if isfield(bc, 'z2') && isnan(bc.z2), bc = rmfield(bc, 'z2'); end
% Loads Panel
lc      = handles.lc;
lc.fixeddof  = [];
lc.loaddof   = [];
lc.loadscale = [];
lc.harr      = [];
if isfield(lc, 'x1') && isnan(lc.x1), lc = rmfield(lc, 'x1'); end
if isfield(lc, 'y1') && isnan(lc.y1), lc = rmfield(lc, 'y1'); end
if isfield(lc, 'z1') && isnan(lc.z1), lc = rmfield(lc, 'z1'); end
if isfield(lc, 'x2') && isnan(lc.x2), lc = rmfield(lc, 'x2'); end
if isfield(lc, 'y2') && isnan(lc.y2), lc = rmfield(lc, 'y2'); end
if isfield(lc, 'z2') && isnan(lc.z2), lc = rmfield(lc, 'z2'); end
% Save files
[FileName, PathName] = uiputfile({'*.top3d', 'Top3d Design (*.top3d)'}, 'Save Top3d Design');
if FileName ~= 0
    save([PathName, FileName], '-mat', 'domain', 'matlcon', 'optpar', 'bc', 'lc');
    fprintf(1, 'Successfully wrote Top3d design file: %s\n', [PathName, FileName]);
end


% --------------------------------------------------------------------
function Toworksapce_Callback(hObject, eventdata, handles)
% hObject    handle to Toworksapce (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    assignin('base', 'xPhys', handles.xopt)
catch err
    if ~isfield(handles, 'xopt')
      errordlg('Please run an optimization first','Export Error');
    else
      rethrow(err)
    end
end

fprintf(1, 'Successfully wrote result to workspace\n');

% --------------------------------------------------------------------
function Tofile_Callback(hObject, eventdata, handles)
% hObject    handle to Tofile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles, 'xopt')
    xPhys = handles.xopt;
    uisave('xPhys', 'xPhys')
else
    errordlg('Please run an optimization first','Export Error');
end

fprintf(1, 'Successfully wrote file\n');


% --------------------------------------------------------------------
function Print_Callback(hObject, eventdata, handles)
% hObject    handle to Print (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function PrintDesign_Callback(hObject, eventdata, handles)
% hObject    handle to PrintDesign (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% + Get file name and path
[FileName, PathName] = uiputfile({'*.bmp;*.eps;*.jpg;*.pdf;*.png;*.tif','All Image Files (*.bmp;*.eps;*.jpg;*.pdf;*.png;*.tif)'}, ...
    'Save Image');

if FileName == 0
    return
end
[~, ~, ext] = fileparts(FileName);

f_tmp = figure('visible', 'off');
copyobj(handles.axes_design, f_tmp);
colormap(flipud(gray));
print(f_tmp, sprintf('-d%s',ext(2:end)), [PathName FileName]);
close(f_tmp);

fprintf(1, 'Successfully wrote file: %s\n', [PathName, FileName]);

% --------------------------------------------------------------------
function PrintResult_Callback(hObject, eventdata, handles)
% hObject    handle to PrintResult (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles, 'xopt')
    errordlg('Please run an optimization first','Export Error');
    return
end

[FileName, PathName] = uiputfile({'*.bmp;*.eps;*.jpg;*.pdf;*.png;*.tiff','All Image Files (*.bmp;*.eps;*.jpg;*.pdf;*.png;*.tif)'}, ...
    'Save Image');

if FileName == 0
    return
end
[~, ~, ext] = fileparts(FileName);

f_tmp = figure('visible', 'off');
copyobj(handles.axes_res, f_tmp);

val = get(handles.plot_type,'Value');
if val == 4
    colormap(jet);
else
    colormap(flipud(gray));
end
print(f_tmp, sprintf('-d%s',ext(2:end)), [PathName FileName]);
close(f_tmp);

fprintf(1, 'Successfully wrote file: %s\n', [PathName, FileName]);
