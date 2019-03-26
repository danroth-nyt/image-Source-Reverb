function varargout = isGUI(varargin)
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @isGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @isGUI_OutputFcn, ...
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


% --- Executes just before isGUI is made visible.
function isGUI_OpeningFcn(hObject, eventdata, handles, varargin)
handles.Fs = 44100;

%absorption coefficients
handles.marble = 0.005;
handles.wood = 0.1;
handles.heavyFabric = 0.42;
handles.lightFabric = 0.14;
handles.glass = 0.03;
handles.plaster = 0.05;
handles.tiles = 0.03;
handles.carpet = 0.37;

%initial reflection coefficients
handles.aL = 1 - handles.plaster;
handles.aR = 1 - handles.plaster;
handles.aF = 1 - handles.heavyFabric;
handles.aB = 1 - handles.plaster;
handles.aC = 1 - handles.tiles;
handles.aG = 1 - handles.carpet;

plotter_Callback(hObject, eventdata, handles); %initialize plot

% Update handles structure
handles.output = hObject;

guidata(hObject, handles);

% initialize IR
calcIR_Callback(hObject, eventdata, handles)

% UIWAIT makes isGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = isGUI_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;


%% Edit Text Boxes

function Lx_Callback(hObject, eventdata, handles)
input = str2num(get(handles.Lx, 'String'));
if input < 1
    input = 1;
end
set(handles.Lx, 'Value', input);
set(handles.Lx, 'String', num2str(input));
guidata(hObject, handles);
plotter_Callback(hObject, eventdata, handles)
calcIR_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function Lx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Lx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Ly_Callback(hObject, eventdata, handles)
input = str2num(get(handles.Ly, 'String'));
if input < 1
    input = 1;
end
set(handles.Ly, 'Value', input);
set(handles.Ly, 'String', num2str(input));
guidata(hObject, handles);
plotter_Callback(hObject, eventdata, handles)
calcIR_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function Ly_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ly (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Lz_Callback(hObject, eventdata, handles)
input = str2num(get(handles.Lz, 'String'));
if input < 1
    input = 1;
end
set(handles.Lz, 'Value', input);
set(handles.Lz, 'String', num2str(input));
guidata(hObject, handles);
plotter_Callback(hObject, eventdata, handles)
calcIR_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function Lz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Lz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function p_Callback(hObject, eventdata, handles)
input = str2num(get(handles.p, 'String'));
if input < 1
    input = 1;
end
set(handles.p, 'Value', input);
set(handles.p, 'String', num2str(input));
guidata(hObject, handles);
plotter_Callback(hObject, eventdata, handles)
calcIR_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function p_CreateFcn(hObject, eventdata, handles)
% hObject    handle to p (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function q_Callback(hObject, eventdata, handles)
input = str2num(get(handles.q, 'String'));
if input < 1
    input = 1;
end
set(handles.q, 'Value', input);
set(handles.q, 'String', num2str(input));
guidata(hObject, handles);
plotter_Callback(hObject, eventdata, handles)
calcIR_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function q_CreateFcn(hObject, eventdata, handles)
% hObject    handle to q (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function r_Callback(hObject, eventdata, handles)
input = str2num(get(handles.r, 'String'));
if input < 1
    input = 1;
end
set(handles.r, 'Value', input);
set(handles.r, 'String', num2str(input));
guidata(hObject, handles);
plotter_Callback(hObject, eventdata, handles)
calcIR_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function r_CreateFcn(hObject, eventdata, handles)
% hObject    handle to r (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function a_Callback(hObject, eventdata, handles)
input = str2num(get(handles.a, 'String'));
if input < 1
    input = 1;
end
set(handles.a, 'Value', input);
set(handles.a, 'String', num2str(input));
guidata(hObject, handles);
plotter_Callback(hObject, eventdata, handles)
calcIR_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function a_CreateFcn(hObject, eventdata, handles)
% hObject    handle to a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function b_Callback(hObject, eventdata, handles)
input = str2num(get(handles.b, 'String'));
if input < 1
    input = 1;
end
set(handles.b, 'Value', input);
set(handles.b, 'String', num2str(input));
guidata(hObject, handles);
plotter_Callback(hObject, eventdata, handles)
calcIR_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function b_CreateFcn(hObject, eventdata, handles)
% hObject    handle to b (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function c_Callback(hObject, eventdata, handles)
input = str2num(get(handles.c, 'String'));
if input < 1
    input = 1;
end
set(handles.c, 'Value', input);
set(handles.c, 'String', num2str(input));
guidata(hObject, handles);
plotter_Callback(hObject, eventdata, handles)
calcIR_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function c_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function N_Callback(hObject, eventdata, handles)
input = str2num(get(handles.N, 'String'));
if input < 1
    input = 1;
end
set(handles.N, 'Value', input);
set(handles.N, 'String', num2str(input));
guidata(hObject, handles);
plotter_Callback(hObject, eventdata, handles)
calcIR_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function N_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function butterOrder_Callback(hObject, eventdata, handles)
input = str2num(get(handles.butterOrder, 'String'));
if input < 1
    input = 1;
end
set(handles.butterOrder, 'Value', input);
set(handles.butterOrder, 'String', num2str(input));
guidata(hObject, handles);
plotter_Callback(hObject, eventdata, handles)
calcIR_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function butterOrder_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function butterCutoff_Callback(hObject, eventdata, handles)
input = str2num(get(handles.butterCutoff, 'String'));
if input < 1
    input = 1;
end
set(handles.butterCutoff, 'Value', input);
set(handles.butterCutoff, 'String', num2str(input));
guidata(hObject, handles);
plotter_Callback(hObject, eventdata, handles)
calcIR_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function butterCutoff_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Material Selection Menus

% --- Executes on selection change in wL.
function wL_Callback(hObject, eventdata, handles)
choice = get(handles.wL, 'Value');
switch choice
    case 1 %marble
        temp = 1 - handles.marble;
    case 2 %wood
        temp = 1 - handles.wood;
    case 3 %heavy fabric
        temp = 1 - handles.heavyFabric;
    case 4 %light fabric
        temp = 1 - handles.lightFabric;
    case 5 %glass
        temp = 1 - handles.glass;
    case 6 %plaster
        temp = 1 - handles.plaster;
    case 7 %tiles
        temp = 1 - handles.tiles;
    case 8 %carpet
        temp = 1 - handles.carpet;
end
handles.aL = temp;
guidata(hObject, handles);
calcIR_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function wL_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in wR.
function wR_Callback(hObject, eventdata, handles)
choice = get(handles.wR, 'Value');
switch choice
    case 1 %marble
        temp = 1 - handles.marble;
    case 2 %wood
        temp = 1 - handles.wood;
    case 3 %heavy fabric
        temp = 1 - handles.heavyFabric;
    case 4 %light fabric
        temp = 1 - handles.lightFabric;
    case 5 %glass
        temp = 1 - handles.glass;
    case 6 %plaster
        temp = 1 - handles.plaster;
    case 7 %tiles
        temp = 1 - handles.tiles;
    case 8 %carpet
        temp = 1 - handles.carpet;
end
handles.aR = temp;
guidata(hObject, handles);
calcIR_Callback(hObject, eventdata, handles)

function wR_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in wFr.
function wFr_Callback(hObject, eventdata, handles)
choice = get(handles.wFr, 'Value');
switch choice
    case 1 %marble
        temp = 1 - handles.marble;
    case 2 %wood
        temp = 1 - handles.wood;
    case 3 %heavy fabric
        temp = 1 - handles.heavyFabric;
    case 4 %light fabric
        temp = 1 - handles.lightFabric;
    case 5 %glass
        temp = 1 - handles.glass;
    case 6 %plaster
        temp = 1 - handles.plaster;
    case 7 %tiles
        temp = 1 - handles.tiles;
    case 8 %carpet
        temp = 1 - handles.carpet;
end
handles.aF = temp;
guidata(hObject, handles);
calcIR_Callback(hObject, eventdata, handles)

function wFr_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in wRe.
function wRe_Callback(hObject, eventdata, handles)
choice = get(handles.wRe, 'Value');
switch choice
    case 1 %marble
        temp = 1 - handles.marble;
    case 2 %wood
        temp = 1 - handles.wood;
    case 3 %heavy fabric
        temp = 1 - handles.heavyFabric;
    case 4 %light fabric
        temp = 1 - handles.lightFabric;
    case 5 %glass
        temp = 1 - handles.glass;
    case 6 %plaster
        temp = 1 - handles.plaster;
    case 7 %tiles
        temp = 1 - handles.tiles;
    case 8 %carpet
        temp = 1 - handles.carpet;
end
handles.aB = temp;
guidata(hObject, handles);
calcIR_Callback(hObject, eventdata, handles)

function wRe_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in ceil.
function ceil_Callback(hObject, eventdata, handles)
choice = get(handles.ceil, 'Value');
switch choice
    case 1 %marble
        temp = 1 - handles.marble;
    case 2 %wood
        temp = 1 - handles.wood;
    case 3 %heavy fabric
        temp = 1 - handles.heavyFabric;
    case 4 %light fabric
        temp = 1 - handles.lightFabric;
    case 5 %glass
        temp = 1 - handles.glass;
    case 6 %plaster
        temp = 1 - handles.plaster;
    case 7 %tiles
        temp = 1 - handles.tiles;
    case 8 %carpet
        temp = 1 - handles.carpet;
end
handles.aC = temp;
guidata(hObject, handles);
calcIR_Callback(hObject, eventdata, handles)

function ceil_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in floor.
function floor_Callback(hObject, eventdata, handles)
choice = get(handles.floor, 'Value');
switch choice
    case 1 %marble
        temp = 1 - handles.marble;
    case 2 %wood
        temp = 1 - handles.wood;
    case 3 %heavy fabric
        temp = 1 - handles.heavyFabric;
    case 4 %light fabric
        temp = 1 - handles.lightFabric;
    case 5 %glass
        temp = 1 - handles.glass;
    case 6 %plaster
        temp = 1 - handles.plaster;
    case 7 %tiles
        temp = 1 - handles.tiles;
    case 8 %carpet
        temp = 1 - handles.carpet;
end
handles.aG = temp;
guidata(hObject, handles);
calcIR_Callback(hObject, eventdata, handles)

function floor_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Push buttons

% --- Executes on button press in calcIR.
function calcIR_Callback(hObject, eventdata, handles)
clear handles.IR;
%parameters
Lx = get(handles.Lx, 'Value'); %length of room, m
Ly = get(handles.Ly, 'Value'); %width of room, m
Lz = get(handles.Lz, 'Value'); %height of room, m
a = get(handles.a, 'Value'); %listener x location, m
b = get(handles.b, 'Value'); %listener y location, m
c = get(handles.c, 'Value'); %listener z location, m
p = get(handles.p, 'Value'); %source x location, m
q = get(handles.q, 'Value'); %source y location, m
r = get(handles.r, 'Value'); %source z location, m
alphaXpos = handles.aL;
alphaXneg = handles.aR;
alphaYpos = handles.aF;
alphaYneg = handles.aB;
alphaZpos = handles.aC;
alphaZneg = handles.aG; % room reflection coefficient, 0-1
cs = 343; %speed of sound m/s
Fs = handles.Fs; %sample rate
N = get(handles.N, 'Value'); %max image source order
dE = 0.2; %distance between ears

%create IR matrix
simtime = (1/cs) * sqrt((N*Lx + p - a)^2 + (N*Ly + q - b)^2 + (N*Lz + r - c)^2);
ir_R = zeros(1, ceil(1.2*simtime*Fs));
ir_L = zeros(1, ceil(1.2*simtime*Fs));

%listener positions
angle = atan((p-a)/(q-b));
aRight = a - ((dE/2)*cos(angle));
aLeft = a + ((dE/2)*cos(angle));
bRight = b + ((dE/2)*sin(angle));
bLeft = b - ((dE/2)*sin(angle));

for d = -N : N %x dimension
    if mod(d,2) == 1
        Aleft = (d+1)*Lx - p - aLeft;
    else
        Aleft = d*Lx + p - aLeft;
    end
    if mod(d,2) == 1
        Aright = (d+1)*Lx - p - aRight;
    else
        Aright = d*Lx + p - aRight;
    end
    [XNeg, XPos] = wallHits(d);
    for e = -N : N %y dimension
        if mod(e,2) == 1
            Bleft = (e+1)*Ly - q - bLeft;
        else
            Bleft = e*Ly + q - bLeft;
        end
        if mod(e,2) == 1
            Bright = (e+1)*Ly - q - bRight;
        else
            Bright = e*Ly + q - bRight;
        end
        [YNeg, YPos] = wallHits(e);
            for f = -N : N %z dimension
                if mod(f,2) == 1
                    C = (f+1)*Lz - r - c;
                else
                    C = f*Lz + r - c;
                end
                [ZNeg, ZPos] = wallHits(f);
                gA = (alphaXneg^(XNeg))*(alphaXpos^(XPos));
                gB = (alphaYneg^(YNeg))*(alphaYpos^(YPos));
                gC = (alphaZneg^(ZNeg))*(alphaZpos^(ZPos));
                Lleft = sqrt(Aleft^2 + Bleft^2 + C^2);
                tLeft = Lleft/cs;
                gLeft = (gA*gB*gC)/Lleft;
                Lright = sqrt(Aright^2 + Bright^2 + C^2);
                tRight = Lright/cs;
                gRight = (gA*gB*gC)/Lright;
                
                indL = ceil(tLeft*Fs);
                ir_L(1, indL) = ir_L(1, indL) + gLeft;
                indR = ceil(tRight*Fs);
                ir_R(1, indR) = ir_R(1, indR) + gRight;
            end
    end
end

if length(ir_L) > length(ir_R)
    ir_L = ir_L(1 : length(ir_R));
else
    ir_R = ir_R(1 : length(ir_L));
end

%plot IR
t = linspace(0, length(ir_L)/handles.Fs, length(ir_L)); %time vector
plot(t, ir_L, 'b', t, ir_R, 'r'), axis tight, grid on;
xlabel('Time (s)'), ylabel('Amplitude');

%HF attenuation
[B1, A1] = butter(get(handles.butterOrder, 'Value'), get(handles.butterCutoff,'Value')/(handles.Fs/2), 'low');
ir_L = filter(B1, A1, ir_L);
ir_R = filter(B1, A1, ir_R);

handles.IR = [ir_L; ir_R]; %store IR
guidata(hObject, handles);

% --- Executes on button press in playIR.
function playIR_Callback(hObject, eventdata, handles)
IRplayer = audioplayer(handles.IR, handles.Fs);
playblocking(IRplayer);

% --- Executes on button press in saveIR.
function saveIR_Callback(hObject, eventdata, handles)
filename = 'myIR.wav';
audiowrite(filename, handles.IR, handles.Fs);

% --- Executes on button press in loadAudio.
function loadAudio_Callback(hObject, eventdata, handles)
[handles.dry, dryFs] = getFile('the audio sample');
if ~isequal(dryFs, handles.Fs)
    handles.dry = resample(handles.dry, handles.Fs, dryFs);
end
if size(handles.dry, 1) > 1
    handles.dry = handles.dry(1,:);
end
set([handles.playDry handles.playWet], 'Enable', 'On');
guidata(hObject, handles);

% --- Executes on button press in playDry.
function playDry_Callback(hObject, eventdata, handles)
dryPlayer = audioplayer(handles.dry, handles.Fs);
playblocking(dryPlayer);

% --- Executes on button press in playWet.
function playWet_Callback(hObject, eventdata, handles)
% convolve dry signal with IR
wet_R = conv(handles.dry, handles.IR(1, :));
wet_L = conv(handles.dry, handles.IR(2, :));
handles.wet = [wet_L; wet_R];
% playback wet audio
player1 = audioplayer(handles.wet, handles.Fs);
playblocking(player1);
set(handles.saveWet, 'Enable', 'On');
guidata(hObject, handles);

% --- Executes on button press in saveWet.
function saveWet_Callback(hObject, eventdata, handles)
filename = 'myProcessedFile.wav';
audiowrite(filename, handles.wet, handles.Fs);

%plot room layout
function plotter_Callback(hObject, eventdata, handles)
axes(handles.axes1), set(handles.axes1, 'Visible', 'On');
cla(handles.axes1), reset(handles.axes1), hold off;
plot(get(handles.p, 'Value'), get(handles.q, 'Value'), 'k', 'Marker', 'o', 'MarkerSize', 10), hold on;
plot(get(handles.a, 'Value'), get(handles.b, 'Value'), 'k', 'Marker', 'x', 'MarkerSize', 10), hold on;
grid on, axis([0 get(handles.Lx, 'Value') 0 get(handles.Ly, 'Value')]);
title('o = source, x = listener');
axes(handles.axes2), set(handles.axes2, 'Visible', 'On');
cla(handles.axes2), reset(handles.axes2), hold off;
xlabel('Time (s)'), ylabel('Amplitude');
guidata(hObject, handles);
