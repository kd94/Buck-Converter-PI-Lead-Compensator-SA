function varargout = BuckConverterControllerDesign(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @BuckConverterControllerDesign_OpeningFcn, ...
                   'gui_OutputFcn',  @BuckConverterControllerDesign_OutputFcn, ...
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

function BuckConverterControllerDesign_OpeningFcn(hObject, eventdata, handles, varargin)

handles.output = hObject;

guidata(hObject, handles);

function varargout = BuckConverterControllerDesign_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;

function Vin_Text_Callback(hObject, eventdata, handles)

function Vin_Text_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Vin_EditBox_Callback(hObject, eventdata, handles)

function Vout_Text_Callback(hObject, eventdata, handles)

function Vout_Text_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Vout_EditBox_Callback(hObject, eventdata, handles)

function Pout_Text_Callback(hObject, eventdata, handles)

function Pout_Text_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Pout_EditBox_Callback(hObject, eventdata, handles)

function DVout_Text_Callback(hObject, eventdata, handles)

function DVout_Text_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function DVout_EditBox_Callback(hObject, eventdata, handles)

function DIL_Text_Callback(hObject, eventdata, handles)

function DIL_Text_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function DIL_EditBox_Callback(hObject, eventdata, handles)

function fsw_Text_Callback(hObject, eventdata, handles)

function fsw_Text_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function fsw_EditBox_Callback(hObject, eventdata, handles)

function L_EditBox_Callback(hObject, eventdata, handles)

function L_EditBox_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Calc_L_Callback(hObject, eventdata, handles)

    Vin  = str2double(get(handles.Vin_EditBox, 'string'));
    Vout = str2double(get(handles.Vout_EditBox, 'string'));
    Pout = str2double(get(handles.Pout_EditBox, 'string'));
    Percent_DIL = (str2double(get(handles.DIL_EditBox, 'string')))/100;
    fsw = str2double(get(handles.fsw_EditBox, 'string'));  
    Tsw = 1 / fsw;
    Iout = Pout / Vout;
    IL = Iout;
    DIL = Percent_DIL * IL;
    D = Vout / Vin;    
    L_Calc = (Vout * (1-D)) / (DIL * fsw);
    set(handles.L_EditBox, 'string',L_Calc);
    msgbox(sprintf('Is The L Value(%f) ok? If Not Please Change It And Click The (Calculate C) Button!',L_Calc)); 
    IL_min = IL - (DIL/2);
    IL_max = IL + (DIL/2);
    save('param.mat','Vin','Vout','Pout','Percent_DIL','fsw','Tsw','DIL','D','IL_min','IL_max');
    guidata(hObject, handles);

function C_EditBox_Callback(hObject, eventdata, handles)

function C_EditBox_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Calc_C_Callback(hObject, eventdata, handles)
    
    load('param.mat');
    Percent_DVout=((str2double(get(handles.DVout_EditBox, 'string'))))/100; %Load transient response %X of Vout
    DVo = Percent_DVout*Vout;
    L_Select = str2double(get(handles.L_EditBox, 'string'));
    C_Calc = (Vout*(1-D))/(8*L_Select*DVo*fsw*fsw);
    set(handles.C_EditBox, 'string',C_Calc);   
    save('param.mat','Percent_DVout','DVo','-append');
    guidata(hObject, handles);

function Bode_TF_BUCK_Callback(hObject, eventdata, handles)
    
    load('param.mat');
    L = str2double(get(handles.L_EditBox, 'string'));
    C = str2double(get(handles.C_EditBox, 'string'));
    rl = str2double(get(handles.rL_EditBox, 'string'));
    rc = str2double(get(handles.rc_EditBox, 'string'));
    rsw = str2double(get(handles.rsw_EditBox, 'string'));
    rd = str2double(get(handles.rd_EditBox, 'string'));
    Vm = 1; %Sawtooth peak

    rx = (D*rsw)+((1-D)*rd);
    Iout = Pout/Vout;
    Rout = Vout/Iout;
    IL = Iout;

    Gvd_num = [(Rout*(((rd-rsw)*IL)+Vin)*rc*C)/((Rout+rc)*L*C) (Rout*(((rd-rsw)*IL)+Vin))/((Rout+rc)*L*C)];
    Gvd_denum = [1 (L+(C*((rc*Rout)+((rx+rl)*(Rout+rc)))))/((Rout+rc)*L*C) (Rout+rx+rl)/((Rout+rc)*L*C)];
    Gvd = tf(Gvd_num,Gvd_denum);
    Gpwm = 1/Vm;
    tf_buck = Gpwm*Gvd

    save('param.mat','L','C','rl','rc','rsw','rd','Iout','Rout','tf_buck','-append');

    bodeplot(handles.axes_bode,tf_buck,{0.1,1e+10})
    grid on
    hold on
    
    [Gm,Pm,Wcg,Wcp] = margin(tf_buck)
    set(handles.GainMargin_EditBox, 'string',Gm);
    set(handles.PhaseMargin_EditBox, 'string',Pm);
    set(handles.PM_freq_EditBox, 'string',Wcp);

    tf_buck_step=feedback(tf_buck,1)
    step(handles.axes_StepResponse,tf_buck_step)

    tf_buck_stepinfo = stepinfo(tf_buck_step);
    set(handles.Overshoot_EditBox, 'string', tf_buck_stepinfo.Overshoot);
    set(handles.Peaktime_EditBox, 'string', tf_buck_stepinfo.PeakTime);
    set(handles.Settlingtime_EditBox, 'string', tf_buck_stepinfo.SettlingTime);
    set(handles.Risetime_EditBox, 'string', tf_buck_stepinfo.RiseTime);
    grid on
    hold on
    guidata(hObject, handles);
    
function GainMargin_EditBox_Callback(hObject, eventdata, handles)

function GainMargin_EditBox_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Overshoot_EditBox_Callback(hObject, eventdata, handles)

function Overshoot_EditBox_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function PhaseMargin_EditBox_Callback(hObject, eventdata, handles)

function PhaseMargin_EditBox_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function PM_freq_EditBox_Callback(hObject, eventdata, handles)

function PM_freq_EditBox_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Risetime_EditBox_Callback(hObject, eventdata, handles)

function Risetime_EditBox_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Settlingtime_EditBox_Callback(hObject, eventdata, handles)

function Settlingtime_EditBox_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Peaktime_EditBox_Callback(hObject, eventdata, handles)

function Peaktime_EditBox_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function G_fsw_EditBox_Callback(hObject, eventdata, handles)

function G_fsw_EditBox_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function G_PI_Text_Callback(hObject, eventdata, handles)

function G_PI_Text_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function G_PI_EditBox_Callback(hObject, eventdata, handles)

function G_PI_EditBox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function G_Lead_Text_Callback(hObject, eventdata, handles)

function G_Lead_Text_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Add_PILead_Comp_Callback(hObject, eventdata, handles)
    StepNum=str2double(char(inputdlg({'ENTER THE NUMBER OF ITERATIONS'})));
    load('param.mat');

    [optimum_fz, optimum_fc, optimum_pm] = PILeadOptimizer(StepNum);
    
    set(handles.fzero_EditBox, 'string',optimum_fz);
    set(handles.fcrossover_EditBox, 'string',optimum_fc);
    set(handles.PM_EditBox, 'string',optimum_pm);

    optimum_wz = 2*pi*optimum_fz;
    Gpi_num = [(1/optimum_wz) 1];
    Gpi_denum = [1 0];
    Gpi = tf(Gpi_num,Gpi_denum);
  
    Gpi_string = evalc('Gpi');
    set(handles.G_PI_EditBox, 'string', Gpi_string);

    G1 = Gpi*tf_buck;
    G1_close_loop = feedback(G1,1);
    bodeplot(handles.axes_bode,G1,{0.1,1e+10})
    step(handles.axes_StepResponse,G1_close_loop)
    grid on
    hold on
    
%     fgc_max = fsw/4;
%     fgc_min = fsw/10;
    optimum_wc = 2*pi*optimum_fc;

    [GM, PM, Wcg, Wcp] = margin(G1); %wcg Magnitude eğrisinde 0dB'ye karşılık gelen frekans (PM'deki), wcp Phase eğrisinde -180 dereceye karşılık gelen frekans (GM'deki) 
    [mag, phase, wout] = bode(G1, {0.1,1e+15});
    mag = squeeze(mag);
    phase = squeeze(phase);
    wout = squeeze(wout);

    G_K1 = interp1(wout,20*log10(mag),optimum_wc); % at Desired_wgc
    PM1 = interp1(wout,phase,optimum_wc); % at Desired_wgc

    K1 = 10^(G_K1/20);
    Kreq = 1/K1;
    PMreq = -180-PM1+optimum_pm;
    K = Kreq*sqrt((1+(sind(PMreq)))/(1-(sind(PMreq))));
    Alpa = optimum_wc*sqrt((1-(sind(PMreq)))/(1+(sind(PMreq))));
    Beta = optimum_wc*sqrt((1+(sind(PMreq)))/(1-(sind(PMreq))));

    G_Lead_num = [K (K*Alpa)];
    G_Lead_denum = [1 Beta];
    G_Lead = tf(G_Lead_num,G_Lead_denum);
    G_Lead_string = evalc('G_Lead');
    set(handles.G_Lead_EditBox, 'string', G_Lead_string);

    Gpi_lead = G_Lead*Gpi*tf_buck;
    Gpi_lead_CloseLoop = feedback(Gpi_lead, 1);
    bodeplot(handles.axes_bode,Gpi_lead,{0.1,1e+10})
    step(handles.axes_StepResponse,Gpi_lead_CloseLoop)
    grid on
    hold on

    [Gm,Pm,Wcg,Wcp] = margin(Gpi_lead);
    set(handles.GainMargin_EditBox, 'string',Gm);
    set(handles.PhaseMargin_EditBox, 'string',Pm);
    set(handles.PM_freq_EditBox, 'string',Wcp);
     
    Gpi_lead_stepinfo = stepinfo(Gpi_lead_CloseLoop);
    set(handles.Overshoot_EditBox, 'string',Gpi_lead_stepinfo.Overshoot);
    set(handles.Peaktime_EditBox, 'string',Gpi_lead_stepinfo.PeakTime);
    set(handles.Settlingtime_EditBox, 'string',Gpi_lead_stepinfo.SettlingTime);
    set(handles.Risetime_EditBox, 'string',Gpi_lead_stepinfo.RiseTime);

function Clc_Axes_Callback(hObject, eventdata, handles)
cla(handles.axes_bode,'reset')
cla(handles.axes_StepResponse,'reset')

function Clc_Data_Callback(hObject, eventdata, handles)
delete('param.mat');
set([handles.Vin_EditBox,...
    handles.Vout_EditBox,...
    handles.Pout_EditBox,...
    handles.DVout_EditBox,...
    handles.DIL_EditBox,...
    handles.fsw_EditBox,...
    handles.L_EditBox,...
    handles.C_EditBox], 'String','');
clear all;
clc;

function fzero_Text_Callback(hObject, eventdata, handles)

function fzero_Text_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function fzero_EditBox_Callback(hObject, eventdata, handles)

function fzero_EditBox_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function fcrossover_Text_Callback(hObject, eventdata, handles)

function fcrossover_Text_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function fcrossover_EditBox_Callback(hObject, eventdata, handles)

function fcrossover_EditBox_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function PM_Text_Callback(hObject, eventdata, handles)

function PM_Text_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function PM_EditBox_Callback(hObject, eventdata, handles)

function PM_EditBox_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function G_Lead_EditBox_Callback(hObject, eventdata, handles)

function G_Lead_EditBox_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function rc_EditBox_Callback(hObject, eventdata, handles)

function rc_EditBox_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function rc_Text_Callback(hObject, eventdata, handles)

function rc_Text_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function rL_EditBox_Callback(hObject, eventdata, handles)

function rL_EditBox_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function rL_Text_Callback(hObject, eventdata, handles)

function rL_Text_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function rd_EditBox_Callback(hObject, eventdata, handles)

function rd_EditBox_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function rd_Text_Callback(hObject, eventdata, handles)

function rd_Text_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function rsw_EditBox_Callback(hObject, eventdata, handles)

function rsw_EditBox_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function rsw_Text_Callback(hObject, eventdata, handles)

function rsw_Text_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
