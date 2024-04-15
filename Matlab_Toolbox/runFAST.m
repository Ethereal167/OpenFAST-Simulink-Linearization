%% runFAST.m
% This script will run a single FAST simulation with the ROSCO Simulink
% controller.  A more in depth version can be found at https://github.com/dzalkind/matlab-toolbox/tree/master/Simulations
clc
clear;
load("StateSpaceList.mat")

% 以下量不会随着工作点的切换而变换
Azimuth_list = StateSpaceList.Azimuth_list;
n_FixFrameInputs = StateSpaceList.n_FixFrameInputs;
n_RotTripletInputs = StateSpaceList.n_RotTripletInputs;
new_seq_inp = StateSpaceList.new_seq_inp;
n_FixFrameOutputs = StateSpaceList.n_FixFrameOutputs;
n_RotTripletOutputs = StateSpaceList.n_RotTripletOutputs;
new_seq_out = StateSpaceList.new_seq_out;
n_FixFrameStates2 = StateSpaceList.n_FixFrameStates2;
n_RotTripletStates2 = StateSpaceList.n_RotTripletStates2;
new_seq_states = StateSpaceList.new_seq_states;
windspeed_list = StateSpaceList.windspeed_list;

% 会随着工作点切换而变化的值
A_total = StateSpaceList.A_total;
B_total = StateSpaceList.B_total;
Bd_total = StateSpaceList.Bd_total;
C_total = StateSpaceList.C_total;
D_total = StateSpaceList.D_total;
Dd_total = StateSpaceList.Dd_total;
yop_List_total = StateSpaceList.yop_List_total;
xop_total = StateSpaceList.xop_total;
uop_total = StateSpaceList.uop_total;
wop_total = StateSpaceList.wop_total;


[this_dir,~,~] = fileparts(mfilename('fullpath'));

% Compile FAST for use with simulink & mex using openfast docs
fast.FAST_SFuncDir     = 'D:\MatlabProject\OpenFast_v2.4.0-v3.0.0\Openfast-Simulink联合仿真模型\NREL5MW_ROSCO_IPC\ROSCO_File\ROSCO-2.7.0\Test_Cases\simulink';  %%%% NEED FOR SIMULINK
fast.FAST_InputFile    = '5MW_Land_Simulink.fst';   % FAST input file (ext=.fst)
fast.FAST_directory    = 'D:\MatlabProject\OpenFast_v2.4.0-v3.0.0\Openfast-Simulink联合仿真模型\NREL5MW_ROSCO_IPC\ROSCO_File\ROSCO-2.7.0\Test_Cases\5MW_Land_Simulink';   % Path to fst directory files

% Simulink Parameters
% Model
simu.SimModel           = fullfile(this_dir,'Simulink','ROSCO');

% Script for loading parameters
simu.ParamScript        = fullfile(this_dir,'Utilities','load_ROSCO_params');

%% Simulink Setup

[ControlScriptPath,ControScript] = fileparts(simu.ParamScript);
addpath(ControlScriptPath);
addpath(fast.FAST_SFuncDir);
addpath('Utilities')
%% Read FAST Files & Load ROSCO Params from DISCON.IN

[Param,Cx] = ReadWrite_FAST(fast);

% Simulation Parameters
simu.TMax   = 100;   % 仿真时长
simu.dt   = Param.FP.Val{contains(Param.FP.Label,'DT')};
[R,F] = feval(ControScript,Param,simu);
% IPC参数
R.IPC_KP_1P = 1*1e-8;
R.IPC_KP_2P = 1*1e-8;
R.IPC_KI_1P = 1e-5;
R.IPC_KI_2P = 1e-5;
R.IPC_ControlMode = 2;

%% Premake OutList for Simulink

OutList = {'Time'};
OutList = [OutList;
    Param.IWP.OutList;
    Param.EDP.OutList;
    Param.ADP.OutList;
    Param.SvDP.OutList;
    ];

OutList = vertcat(OutList{1:5});
OutList = erase(OutList, """");


%% Exectute FAST

% Using Simulink/S_Func
FAST_InputFileName = [fast.FAST_directory,filesep,fast.FAST_InputFile];
TMax               = simu.TMax;

SimulinkModel = simu.SimModel;
open([SimulinkModel '.slx'])

% Out         = sim(SimulinkModel, 'StopTime', num2str(simu.TMax));
% sigsOut     = get(Out,'sigsOut');   %internal control signals
% 
% %% Get OutData
% 
% SFuncOutStr = '.SFunc';
% 
% % Try text first, then binary
% [OutData,OutList] = ReadFASTtext([fast.FAST_directory,filesep,fast.FAST_InputFile(1:end-4),SFuncOutStr,'.out']);
% if isempty(OutData)
%     [OutData,OutList] = ReadFASTbinary([fast.FAST_directory,filesep,fast.FAST_InputFile(1:end-4),SFuncOutStr,'.outb']);
% end
% 
% % Dump data to structure
% for i = 1:length(OutList)
%     simout.(OutList{i}) = OutData(:,i);
% end
% 
% %% Plot
% Pl_FastPlots(simout)
