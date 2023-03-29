% specifics the initial configuration of the cortical circuit
% This code sets up vectors of parameters for the model.

params.N_params=40;

% Time axis
params.tmin=0;
params.tmax=1;% s  % run time (seconds)
params.Fs=2048;%Hz  % sample rate (samples/second)
Fs=params.Fs;
params.dt=1/Fs; %  time step 

% external input
% ~~~~~~~~~~~~~~~~~
%params.noise_var = 0.574;
params.noise_var =0.1e-1;
%params.noise_var =1.8e-1;
params.He=1e-3;
% column parameters
% ~~~~~~~~~~~~~~~~~
params.input_index1=[1:22];% indexes for synapses that receive noise
%params.input_index2=[22+2:22+3,22+6:22+7,22+10:22+11,22+14:22+15   ,22+18:22+19,22+23:22+40];% indexes for synapses that receive input from other cortical areas and thalamic input
%params.input_index2=[22+3,22+14:22+15,22+18,22+21,22+24,22+28,22+33,22+35:22+40];% indexes for synapses that receive input from other cortical areas and thalamic input
params.input_index2=[22+2:22+3,22+6:22+8,22+10:22+11,22+14:22+15,22+18:22+19,22+21:22+22+7,22+22+9,22+22+11:22+22+14];% indexes for synapses that receive input from other cortical areas and thalamic input

% synapse parameters
% ~~~~~~~~~~~~~~~~~~
params.tau_e=0.02;



params.ConnectivityConst=0.5e-4;%connectivity constant

params.c=0.35e-9;% membrane capacitance

% Input
% ~~~~~~~~~~~~~~~~~~

params.Ip_o_inh=1; %Input from othere cortical areas to inhibitory neurons
params.Ip_o_exc=1; %Input from othere cortical areas to inhibitory neurons

params.Ip_t_inh=1;  %Thalamic input to inhibitory neurons
params.Ip_t_exc=1;  %Thalamic input to inhibitory neurons

% connectivity weights
% ~~~~~~~~~~~~~~~~~~
params.wi=0.2;
params.we=0.1;

% maximum firing rate
% ~~~~~~~~~~~~~~~~~~
params.alpha_e=0.7;
params.alpha_i=1;