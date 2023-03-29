function [Y_new,H] = whole_column_matrix_euler(Y,Q)


%n=params.n; % index of the sample
Fs=2048;%Hz  % sample rate (samples/second)

dt=1/Fs; %  time step 

ConnectivityConst=0.5e-4;%connectivity constant

c=0.35e-9;% membrane capacitance


tau_e2=0.02; % membrane time constant for pyramidal neurons in layer 2/3
tau_i2=0.04; % membrane time constant for inhibitory neurons in layer 2/3

tau_e4=0.02; % membrane time constant for pyramidal neurons in layer 4
tau_i4=0.04; % membrane time constant for inhibitory neurons in layer 4

tau_e5=0.02; % membrane time constant for pyramidal neurons in layer 5
tau_i5=0.04; % membrane time constant for inhibitory neurons in layer 5


tau_ee=1*8.5e-3; % synaptic time constant for excitatory-to-excitatory synapses
tau_ie=1*4.3e-3; % synaptic time constant for inhibitory-to-excitatory synapses
tau_ei=1*13e-3; % synaptic time constant for excitatory-to-inhibitory synapses
tau_ii=1*9e-3; % synaptic time constant for inhibitory-to-inhibitory synapses



wi=0.2;% weight for the connection from external input to the inhibitory neurons
we=0.1;% weight for the connection from external input to the excitatory neurons

% connectivity weights between different types of neurons
Ce2e2=0.06*ConnectivityConst;
Ce2i2=0.35*ConnectivityConst;
Ce2e4=0.11*ConnectivityConst;
Ce2i4=0.29*ConnectivityConst;

Ci2i2=1.06*ConnectivityConst;
Ci2e2=0.88*ConnectivityConst;
Ci2e4=0.21*ConnectivityConst;
Ci2i4=0.59*ConnectivityConst;
Ci2e5=0.13*ConnectivityConst;

Ce4e4=0.2*ConnectivityConst;
Ce4i4=0.48*ConnectivityConst;

Ci4e4=1.14*ConnectivityConst;
Ci4i4=1.47*ConnectivityConst;
Ci4e2=0.25*ConnectivityConst;
Ci4i2=0.18*ConnectivityConst;
Ci4e5=0.13*ConnectivityConst;
Ci4i5=0.11*ConnectivityConst;

Ce5e5=0.087*ConnectivityConst;
Ce5i5=0.49*ConnectivityConst;

Ci5e5=0.27*ConnectivityConst;
Ci5i5=1.51*ConnectivityConst;
Ci5e2=0.16*ConnectivityConst;
Ci5e4=0.24*ConnectivityConst;
Ci5i4=0.17*ConnectivityConst;

Ci=(wi)*ConnectivityConst;
Ce=(we)*ConnectivityConst;

% sigmoid parameters
v0=6;% mV
e0=15;%1/s
r=0.3; %1/mV

% Input
% ~~~~~~~~~~~~~~~~~~

Ip_o_inh=1; %Input from othere cortical areas to inhibitory neurons
Ip_o_exc=1; %Input from othere cortical areas to inhibitory neurons

Ip_t_inh=1;  %Thalamic input to inhibitory neurons
Ip_t_exc=1;  %Thalamic input to inhibitory neurons

% maximum firing rate
% ~~~~~~~~~~~~~~~~~~
alpha_e=0.7;
alpha_i=1;



N_states=22;
N_inputs=4;
N_beta=36+N_inputs;

Theta0=zeros(length(Y),1);
Theta0(N_states+2:N_states+N_beta)=ones(39,1);

decay_rate=0.0005;

CC=zeros(N_states+N_beta,N_states+N_beta);
RTAU=zeros(N_states+N_beta,1);

H=zeros(1,N_states+N_beta);


% a=[ repmat([1 0 0 0 0 0],5,1);...
%     repmat([0 1 0 0 0 0],5,1);...
%     repmat([0 0 1 0 0 0],7,1);...
%     repmat([0 0 0 1 0 0],7,1);...
%     repmat([0 0 0 0 1 0],5,1);...
%     repmat([0 0 0 0 0 1],4,1);...
%     ];
% A = [zeros(33,33)  a zeros(33,33); zeros(6+33,33+6+33)];
% B = [zeros(33,33+6) eye(33,33);zeros(33+6,66+6)];


syn_index = 0;
%%%%%%%%%%% Pyramidal layer 2/3
syn_index = syn_index + 1; %Synapse 1
RTAU(syn_index)=1/tau_e2;
presyn_inputs =[2 4 6 11 14];
CC(syn_index,presyn_inputs) = 1;
CC(syn_index,2)=Y(N_states+2)*Ce2e2;
CC(syn_index,4)=Y(N_states+4)*Ce;
CC(syn_index,6)=-Y(N_states+6)*Ce2i2;
CC(syn_index,11)=Y(N_states+11)*Ce2e4;
CC(syn_index,14)=-Y(N_states+14)*Ce2i4;



%%currents

syn_index = syn_index + 1;
RTAU(syn_index)=1/tau_ee;
H(syn_index)=100;

syn_index = syn_index + 1;
RTAU(syn_index)=1/tau_ie;



syn_index = syn_index + 1; %synapse 4 (Input from other cortical areas)
RTAU(syn_index)=1/tau_ee;
H(syn_index)=200;

%%%%%%%%%%% Inhibitory layer 2/3
syn_index = syn_index + 1; %Synapse 5
RTAU(syn_index)=1/tau_i2;
presyn_inputs =[3 7 8 10 15 18];
CC(syn_index,presyn_inputs) = 1; 
CC(syn_index,3)=Y(N_states+3)*Ci2e2;
CC(syn_index,7)=-Y(N_states+7)*Ci2i2;
CC(syn_index,8)=Y(N_states+8)*Ci;
CC(syn_index,10)=Y(N_states+10)*Ci2e4;
CC(syn_index,15)=-Y(N_states+15)*Ci2i4;
CC(syn_index,18)=Y(N_states+18)*Ci2e5;

%%currents
syn_index = syn_index + 1;
RTAU(syn_index)=1/tau_ei;
H(syn_index)=-100;

syn_index = syn_index + 1;
RTAU(syn_index)=1/tau_ii;


syn_index = syn_index + 1; % synapse 8 (Input from other cortical areas)
RTAU(syn_index)=1/tau_ie;


%%%%%%%%%%% Pyramidal layer 4
syn_index = syn_index + 1; %Synapse 9
RTAU(syn_index)=1/tau_e4;
presyn_inputs =[11 12 14];
CC(syn_index,presyn_inputs) = 1; 
CC(syn_index,11)=Y(N_states+22+1)*Ce4e4;
CC(syn_index,12)=Y(N_states+12)*Ce;
CC(syn_index,14)=-Y(N_states+22+2)*Ce4i4;
H(syn_index)=0;
%%currents
syn_index = syn_index + 1;
RTAU(syn_index)=1/tau_ie;


syn_index = syn_index + 1;
RTAU(syn_index)=1/tau_ee;
H(syn_index)=100;


syn_index = syn_index + 1; %Synpase 12 (Input from thalamus)
RTAU(syn_index)=1/tau_ee;
H(syn_index)=100;
%%%%%%%%%%% Inhibitory layer 4
syn_index = syn_index + 1; %Synapse 13 
RTAU(syn_index)=1/tau_i4;
presyn_inputs =[3 7 10 15 16 18 22];
CC(syn_index,presyn_inputs) = 1; 
CC(syn_index,3)=Y(N_states+22+3)*Ci4e2;
CC(syn_index,7)=-Y(N_states+22+4)*Ci4i2;
CC(syn_index,10)=Y(N_states+22+5)*Ci4e4;
CC(syn_index,15)=-Y(N_states+22+6)*Ci4i4;
CC(syn_index,16)=Y(N_states+16)*Ci;
CC(syn_index,18)=Y(N_states+22+7)*Ci4e5;
CC(syn_index,22)=-Y(N_states+22)*Ci4i5;


%%currents
syn_index = syn_index + 1;
RTAU(syn_index)=1/tau_ei;
H(syn_index)=-100;

syn_index = syn_index + 1;
RTAU(syn_index)=1/tau_ii;



syn_index = syn_index + 1; % synapse 16 (Input from thalamus)
RTAU(syn_index)=1/tau_ie;

%%%%%%%%%%% Pyramidal layer 5
syn_index = syn_index + 1; %Synapse 17 
RTAU(syn_index)=1/tau_e5;
presyn_inputs =[4 19 21];
CC(syn_index,presyn_inputs) = 1;
CC(syn_index,4)=Y(N_states+22+8)*Ce;
CC(syn_index,19)=Y(N_states+19)*Ce5e5;
CC(syn_index,21)=-Y(N_states+21)*Ce5i5;

H(syn_index)=0;

%%currents
syn_index = syn_index + 1;
RTAU(syn_index)=1/tau_ie;


syn_index = syn_index + 1;
RTAU(syn_index)=1/tau_ee;
H(syn_index)=100;


%%%%%%%%%%% Inhibitory layer 5
syn_index = syn_index + 1; %Synapse 20 [2 4 6 11 14][3 7 8 10 15 18][3 7 10 15 16 18 22][4 19 21];;;
RTAU(syn_index)=1/tau_i5;
presyn_inputs =[3 8 10 15 18 22];
CC(syn_index,presyn_inputs) = 1;
CC(syn_index,3)=Y(N_states+22+9)*Ci5e2;
CC(syn_index,8)=Y(N_states+22+10)*Ci;
CC(syn_index,10)=Y(N_states+22+11)*Ci5e4;
CC(syn_index,15)=-Y(N_states+22+12)*Ci5i4;
CC(syn_index,18)=Y(N_states+22+13)*Ci5e5;
CC(syn_index,22)=-Y(N_states+22+14)*Ci5i5;

%%currents
syn_index = syn_index + 1;
RTAU(syn_index)=1/tau_ei;
H(syn_index)=-100;

syn_index = syn_index + 1;
RTAU(syn_index)=1/tau_ii;

% model noise 
W=zeros(N_states+N_beta,1);% a vector for model noise
w = mvnrnd(0,Q(1),62);% noise for the input
W(1:N_states+N_beta)=w;

SS=[repmat(zeros(1,1),1,1);...
       repmat(alpha_e.*S(Y(1),v0,e0,r),2,1);...
    repmat(Y(N_states+36+1)*Ip_o_exc,1,1);...
    repmat(zeros(1,1),1,1);...
    repmat(alpha_i.*S(Y(5),v0,e0,r),2,1);...
    repmat(Y(N_states+36+2)*Ip_o_inh,1,1);...
    repmat(zeros(1,1),1,1);...
    repmat(alpha_e.*S(Y(9),v0,e0,r),2,1);...
    repmat(Y(N_states+36+3)*Ip_t_exc,1,1);...
    repmat(zeros(1,1),1,1);...
    repmat(alpha_i.*S(Y(13),v0,e0,r),2,1);...
    repmat(Y(N_states+36+4)*Ip_t_inh,1,1);...
    repmat(zeros(1,1),1,1);...
    repmat(alpha_e.*S(Y(17),v0,e0,r),2,1);...
    repmat(zeros(1,1),1,1);...
    repmat(alpha_i.*S(Y(20),v0,e0,r),2,1);...
    repmat(zeros(N_beta,1),1,1)];


Y_new=Y+(-RTAU.*Y+1/c.*CC*Y+SS).*dt+W+dt*decay_rate.*(Y.*Theta0(:,1)-Theta0(:,1));

function out=S(v,v0,e0,r)
    out=2*(e0)/(1+exp(r*(v0-v)));
end

end