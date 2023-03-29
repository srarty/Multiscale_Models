
clc 
clear 
close all
% Handle to model function
func_f  = @neural_mass_model;

Parameters % parameters of the neural mass model




N_states=22;
N_inputs=4;
N_Beta=36+N_inputs;

N_aug=N_states+N_Beta; % augmented state vector
input_index1=params.input_index1;
input_index2=params.input_index2;
%load EEG_CN2

% load EEG_T1Day6
% BACKGROUND{1}=Background_signal;
% load EEG_T1Day7
% BACKGROUND{2}=Background_signal;
% load EEG_T1Day9
% BACKGROUND{3}=Background_signal;
%load EEG_CN3
%EEG_CN={EEG_CN1_1,EEG_CN1_2,EEG_CN1_3,EEG_CN1_4,EEG_CN1_5,EEG_CN1_6,EEG_CN1_7,EEG_CN1_8,EEG_CN1_9,EEG_CN1_10};

for jj=1:40
    %load(['EEG_T1Day',num2str(jj) '.mat'])
        load(['C:\Users\pzarei\OneDrive - The University of Melbourne\MATLAB codes\2023\March\read Data\CD2\EEG_CD2Day',num2str(jj) '.mat'])

    %Background_signal=BACKGROUND{jj};
% load 10995_20151213.mat
% output=chRData(10800000:11200000)';
%Background_signal=BACKGROUND{jj};
parfor ii=1:20
EEG_cn=Background_signal{ii};

output=EEG_cn(2048*3:length(EEG_cn)-2048*3)';
if isempty(output)
    continue
end
trend_line = smoothdata(output,'sgolay',2*2048);
output1 = output - 1*trend_line;
figure
plot(output1)

DC_offset=5.78;




% % scale data
% new_min = -10;
% new_max = 10;

% % scale data for non seizure
new_min = -0.5;
new_max = 0.5;

maxEEG=max(output1);
minEEG=min(output1);
dataRange=maxEEG-minEEG;
output2=output1;
output1(isnan(output1))=0;
m=mean(output1);
output1=output2;
output1 = (new_max - new_min) * (output1 - m) ./ dataRange+DC_offset;

N_samples=length(output1);

tmax=floor(N_samples/2048); % s  % run time (seconds)
tmin=0;
dt=1/2048;
t=tmin:dt:tmax;% vector for all time points






% add noise to measurements
%R = 0.1;
%R = .0001;
R=0.001;
%R=0;
% obs_noise = sqrt(R).*randn(1,N_samples); %sqrt(R).*randn(N_obs,N_samples);
% output1= output1 + obs_noise;

% Define the variance of the noise input
Q = zeros(1,N_aug);
%params.noise_var=0;
sigma = 0.1e-1;
He = 1e-3;             %noise coefficient
taue = 0.02;         % excitatory time constant
input_index1=[1:22];% indexes for synapses that receive noise
input_index2=[22+2:22+3,22+6:22+8,22+10:22+11,22+14:22+15,22+18:22+19,22+22:22+22+7,22+22+9,22+22+11:22+22+14];% indexes for synapses that receive input from other cortical areas and thalamic input

Q(input_index1) = 0.1*(sqrt(dt)*sigma*He/taue)^2;
Q(input_index2) = 0.1*(sqrt(dt)*sigma*He/taue)^2;

y0=zeros(N_aug,1);
y0(N_states+1:N_states+N_Beta)=ones(N_Beta,1);
[~,H]=whole_column_matrix_euler(y0,Q);
% Process noise, Q
estQ = zeros(1,N_aug);
estQ(1:N_aug) = Q;        % use the known process noise variance
estQ = diag(estQ);



%% Apply Kalman Filter
%~~~~~~~~~~~~~~~~~~~~~

% FILTER PARAMETERS
alpha = 1;                  % between 0 and 1. Use 1 for large N_aug
beta = 2;                   % 2 is optimal for Gaussians
kappa = 3;                  % ***Don't change***** 3 or 4
mat = 0;                    % zero

y_hat = zeros(N_aug,N_samples);     % initialise augmented state vector
P_hat = zeros(N_aug,N_samples);     % initialise model covariance matrix (to store diagonals)

M = y_hat(:,1);   % first guess of states/params


M(N_states+1:end)=ones(1,N_Beta);
% M(N_states+2)=1;%*Ce2e2;
% M(N_states+4)=%Ce;
% M(N_states+6)=1;%Ce2i2;
% M(N_states+11)%Ce2e4;
% M(N_states+14)%Ce2i4;
% M(N_states+2)%Ce2e2;
% M(N_states+4)%Ce;
% M(N_states+6)%Ce2i2;
% M(N_states+11)%Ce2e4;
% M(N_states+14)%Ce2i4;
% 
% M(N_states+3)=0.88;%Ci2e2;
% M(N_states+7)=1;%Ci2i2;
% M(N_states+8)%Ci;
% M(N_states+10)%Ci2e4;
% M(N_states+15)%Ci2i4;
% M(N_states+18)%Ci2e5;
% 
% M(N_states+22+1)%Ce4e4;
% M(N_states+12)%Ce;
% M(N_states+22+2)%Ce4i4;
% 
% M(N_states+22+3)%Ci4e2;
% M(N_states+22+4)%Ci4i2;
% M(N_states+22+5)%Ci4e4;
% M(N_states+22+6)%Ci4i4;
% M(N_states+16)%Ci;
% M(N_states+22+7)%Ci4e5;
% M(N_states+22)%Ci4i5;
% 
% M(N_states+22+8)%Ce;
% M(N_states+19)%Ce5e5;
% M(N_states+21)%Ce5i5;
% 
% M(N_states+22+9)%Ci5e2;
% M(N_states+22+10)%Ci;
% M(N_states+22+11)%Ci5e4;
% M(N_states+22+12)%Ci5i4;
% M(N_states+22+13)%Ci5e5;
% M(N_states+22+14)%Ci5i5;
% 

y_hat(:,1)=M ;
% run the model forward
y_out=zeros(N_aug,N_samples);
for n=2:N_samples

    y_out(:,n) = func_f(y_out(:,n-1),Q);
    y_out(:,n)=y_out(:,n);
end

%%first guess of P for the states is variance of the simulation
y_var = var(y_out,0,2);
P0 = y_var;
P = diag(P0);     % first estimate of covariance
%%Apply bounds to augmented state vector
% *************************************

constraints =repmat([0 100],N_Beta,1);  
   
key =N_states+1:N_Beta+N_states;         
for n = 1:N_samples

    %apply filter
    if isnan(output1(:,n))
        y_hat(:,n)=NaN(62,1);
        continue
    end
    [M,P,S,C] = ukf_predict1(M,P,func_f,estQ,alpha,beta,kappa,mat,constraints,key);
    [M,P,K] = ukf_update1(M,P,output1(:,n),H,R,alpha,beta,kappa,mat,constraints,key);
    %store next estimate and covariance
    y_hat(:,n) = M;
    P_hat(:,n) = diag(P);

end

% save y_hat y_hat
 
y_output=H*y_hat;
% Plots
% figure
% plot(t,output1(1:tmax*Fs+1),'r--')
% hold on
% plot(t(1,10:end),y_output(1,10:tmax*Fs+1),'b')
% xlabel('time (s)')
% title(' EEG')
Y_HAT{ii}=y_hat;
Y_OUT{ii}=y_output;
OUTPUT1{ii}=output1;


% figure
% for i=1:8
%     subplot(8,1,i),plot(t(1,1000:end),y_hat(i,(1000:tmax*Fs+1)))
% end
% subplot(8,1,8), title('input to inhibitory neurons ')
% subplot(8,1,4),title('input to excitatory neurons ')
% subplot(8,1,1),title('membrane potential of excitatory neurons (2/3)')
% subplot(8,1,5),title('membrane potential of inhibitory neurons (2/3)')
% subplot(8,1,2),title('E-E current')
% subplot(8,1,3),title('E-I current')
% subplot(8,1,6),title('I-E current')
% subplot(8,1,7),title('I-I current')
% 
% 
% figure
% for i=1:8
%     subplot(8,1,i),plot(t(1000:end),y_hat(i+8,1000:tmax*Fs+1))
%     
% end
% subplot(8,1,8), title('thalamic input to inhibitory neurons ')
% subplot(8,1,4),title('thalamic input to excitatory neurons ')
% subplot(8,1,1),title('membrane potential of excitatory neurons (4)')
% subplot(8,1,5),title('membrane potential of inhibitory neurons (4)')
% subplot(8,1,2),title('E-I current')
% subplot(8,1,3),title('E-E current')
% subplot(8,1,6),title('I-E current')
% subplot(8,1,7),title('I-I current')
% 
% figure
% for i=1:6
%     subplot(6,1,i),plot(t(1000:end),y_hat(i+8+8,1000:tmax*Fs+1))
%     
% end
% subplot(6,1,1),title('membrane potential of excitatory neurons (5)')
% subplot(6,1,4),title('membrane potential of inhibitory neurons (5)')
% subplot(6,1,2),title('E-I current')
% subplot(6,1,3),title('E-E current')
% subplot(6,1,5),title('I-E current')
% subplot(6,1,6),title('I-I current')
% 
% 
% figure
% for i=1:10
%     subplot(10,1,i)
%     plot(t(1:end),y_hat(i+22,1:tmax*Fs+1))
%  
% end
% figure
% subplot(7,1,1),plot(t(1,1:end),output1(1,1:tmax*Fs+1),'b')
% title(' EEG')
% 
% 
% for i=2:7
%     subplot(7,1,i)
%     plot(t(1:end),y_hat(i+22+10,1:tmax*Fs+1))
% 
%     %     ylim([0 1.5])
% end
% 
% figure
% subplot(9,1,1),plot(t(1,1:end),output1(1,1:tmax*Fs+1),'b')
% title(' EEG')
% 
% 
% for i=2:9
%     subplot(9,1,i)
%     plot(t(1:end),y_hat(i+22+17,1:tmax*Fs+1))
% 
%     %     ylim([0 1.5])
% end
% %%
% figure
% subplot(7,1,1),plot(t(1,1:end),output1(1,1:tmax*Fs+1),'b')
% title(' EEG')
% 
% 
% subplot(7,1,2),plot(t(1,1:end),y_hat(2+22,1:tmax*Fs+1))
% title('Beta for excitatory to excitatory connections (e2e2)')
% 
% 
% subplot(7,1,3),plot(t(1,1:end),y_hat(3+22,1:tmax*Fs+1))
% title('Beta for excitatory to inhibitory connections (i2e2)')
% 
% 
% subplot(7,1,4),plot(t(1,1:end),y_hat(6+22,1:tmax*Fs+1))
% title('Beta for inhibitory to excitatory connections (e2i2)')
% 
% 
% subplot(7,1,5),plot(t(1,1:end),y_hat(7+22,1:tmax*Fs+1))
% title('Beta for inhibitory to inhibitory connections (i2i2)')
% 
% 
% subplot(7,1,6),plot(t,y_hat(22+36+1,1:tmax*Fs+1))
% title('Firing rate of external input to excitatory neurons (other)')
% 
% 
% subplot(7,1,7),plot(t,y_hat(22+36+2,1:tmax*Fs+1))
% title('Firing rate of external input to inhibitory neurons (other)')
% 
% 
% %%
% 
% 
% % a1=max([max(y_hat(5+22+22,1:end)),max(y_hat(1+22+22,1:end)),max(y_hat(2+22+22,1:end)),max(y_hat(6+22+22,1:end)),max(y_hat(22+36+3,:)),max(y_hat(22+36+4,:))]);
% % b1=min([min(y_hat(5+22+22,1:end)),min(y_hat(1+22+22,1:end)),min(y_hat(2+22+22,1:end)),min(y_hat(6+22+22,1:end)),min(y_hat(22+36+3,:)),min(y_hat(22+36+4,:))]);
% % ylim([b1 a1]);
% 
% figure
% subplot(6,1,1),plot(t(1,1:end),y_hat(5+22+22,1:tmax*Fs+1),'b')
% title('Beta for excitatory to inhibitory connections (i4e4)')
% 
% 
% subplot(6,1,2),plot(t(1,1:end),y_hat(1+22+22,1:tmax*Fs+1))
% title('Beta for excitatory to excitatory connections (e4e4)')
% 
% 
% subplot(6,1,3),plot(t(1,1:end),y_hat(2+22+22,1:tmax*Fs+1))
% title('Beta for inhibitory to excitatory connections (e4i4)')
% 
% subplot(6,1,4),plot(t(1,1:end),y_hat(6+22+22,1:tmax*Fs+1))
% title('Beta for inhibitory to inhibitory connections (i4i4)')
% 
% 
% subplot(6,1,5),plot(t,y_hat(22+36+3,1:tmax*Fs+1))
% title('Firing rate of thalamic input to excitatory neurons')
% 
% 
% subplot(6,1,6),plot(t,y_hat(22+36+4,1:tmax*Fs+1))
% title('Firing rate of thalamic input to inhibitory neurons')
% 
% 
% % a1=max(max(y_hat(5+22+22,1:end)),max(y_hat(1+22+22,1:end)),max(y_hat(2+22+22,1:end)),max(y_hat(6+22+22,1:end)),max(y_hat(22+36+3,:)),max(y_hat(22+36+4,:)));
% % b1=min(min(y_hat(5+22+22,1:end)),min(y_hat(1+22+22,1:end)),min(y_hat(2+22+22,1:end)),min(y_hat(6+22+22,1:end)),min(y_hat(22+36+3,:)),min(y_hat(22+36+4,:)));
% 
% %%
% figure
% subplot(4,1,1),plot(t(1,1:end),y_hat(13+22+22,1:tmax*Fs+1),'b')
% title('Beta for excitatory to inhibitory connections (i5e5)')
% 
% subplot(4,1,2),plot(t(1,1:end),y_hat(19+22,1:tmax*Fs+1))
% title('Beta for excitatory to excitatory connections (e5e5)')
% 
% 
% 
% 
% subplot(4,1,3),plot(t(1,1:end),y_hat(21+22,1:tmax*Fs+1))
% title('Beta for inhibitory to excitatory connections (e5i5)')
% 
% 
% subplot(4,1,4),plot(t(1,1:end),y_hat(14+22+22,1:tmax*Fs+1))
% title('Beta for inhibitory to inhibitory connections (i5i5)')
% 
% %%
% 
% 
% figure %Layer 2 and Layer 4
% subplot(6,1,1),plot(t(1,1:end),y_hat(10+22,1:tmax*Fs+1),'b')
% title('Beta for excitatory to inhibitory connections (i2e4)')
% 
% 
% subplot(6,1,2),plot(t(1,1:end),y_hat(11+22,1:tmax*Fs+1))
% title('Beta for excitatory to excitatory connections (e2e4)')
% 
% 
% subplot(6,1,3),plot(t(1,1:end),y_hat(14+22,1:tmax*Fs+1))
% title('Beta for inhibitory to excitatory connections (e2i4)')
% 
% 
% subplot(6,1,4),plot(t(1,1:end),y_hat(15+22,1:tmax*Fs+1))
% title('Beta for inhibitory to inhibitory connections (i2i4)')
% 
% 
% 
% subplot(6,1,5),plot(t(1,1:end),y_hat(3+22+22,1:tmax*Fs+1),'b')
% title('Beta for excitatory to inhibitory connections (i4e2)')
% 
% 
% subplot(6,1,6),plot(t(1,1:end),y_hat(4+22+22,1:tmax*Fs+1))
% title('Beta for inhibitory to inhibitory connections (i4i2)')
% 
% 
% 
% %%
% figure %Layer 5 and other layers
% subplot(6,1,1),plot(t(1,1:end),y_hat(18+22,1:tmax*Fs+1),'b')
% title('Beta for excitatory to inhibitory connections (i2e5)')
% 
% 
% subplot(6,1,2),plot(t(1,1:end),y_hat(22+7+22,1:tmax*Fs+1))
% title('Beta for excitatory to inhibitory connections (i4e5)')
% 
% 
% subplot(6,1,3),plot(t(1,1:end),y_hat(22+11+22,1:tmax*Fs+1))
% title('Beta for excitatory to inhibitory connections (i5e4)')
% 
% 
% subplot(6,1,4),plot(t(1,1:end),y_hat(12+22+22,1:tmax*Fs+1))
% title('Beta for inhibitory to inhibitory connections (i5i4)')
% 
% 
% 
% subplot(6,1,5),plot(t(1,1:end),y_hat(9+22+22,1:tmax*Fs+1),'b')
% title('Beta for excitatory to inhibitory connections (i5e2)')
% 
% 
% subplot(6,1,6),plot(t(1,1:end),y_hat(22+22,1:tmax*Fs+1))
% title('Beta for inhibitory to inhibitory connections (i4i5)')

end
save(['C:\Users\pzarei\OneDrive - The University of Melbourne\MATLAB codes\2023\March\read Data\CD    2\Y_HAT_CD2Day' num2str(jj) '.mat'],'Y_HAT', '-v7.3')
save(['C:\Users\pzarei\OneDrive - The University of Melbourne\MATLAB codes\2023\March\read Data\CD2\OUTPUT1_CD2Day' num2str(jj) '.mat'],'OUTPUT1')
save(['C:\Users\pzarei\OneDrive - The University of Melbourne\MATLAB codes\2023\March\read Data\CD2\Y_OUT_CD2Day' num2str(jj) '.mat'],'Y_OUT')
% if jj==1
%     save('Y_HAT_T1Day6.mat', 'Y_HAT', '-v7.3')
%     save OUTPUT1_T1Day6 OUTPUT1
%     save Y_OUT_T1Day6 Y_OUT
% elseif jj==2
%     save('Y_HAT_T1Day7.mat', 'Y_HAT', '-v7.3')
%     save OUTPUT1_T1Day7 OUTPUT1
%     save Y_OUT_T1Day7 Y_OUT
% elseif jj==3
%     save('Y_HAT_T1Day9.mat', 'Y_HAT', '-v7.3')
%     save OUTPUT1_T1Day9 OUTPUT1
%     save Y_OUT_T1Day9 Y_OUT
%     
% end
end

