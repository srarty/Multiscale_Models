%Models a neural mass column
%
% ************ INPUTS *************
% y - state vector of membrane voltages and parameters
% params - struct containing model parameters

%************* OUTPUTS ***************
% y_out - next state vector

function out = neural_mass_model(y, Q)


% Time axis
% tmin=params.tmin;
% tmax=params.tmax;  % s  % run time (seconds)
% dt=params.dt;      %  time step 
% n=params.n;
%t=tmin:dt:tmax;% vector for all time points
y0=y;


% options=odeset('RelTol',1e-3,'AbsTol',1e-4);
% [~,ysol]=ode45(@(t,Y) whole_column_matrix(Y,params),t,y0);
[Y_new,H]=whole_column_matrix_euler(y0,Q);
% [~,H]=whole_column_matrix(y0,params);
% output1=H*ysol';
out=Y_new;


end