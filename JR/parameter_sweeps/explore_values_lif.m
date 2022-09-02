%% Goes through files generated by the LIF model in Spartan.
% Files saved in "...\Spartan\sweep\" are simulated by varying u (external
% input) and alpha_xy (the gain, called 'j' in the LIF). _xy refers to the
% synapses, where x is the pre-synaptic population and y is the
% post-synaptic population.
%
% Artemio - August 2022
% function explore_values_lif(TYPE)
TYPE = 'ii'; % type = {'ee', 'ii', 'ie', 'ei'} % xy, x=pre-synaptic, y=post-synaptic

root = 'C:\Users\artemios\Documents\Multiscale_Models_Data\Spartan\sweep\';
d = dir([root '*_' TYPE '*.mat']);

% Declare empty arrays
u = nan(size(d));
alpha = nan(size(d));
for i = 1:length(d)
    % Retrieve current name
    n = d(i).name;
    
    % Find the positino of u and alpha in the file name
    u_idx = strfind(n,'_u') + 2;
    alpha_idx = strfind(n, strcat('_', TYPE)) + length(TYPE) + 1;
    end_idx = strfind(n, '.mat') - 1;
    
    % Retrieve u and alpha
    u(i) = str2double(n(u_idx : alpha_idx -(length(TYPE) + 2)));
    alpha(i) = str2double(n(alpha_idx : end_idx));
    
end

[u_unique,~,u_order] = unique(round(u,3));
[alpha_unique,~,alpha_order] = unique(round(alpha,3));
f_py = nan(length(u_unique), length(alpha_unique));
f_in = nan(length(u_unique), length(alpha_unique));
for i = 1:length(d)
    % Load the file
    data = load([root d(i).name]);

    f_py(u_order(i), alpha_order(i)) = mean(data.R_py(2000:end-1500));
    f_in(u_order(i), alpha_order(i)) = mean(data.R_in(2000:end-1500));
end
%% Plot spike rates mesh
% Set label and View
if TYPE(1)=='e', angle=([0 -90]); else, angle=([0 90]); end
    
figure('Position', [325 404 1112 364]);
ax = subplot(1,2,1);

mesh(u_unique, -1*alpha_unique, f_py', 'FaceColor', 'flat', 'EdgeColor', 'none')
xlabel('u');%('Input rate');
ylabel(['j_{' TYPE '}']);
zlabel('Firing rate (Py)');
c = colorbar;
c.Label.String = 'Mean firing rate (Hz)';
caxis([0 1.5]);
c.Limits = [0 1.5];
xlim([0 5]);
% zlim([0 1.5]);
% title('Pyramidal')
ax.View = (angle);
%%
ax = subplot(1,2,2);
mesh(u_unique, -1*alpha_unique, f_in', 'FaceColor', 'flat', 'EdgeColor', 'none')
xlabel('u');%('Input rate');
ylabel(['j_{' TYPE '}']);
zlabel('Firing rate (In)');
% title('Inhibitory');
colormap jet
c = colorbar;
c.Label.String = 'Mean firing rate (Hz)';
caxis([0 3.5]);
c.Limits = [0 3.5];
xlim([0 5])
ax.View = (angle);