clear

%% LIF
% folder = 'C:\Users\artemios\Documents\Multiscale_Models_Data\2023\excitability_spartan\no_current_pulse\'; 
folder = 'C:\Users\artemios\Documents\Multiscale_Models_Data\2023\excitability_spartan\'; 

var_vec = {'py_j_AMPA', 'in_j_AMPA', 'py_j_GABA_', 'py_j_GABAb_', 'in_j_GABA_'};
% var_vec = {'py_j_AMPA'};

for ii = 1:length(var_vec)
    u_vector = 0.1:0.1:2.0;
    var_ = var_vec{ii};
    for jj = 1:numel(u_vector)
        j = u_vector(jj);
        d = dir([folder '*' var_ '*_u' sprintf('%0.1f',j) '*']);

        if isempty(d),continue; end

        % Sort files
        u_value = [];
        for i = 1:length(d), u_value(i) = str2double(d(i).name( strfind(d(i).name,'_u')+2 : strfind(d(i).name,'_u') + 4) ); end
        [~, idx] = sort(u_value);
        d = d(idx);

        for i = 1:length(d)
            % Load file
            a = load([folder filesep d(i).name]);
            y_lif = a.LFP_V;
            t_lif = a.lfp_dt:a.lfp_dt:length(a.LFP_V)*a.lfp_dt;
            % Analyze excitabilty
            tr(i,jj) = analyze_excitability(y_lif', t_lif', 4899, -1.1, 2500, false);
            fr_in(i,jj) = mean(a.R_in(1000:4000));
            fr_py(i,jj) = mean(a.R_py(1000:4000));
            values_(i,jj) = a.value_;
        end
        
        [~, idx] = sort(values_(:,jj));
        tr(:,jj) = tr(idx,jj);
        fr_in(:,jj) = fr_in(idx,jj);
        fr_py(:,jj) = fr_py(idx,jj);
        values_(:,jj) = values_(idx,jj);
        
    end


    %% Plot LIF
    f_handle = figure;
    
    subplot(121)
    % mesh(values_(:,1), u_vector, fr_py', 'FaceColor', 'flat', 'EdgeColor', 'none')
    imagesc(values_(:,1), u_vector, fr_py');
    xlabel(var_);
    ylabel('u');
    zlabel('fr_{py}');
    hold
    title('LIF');
    colormap cool
    c = colorbar; % caxis(limits); c.Limits = limits; zlim(limits);
        
    subplot(122)
    imagesc(values_(:,1), u_vector, fr_in');
    % xlabel('\tau_{m_{i}}');
    % ylabel('\tau_{m_{e}}');
    xlabel(var_);%('Input rate');
    ylabel('u');
    zlabel('fr_{in}');
    hold
    % plot3(0.01638,0.008115,25.6339,'rx','LineWidth',3)
    %     title('NMM | u = 9');
    title('LIF');
    colormap cool
    c = colorbar;
%     caxis(limits);
%     c.Limits = limits;
%     zlim(limits);
%%    
    figure;
    mesh(values_(:,1), u_vector, tr', 'FaceColor', 'flat', 'EdgeColor', 'none')
    % xlabel('\tau_{m_{i}}');
    % ylabel('\tau_{m_{e}}');
    xlabel(var_);%('Input rate');
    ylabel('u');
    zlabel('t_{recovery}');
    hold
    % plot3(0.01638,0.008115,25.6339,'rx','LineWidth',3)
    %     title('NMM | u = 9');
    title('LIF');
    colormap cool
    c = colorbar;
%     c.Label.String = z_label;
%     caxis(limits);
%     c.Limits = limits;
%     zlim(limits);

end

%% NMM
% TODO