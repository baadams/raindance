function [table_1,table_2] = bin_it(mar,E,E_err,ksn,ksn_err,ksn_q,ksn_q_err,m,n,Klp,divisions,varargin)
%-------------------------------------------------------------------------%
% Function written by Byron A. Adams - Updated: 11 Apr 2019
%-------------------------------------------------------------------------%
%
% Description:
% Function finds the best-fit power law between E and ksn within climate 
% bins. the stream-power law is also used to predict the ksn-E relationship 
% given the mean annual precipitation bin value.
%
% Usage:
% [table_1,table_2] = bin_it(mar,E,E_err,ksn,ksn_err,ksn_q,ksn_q_err,m,n,Klp,divisions);
% [table_1,table_2] = bin_it(mar,E,E_err,ksn,ksn_err,ksn_q,ksn_q_err,m,n,Klp,divisions,'name',value,...);
%
% Required Inputs:
% mar - mean annual rainfall/precipitation (m/yr)
% E - erosion rate (m/Myr)
% E_err - 1 sigma uncertainites on the erosion rate data (m/Myr)
% ksn - normalized channel steepness data
% ksn_err - 1 sigma uncertainites on the channel steepness data
% ksn-q - 
% ksn-q_err - 
% m - water discharge exponent for stream-power law (unitless)
% n - slope exponent for stream-power law (unitless)
% Klp - partial coefficient of erosion (L^1-3m T^m-1). Note that L is
%       in meters and T is in years
% divisions - precipitation division points (bin edges)(m/yr)
%
% Optional Inputs:
% num - number of monte carlo simulations 
% E_max - erosion limit for plots (m/Myr)
% ksn_max - ksn limit for plots
% n_SP - slope exponent to be used in the stream-power law
% print - set to 'y' to save data tables and plots
%
% Outputs:
% E - ksn and E - ksn-q plots
% tables of E - ksn and E - ksn-q relationship parameters
%
%-------------------------------------------------------------------------%
% tashi delek
%-------------------------------------------------------------------------%
%
% parse inputs
	p = inputParser;         
	p.FunctionName = 'bin_it';
    %
	addRequired(p,'mar',@(x) isnumeric(x));
	addRequired(p,'E',@(x) isnumeric(x));
    addRequired(p,'E_err',@(x) isnumeric(x));
    addRequired(p,'ksn',@(x) isnumeric(x));
    addRequired(p,'ksn_err',@(x) isnumeric(x));
    addRequired(p,'ksn_q',@(x) isnumeric(x));
    addRequired(p,'ksn_q_err',@(x) isnumeric(x));
    addRequired(p,'m',@(x) isscalar(x));
    addRequired(p,'n',@(x) isnumeric(x));
    addRequired(p,'Klp',@(x) isscalar(x));
    addRequired(p,'divisions',@(x) isnumeric(x));
    %
	addParameter(p,'num',1e6,@(x) isscalar(x));
	addParameter(p,'E_max',max(E),@(x) isscalar(x));
	addParameter(p,'ksn_max',max(ksn),@(x) isscalar(x));
    addParameter(p,'n_SP',2.22,@(x) isscalar(x));
    addParameter(p,'print','n',@(x) ischar(x));
    %
	parse(p,mar,E,E_err,ksn,ksn_err,ksn_q,ksn_q_err,m,n,Klp,divisions,varargin{:});
	mar = p.Results.mar;
	E = p.Results.E;
    E_err = p.Results.E_err;
    ksn = p.Results.ksn;
    ksn_err = p.Results.ksn_err;
    ksn_q = p.Results.ksn_q;
    ksn_q_err = p.Results.ksn_q_err;
    m = p.Results.m;
    n = p.Results.n;
    Klp = p.Results.Klp;
    divisions = p.Results.divisions;
    %
	num = p.Results.num;
	E_max = p.Results.E_max;
	ksn_max = p.Results.ksn_max;
    n_SP = p.Results.n_SP;
    print = p.Results.print;
%
% calcualte the size of each MAR bin
    for s = 1:length(divisions)-1
        bin_size(s) = divisions(s + 1) - divisions(s);
    end
%
% divide data according to divsions
    group = discretize(mar,divisions);
%
% calculate bin centers based on divisions
    for k = 1:(length(divisions) - 1)
        centers(k) = (divisions(k + 1) + divisions(k))/2; %#ok<*AGROW>
    end
%
% assign each data point a bin based on the bin edges above
    for i = 1:length(group)
        mar_group(i) = centers(group(i));
    end
%
% intialize matricies
    C_LSE = zeros(length(centers),1);
    n_LSE = zeros(length(centers),1);
    MSWD_LSE = zeros(length(centers),1);
    MSWD_std_LSE = zeros(length(centers),1);
    MSWD_std_SP = zeros(length(centers),1);
    MSWD_SP = zeros(length(centers),1);
    name = cell(length(centers),1);
    C_LSE_q = zeros(length(centers),1);
    n_LSE_q = zeros(length(centers),1);
    MSWD_LSE_q = zeros(length(centers),1);
    MSWD_std_LSE_q = zeros(length(centers),1);
    name_q = cell(length(centers),1);
%
% calculate the complete coefficient of erosion values of the stream-power 
% model in units of m and Myr
    Klpr = (1./(Klp.*1e6.*centers)).^(1/n_SP);
%
% plot a histogram to see the distribution of data into MAR bins
    figure(1)
    histogram(mar,divisions);
    ylabel('Counts')
    xlabel('Mean annual rainfall (m/yr)')
    xlim([0 max(divisions)])
%
% loop through each MAR bin to divide samples and calculate regressions.
% create a plot of each bin regression and stream-power model
    for j = 1:max(group)
        ind = find(group == j);
        name{j} = ['ksn_' num2str(centers(j))];
        %
        E_new = E(ind);
        E_err_new = E_err(ind);
        ksn_new = ksn(ind);
        ksn_err_new = ksn_err(ind);
        mar_grp_new = mar_group(ind);
        mar_new = mar(ind);
        %
        % plot nested histograms to show the distrubution in each bin
            figure(1)
            hold on
            num_bins = bin_size(j)/0.1;
            edges = (1:1:num_bins)*0.1 + divisions(j);
            histogram(mar_new,edges)
        %
        % calculate the best-fit power-law realtionship for the MAR bin
        	[C_LSE(j),n_LSE(j),MSWD_LSE(j),MSWD_std_LSE(j)] = MC_York(1./n,E_new,ksn_new,E_err_new,ksn_err_new,num,'x_max',E_max,'fig_num',j+1,'colour',[0.5 0.5 0.5]);
        %
        % calculate the stream-power model for the MAR bin
         	stream_power(m,n_SP,Klp,mar_grp_new(1),E_max,j+1,[0 0 0]);
        %
        % calculate the MSWD of the stream-power law
        	chi_sq = sum(((E_new - (1/Klpr(j)^n_SP).*(ksn_new.^n_SP)).^2)./((E_err_new.^2) + ((1/Klpr(j)^n_SP).*n_SP.*ksn_new.^(n_SP - 1)).^2.*(ksn_err_new.^2)));
        	MSWD_SP(j) = chi_sq/(length(E_new) - 2);
        	MSWD_std_SP(j) = sqrt(2/(length(E_new) - 2));
        %
        % plot original sample data and add graph details
        	figure(j+1)
        	hold on
        	errorbar(E_new,ksn_new,ksn_err_new*2,ksn_err_new*2,E_err_new*2,E_err_new*2,'.k')
        	scatter(E_new,ksn_new,30,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0 0 0])
            %
            xlabel('Erosion rate (m/Myr)')
            xlim([0 E_max])
            ylabel('k_s_n')
            ylim([0 ksn_max]);
            legend({'LSE' 'SPM'});
            title(['MAR = ' num2str(centers(j)) ' ' '(N = ' num2str(length(E_new)) ')'])
    end
%
% plot all original data, regressions, and stream-power curves
	figure(j+2)
	hold on
	h3 = gscatter(E,ksn,mar_group,[],'.',20);
    xlabel('Erosion rate (m/Myr)')
	xlim([0 E_max])
    ylabel('k_s_n')
	ylim([0 ksn_max])
    %
	for l = 1:max(group)
        cgroup = h3(l);
        ind = find(group == l);
        E_new = E(ind);
        E_err_new = E_err(ind);
        ksn_new = ksn(ind);
        ksn_err_new = ksn_err(ind);
        mar_grp_new = mar_group(ind);
        %
        % plot regression curves
            [~,~,~,~] = MC_York(1./n,E_new,ksn_new,E_err_new,ksn_err_new,num,'x_max',E_max,'fig_num',j+2,'colour',cgroup.Color);
        %
        % plot stream-power curves
            stream_power(m,n_SP,Klp,mar_grp_new(1),E_max,j+2,cgroup.Color);
	end
%
% plot all ksn-q data, regressions of each bin, and a regression and
% stream-poer curve for all data
	figure(j+3)
	hold on
    % use Klp to predict ksn-q for range of E. this is equivalent to using
    % the stream-power law
    	E_model = 0:10:E_max;
    	ksn_q_model = (E_model./(Klp*1e6)).^(1/n_SP);
    %
    % calculate the MSWD of the stream-power law
    	chi_sq = sum(((E - (Klp*1e6).*(ksn_q.^n_SP)).^2)./((E_err.^2) + (Klp*1e6.*n_SP.*ksn_q.^(n_SP - 1)).^2.*(ksn_q_err.^2)));
    	MSWD_SP(j+1) = chi_sq/(length(E) - 2);
    	MSWD_std_SP(j+1) = sqrt(2/(length(E) - 2));
    %
    % plot the stream-power model
        h3 = gscatter(E,ksn_q,mar_group,[],'.',20);
        line(E_model,ksn_q_model,'Color',[0 0 0],'LineWidth',1,'DisplayName','SPM (all)')
        xlabel('Erosion rate (m/Myr)')
        xlim([0 E_max])
        ylabel('k_s_n- q')
        ylim([0 max(ksn_q)])
	%
    % calculate the best-fit relationship for all data
    	name{j+1} = 'ksn_q_all';
        [C_LSE(j+1),n_LSE(j+1),MSWD_LSE(j+1),MSWD_std_LSE(j+1)] = MC_York(1./n,E,ksn_q,E_err,ksn_q_err,num,'x_max',E_max,'fig_num',j+3,'colour',[0.5 0.5 0.5]);
    % 
    % split ksn-q data according to bins
        for l = 1:max(group)
            cgroup = h3(l);
            ind = find(group == l);
            E_new = E(ind);
            E_err_new = E_err(ind);
            ksn_q_new = ksn_q(ind);
            ksn_q_err_new = ksn_q_err(ind);
            %
            % plot regression curves
                name_q{l} = ['ksn_q_' num2str(centers(l))];
                [C_LSE_q(l),n_LSE_q(l),MSWD_LSE_q(l),MSWD_std_LSE_q(l)] = MC_York(1./n,E_new,ksn_q_new,E_err_new,ksn_q_err_new,num,'x_max',E_max,'fig_num',j+3,'colour',cgroup.Color);
        end
%
% plot original data and a regression and stream-power curve for all data
	figure(j+4)
	hold on
    % use Klp to predict ksn for range of E. this is equivalent to using
    % the stream-power law
    	ksn_model = (E_model./(Klp*1e6*median(mar))).^(1/n_SP);
    %
    % calculate the MSWD of the stream-power model
    	chi_sq = sum(((E - (Klp*1e6*median(mar)).*(ksn.^n_SP)).^2)./((E_err.^2) + ((Klp*1e6*median(mar)).*n_SP.*ksn.^(n_SP - 1)).^2.*(ksn_err.^2)));
    	MSWD_SP(j+2) = chi_sq/(length(E) - 2);
        MSWD_std_SP(j+2) = sqrt(2/(length(E) - 2));
    %
    % plot original data
        errorbar(E,ksn,ksn_err,ksn_err,E_err,E_err,'.k')
        scatter(E,ksn,30,mar,'filled');
        line(E_model,ksn_model,'Color',[0 0 0],'LineWidth',1,'DisplayName','SPM (all)')
    %
    % calculate the best-fit relationship for all data
        name{j+2} = 'ksn_all';
        [C_LSE(j+2),n_LSE(j+2),MSWD_LSE(j+2),MSWD_std_LSE(j+2)] = MC_York(1./n,E,ksn,E_err,ksn_err,num,'x_max',E_max,'fig_num',j+4,'colour',[0.5 0.5 0.5]);
    %
    % plot details
        xlabel('Erosion rate (m/Myr)')
        xlim([0 E_max])
        ylabel('k_s_n')
        ylim([0 ksn_max])
        caxis([0 max(mar)])
        colormap(mar_map)
        h = colorbar;
        ylabel(h,'Mean annual rainfall (m/yr)')
        legend('show')
%
% collate data in tables
    table_1 = table(name,C_LSE,n_LSE,MSWD_LSE,MSWD_std_LSE,MSWD_SP,MSWD_std_SP);
    table_2 = table(name_q,C_LSE_q,n_LSE_q,MSWD_LSE_q,MSWD_std_LSE_q);
%
% save tables and figures
    if print == 'y'
        save 'table_1.mat' 'table_1'
        save 'table_2.mat' 'table_2'
        %
        for i = 1:j+4
            saveas(i,['fig_' num2str(i)],'epsc')
            saveas(i,['fig_' num2str(i)],'png')
        end
    end
%