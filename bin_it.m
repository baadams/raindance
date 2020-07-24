function [t] = bin_it(mar,E,E_err,ksn,ksn_err,ksn_q,ksn_q_err,m,n,Klp,divisions,varargin)
%-------------------------------------------------------------------------%
% Function written by Byron A. Adams - Updated: 11 Apr 2019
%-------------------------------------------------------------------------%
%
% Description:
% This script was written to explore the sensitivity of fluvial relief and
% erosion rates to rainfall rates. It relies on a priory information about
% probable values for stream-power model parameters - specifially, m, n and
% Klp. It requires commonly measured metrics such as MAR, ksn, E, but also
% ksn-q (a metrics that uses dischare instead of area to calculate channel
% steepness). It is useful to look at data bined by mean annual rainfall,
% and therefore, the user is also asked for guestimates of such bins. 
%
% Usage:
% bin_it(mar,E,E_err,ksn,ksn_err,ksn_q,ksn_q_err,m,n,Klp,divisions);
% [table_1] = bin_it(mar,E,E_err,ksn,ksn_err,ksn_q,ksn_q_err,m,n,Klp,divisions,'name',value,...);
%
% Required Inputs:
% mar - mean annual rainfall (m/yr)
% E - erosion rate (m/Myr)
% E_err - 1 sigma uncertainites on the erosion rate data (m/Myr)
% ksn - normalized channel steepness data
% ksn_err - 1 sigma uncertainites on the channel steepness data
% ksn-q - discharge-based channel steepness data 
% ksn-q_err - 1 sigma uncertainites on the channel steepness data
% m - water discharge exponent for stream-power law (unitless) (scalar)
% n - slope exponent for stream-power law (unitless) (scalar)
% Klp - partial coefficient of erosion (L^1-3m T^m-1). Note that L is
%       in meters and T is in years. (scalar)
% divisions - precipitation division points (bin edges)(m/yr) (array). 
%             E.g. [0 1.5 2.5 3.5 6]
%
% Optional Inputs:
% num - number of monte carlo simulations for regressions
% E_max - erosion limit for plots (m/Myr)
% ksn_max - ksn limit for plots
% n_SP - slope exponent to be used in the stream-power law. default is 2.
% print - set to 'y' to save data tables and plots
%
% Outputs:
% E - ksn and E - ksn-q plots
% tables of E - ksn and E - ksn-q relationship parameters
%
%-------------------------------------------------------------------------%
% tashi delek!
%-------------------------------------------------------------------------%

% parse inputs
	p = inputParser;         
	p.FunctionName = 'bin_it';
    
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
    
	addParameter(p,'num',1e6,@(x) isscalar(x));
	addParameter(p,'E_max',max(E),@(x) isscalar(x));
	addParameter(p,'ksn_max',max(ksn),@(x) isscalar(x));
    addParameter(p,'n_SP',2,@(x) isscalar(x));
    addParameter(p,'print','n',@(x) ischar(x));
    
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
    
	num = p.Results.num;
	E_max = p.Results.E_max;
	ksn_max = p.Results.ksn_max;
    n_SP = p.Results.n_SP;
    print = p.Results.print;

% calcualte the size of each rainfall bin
    for s = 1:length(divisions)-1
        bin_size(s) = divisions(s + 1) - divisions(s);
    end

% divide data according to divisions
    group = discretize(mar,divisions);

% calculate bin centers based on divisions
    for k = 1:(length(divisions) - 1)
        centers(k) = (divisions(k + 1) + divisions(k))/2; %#ok<*AGROW>
    end

% assign each data point a bin based on the bin edges above
    for i = 1:length(group)
        mar_group(i) = centers(group(i));
    end

% intialize matricies
    C_LSE_fix = NaN(length(centers)*2+2,1);
    K_LSE_fix = NaN(length(centers)*2+2,1);
    phi_LSE_fix = NaN(length(centers)*2+2,1);
    n_LSE_fix = NaN(length(centers)*2+2,1);
    MSWD_LSE_fix = NaN(length(centers)*2+2,1);
    MSWD_std_LSE_fix = NaN(length(centers)*2+2,1);
    C_LSE_free = NaN(length(centers)*2+2,1);
    K_LSE_free = NaN(length(centers)*2+2,1);
    n_LSE_free = NaN(length(centers)*2+2,1);
    phi_LSE_free = NaN(length(centers)*2+2,1);
    MSWD_LSE_free = NaN(length(centers)*2+2,1);
    MSWD_std_LSE_free = NaN(length(centers)*2+2,1);
    MSWD_std_SP = NaN(length(centers)*2+2,1);
    MSWD_SP = NaN(length(centers)*2+2,1);
    K_SP = NaN(length(centers)*2+2,1);
    name = cell(length(centers)*2+2,1);

% calculate the complete coefficient of erosion of the stream-power 
% model in units of m and Myr
    K = (1./((Klp.*1e6.*centers).^m)).^(1/n_SP);

% plot a histogram to see the distribution of data into rainfall bins
    figure(1)
    histogram(mar,divisions);
    ylabel('Counts','FontSize',14,'FontName','Arial','FontWeight','bold')
    xlabel('Mean annual rainfall (m yr^-^1)','FontSize',14,'FontName','Arial','FontWeight','bold')
    xlim([0 max(divisions)])
    set(gca,'fontsize',12,'FontName','Arial')

% loop through each rainfall bin to divide samples and regress
    for j = 1:max(group)
        ind = find(group == j);
        name{j} = ['ksn_' num2str(centers(j))];
        E_new = E(ind);
        E_err_new = E_err(ind);
        ksn_new = ksn(ind);
        ksn_err_new = ksn_err(ind);
        mar_grp_new = mar_group(ind);
        mar_new = mar(ind);
        
        % plot nested histograms to show the distrubution in each bin
            figure(1)
            hold on
            num_bins = bin_size(j)/0.1;
            edges = (1:1:num_bins)*0.1 + divisions(j);
            histogram(mar_new,edges)
        
        % plot regressions
            figure(2)
            subplot(2,2,j);
            hold on
        
        % calculate the best-fit power-law for the rainfall bin with n free
        	[C_LSE_free(j),n_LSE_free(j),MSWD_LSE_free(j),MSWD_std_LSE_free(j),~,x_plot,y_plot] = MC_York(0,E_new,ksn_new,E_err_new,ksn_err_new,'num',num,'x_max',E_max);
            plot(x_plot,y_plot,'LineStyle','-','Color',[0.5 0.5 0.5],'LineWidth',2)
            K_LSE_free(j) = (C_LSE_free(j))^-n_LSE_free(j)*1e-6;
            phi_LSE_free(j) = 1/n_LSE_free(j);
            
        % calculate the best-fit power-law for the rainfall bin with n fixed
        	[C_LSE_fix(j),n_LSE_fix(j),MSWD_LSE_fix(j),MSWD_std_LSE_fix(j),~,x_plot,y_plot] = MC_York(n,E_new,ksn_new,E_err_new,ksn_err_new,'num',num,'x_max',E_max);
            plot(x_plot,y_plot,'LineStyle','-','Color',[0 0 0],'LineWidth',2)
            K_LSE_fix(j) = (C_LSE_fix(j))^-n_LSE_fix(j)*1e-6;
            phi_LSE_fix(j) = 1/n_LSE_fix(j);
        
        % calculate the stream-power model for the rainfall bin
         	[x_plot,y_plot] = SPM(m,n_SP,Klp,mar_grp_new(1),E_max);
            plot(x_plot,y_plot,'LineStyle','--','Color',[0 0 0],'LineWidth',2)
        
        % calculate the MSWD of the stream-power law
        	chi_sq = sum(((E_new - (1/K(j)^n_SP).*(ksn_new.^n_SP)).^2)./((E_err_new.^2) + ((1/K(j)^n_SP).*n_SP.*ksn_new.^(n_SP - 1)).^2.*(ksn_err_new.^2)));
        	MSWD_SP(j) = chi_sq/(length(E_new) - 2);
        	MSWD_std_SP(j) = sqrt(2/(length(E_new) - 2));
            K_SP(j) = Klp*(centers(j)^m);
        
        % plot original sample data and add graphical details
        	figure(2)
            subplot(2,2,j)
        	hold on
        	e = errorbar(E_new,ksn_new,ksn_err_new*2,ksn_err_new*2,E_err_new*2,E_err_new*2,'.k');
            e.LineWidth = 1;
            e.CapSize = 0;
        	scatter(E_new,ksn_new,30,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0 0 0])
            xlabel('Erosion rate (m Myr^-^1)','FontSize',14,'FontWeight','bold','FontName','Arial')
            xlim([0 E_max])
            ylabel('Channel steepness (k_s_n)','FontSize',14,'FontName','Arial','FontWeight','bold')
            ylim([0 ksn_max]);
            legend({'free exponent' 'fixed exponent' 'stream-power model'},'FontSize',10,'FontName','Arial','Location','southeast');
            title(['MAR = ' num2str(centers(j)) ' ' '(N = ' num2str(length(E_new)) ')'],'FontSize',14,'FontWeight','bold','FontName','Arial')
            set(gca,'fontsize',12,'FontName','Arial')
            box on
    end

% plot all original data and stream-power curves
	for l = 1:max(group)
        ind = group == l;
        mar_grp_new = mar_group(ind);
        colours = colormap(mar_map(max(group)));
        
        % plot group data
            figure(3)
            subplot(2,1,1)
            hold on
        
        % plot stream-power curves
            [x_plot,y_plot] = SPM(m,n_SP,Klp,mar_grp_new(1),E_max);
            plot(x_plot,y_plot,'LineStyle','--','Color',colours(l,:),'LineWidth',2,'DisplayName',['SPM' ' ' num2str(mar_grp_new(1))])
	end
    
    % use Klp to predict ksn for range of E. this is equivalent to using
    % the stream-power law
    	E_model = 0:10:E_max;
        ksn_model = (E_model./(Klp*1e6*median(mar).^m)).^(1/n_SP);
        figure(3)
        subplot(2,1,1)
        hold on
        line(E_model,ksn_model,'Color',[0 0 0],'LineStyle','--','LineWidth',2,'DisplayName','SPM (all)')
    
    % calculate the MSWD of the stream-power model
    	chi_sq = sum(((E - ((Klp*1e6*median(mar)).^m).*(ksn.^n_SP)).^2)./((E_err.^2) + (((Klp*1e6*median(mar)).^m).*n_SP.*ksn.^(n_SP - 1)).^2.*(ksn_err.^2)));
    	MSWD_SP(j+1) = chi_sq/(length(E) - 2);
        MSWD_std_SP(j+1) = sqrt(2/(length(E) - 2));
        K_SP(j+1) = Klp*(median(mar)^m);
   
    % calculate the best-fit relationship for all data
        name{j+1} = 'ksn_all';
        [C_LSE_fix(j+1),n_LSE_fix(j+1),MSWD_LSE_fix(j+1),MSWD_std_LSE_fix(j+1),~,x_plot,y_plot] = MC_York(n,E,ksn,E_err,ksn_err,'num',num,'x_max',E_max);
        line(x_plot,y_plot,'Color',[0 0 0],'LineWidth',2,'DisplayName','LSE fix (all)')
        K_LSE_fix(j+1) = (C_LSE_fix(j+1))^-n_LSE_fix(j+1)*1e-6;
        phi_LSE_fix(j+1) = 1/n_LSE_fix(j+1);
        
        [C_LSE_free(j+1),n_LSE_free(j+1),MSWD_LSE_free(j+1),MSWD_std_LSE_free(j+1),~,~,~] = MC_York(0,E,ksn,E_err,ksn_err,'num',num,'x_max',E_max);
        K_LSE_free(j+1) = (C_LSE_free(j+1))^-n_LSE_free(j+1)*1e-6;
        phi_LSE_free(j+1) = 1/n_LSE_free(j+1);
        
    % make the plot pretty
        scatter(E,ksn,30,mar,'filled','MarkerEdgeColor',[0 0 0]);
        colormap(mar_map)
        set(gca,'fontsize',14,'FontName','Arial')
        xlabel('Erosion rate (m Myr^-^1)','FontSize',18,'FontWeight','bold','FontName','Arial')
        xlim([0 E_max])
        ylabel('Channel steepness (k_s_n)','FontSize',18,'FontWeight','bold','FontName','Arial')
        ylim([0 ksn_max])
        box on
        legend('show','Location','southeast','NumColumns',3,'FontSize',8)

% plot all ksn-q data, regressions of each bin, and a regression and
% stream-power curve for all data
	figure(3)
    subplot(2,1,2)
	hold on
    
    % use Klp to predict ksn-q for range of E. this is equivalent to using
    % the stream-power law
    	ksn_q_model = (E_model./(Klp*1e6)).^(1/n_SP);
        line(E_model,ksn_q_model,'Color',[0 0 0],'LineWidth',2,'LineStyle','--','DisplayName','SPM (all)')
        set(gca,'fontsize',14,'FontName','Arial')
        xlabel('Erosion rate (m Myr^-^1)','FontSize',18,'FontWeight','bold','FontName','Arial')
        xlim([0 E_max])
        ylabel('Channel steepness (k_s_n-q)','FontSize',18,'FontWeight','bold','FontName','Arial')
        ylim([0 600])
        box on
        legend('show','Location','southeast','NumColumns',3,'FontSize',8)
    
    % calculate the MSWD of the stream-power law
    	chi_sq = sum(((E - (Klp*1e6).*(ksn_q.^n_SP)).^2)./((E_err.^2) + (Klp*1e6.*n_SP.*ksn_q.^(n_SP - 1)).^2.*(ksn_q_err.^2)));
    	MSWD_SP(j+2) = chi_sq/(length(E) - 2);
    	MSWD_std_SP(j+2) = sqrt(2/(length(E) - 2));
	
    % calculate the best-fit relationship for all data
    	name{j+2} = 'ksn_q_all';
        [C_LSE_fix(j+2),n_LSE_fix(j+2),MSWD_LSE_fix(j+2),MSWD_std_LSE_fix(j+2),~,x_plot,y_plot] = MC_York(n,E,ksn_q,E_err,ksn_q_err,'num',num,'x_max',E_max);
        line(x_plot,y_plot,'Color',[0 0 0],'LineWidth',2,'DisplayName','LSE fix (all)')
        K_LSE_fix(j+2) = (C_LSE_fix(j+2))^-n_LSE_fix(j+2)*1e-6;
        phi_LSE_fix(j+2) = 1/n_LSE_fix(j+2);
        
        [C_LSE_free(j+2),n_LSE_free(j+2),MSWD_LSE_free(j+2),MSWD_std_LSE_free(j+2),~,~,~] = MC_York(0,E,ksn_q,E_err,ksn_q_err,'num',num,'x_max',E_max);
        K_LSE_free(j+2) = (C_LSE_free(j+2))^-n_LSE_free(j+2)*1e-6;
        phi_LSE_free(j+2) = 1/n_LSE_free(j+2);
    
    % split ksn-q data according to bins
        for l = 1:max(group)
            ind = find(group == l);
            E_new = E(ind);
            E_err_new = E_err(ind);
            ksn_q_new = ksn_q(ind);
            ksn_q_err_new = ksn_q_err(ind);
            mar_grp_new = mar_group(ind);
            
            % plot group data
                scatter(E_new,ksn_q_new,30,'MarkerFaceColor',colours(l,:),'MarkerEdgeColor',[0 0 0],'DisplayName',['data' ' ' num2str(mar_grp_new(1))]);
            
            % plot regression curves
                name{l+j+2} = ['ksn_q_' num2str(centers(l))];
                [C_LSE_fix(j+l+2),n_LSE_fix(j+l+2),MSWD_LSE_fix(j+l+2),MSWD_std_LSE_fix(j+l+2),~,x_plot,y_plot] = MC_York(n,E_new,ksn_q_new,E_err_new,ksn_q_err_new,'num',num,'x_max',E_max);
                plot(x_plot,y_plot,'LineStyle','-','Color',colours(l,:),'LineWidth',2,'DisplayName',['LSE fix' ' ' num2str(mar_grp_new(1))])
                K_LSE_fix(l+j+2) = (C_LSE_fix(l+j+2))^-n_LSE_fix(l+j+2)*1e-6;
                phi_LSE_fix(l+j+2) = 1/n_LSE_fix(l+j+2);
                
                [C_LSE_free(j+l+2),n_LSE_free(j+l+2),MSWD_LSE_free(j+l+2),MSWD_std_LSE_free(j+l+2),~,~,~] = MC_York(0,E_new,ksn_q_new,E_err_new,ksn_q_err_new,'num',num,'x_max',E_max);
                K_LSE_free(l+j+2) = (C_LSE_free(l+j+2))^-n_LSE_free(l+j+2)*1e-6;
                phi_LSE_free(l+j+2) = 1/n_LSE_free(l+j+2);
                
           % calculate the MSWD of the stream-power law
                chi_sq = sum(((E_new - (Klp*1e6).*(ksn_q_new.^n_SP)).^2)./((E_err_new.^2) + (Klp*1e6.*n_SP.*ksn_q_new.^(n_SP - 1)).^2.*(ksn_q_err_new.^2)));
                MSWD_SP(l+j+2) = chi_sq/(length(E_new) - 2);
                MSWD_std_SP(l+j+2) = sqrt(2/(length(E_new) - 2));
        end
        
% collate data in tables
    t = table(name,C_LSE_fix,phi_LSE_fix,K_LSE_fix,n_LSE_fix,MSWD_LSE_fix,MSWD_std_LSE_fix,C_LSE_free,phi_LSE_free,K_LSE_free,n_LSE_free,MSWD_LSE_free,MSWD_std_LSE_free,K_SP,MSWD_SP,MSWD_std_SP);

% save tables and figures
    if print == 'y'
        save 'summary_table.mat' 't'
        writetable(t,'summary_table.csv','Delimiter',',')
        
        for i = 1:3
            saveas(i,['fig_' num2str(i)],'epsc')
            saveas(i,['fig_' num2str(i)],'png')
        end
    end
