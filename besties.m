function [chi2_best,n_best,Klp_best] = besties(mar,E,E_err,ksn,ksn_err,theta,n,nbins,cutoff,varargin)
%-------------------------------------------------------------------------%
% Function written by Byron A. Adams - Updated: 3 Aug 2019
%-------------------------------------------------------------------------%
%
% Description:
% Function finds the parital erosional efficiency component Klp, n, and m 
% values that best fit E vs ksn observed data, and the stream-power model.
% The script calculates the erosional efficiencies predicted from 
% regressions of channel steepness and erosion rate data. Before 
% regressions are carried out, the data are parsed into mean annual 
% rainfall bins. The rainfall bins are calculated based on R_max, R_min and 
% nbins, where nbins is the maximum number of bins that will be used to 
% divded the R_max - R_min range. The function then calculates further 
% regressions where the rainfall range is divided by nbins - 1, nbins - 2, 
% and so on until nbins = 1. This method creates a large number of semi-
% random bins of varying size and center value. This process is repeated 
% for a combintation of n and Klp values provided by the user. The chi 
% squared value of the predicted and observed ksn values is calculated. The 
% best-fit n andKlp values are found where the chi^2 is minimized. m is 
% assumed to be equal to theta*n.
%
% Usage:
% besties(mar,E,E_err,ksn,ksn_err,ksn_q,ksn_q_err,m,n,Klp,nbins,cutoff);
% [chi2_ksn,chi2_best,n_best,Klp_best] = besties(mar,E,E_err,ksn,ksn_err,ksn_q,ksn_q_err,m,n,Klp,nbins,cutoff,'num',1e6...);
%
% Required Inputs:
% mar - mean annual rainfall (m/yr)
% E - erosion rate (m/Myr)
% E_err - 1 sigma uncertainites on the erosion rate data (m/Myr)
% ksn - normalized channel steepness data
% ksn_err - 1 sigma uncertainites on the channel steepness data
% theta - reference concavity (unitless)(scalar)
% n - an array of slope exponents for stream-power law (unitless)
% nbins - maximum number or rainfall bins to created within the rainfall
%         range
% cutoff - this is the minimum number of data that will be included i
%          reported regression (scalar)
%
% Optional Inputs:
% num - number of monte carlo simulations used in regressions
% R_max - maximum rainfall value to test (m/yr)
% R_min - minimum rainfall value to test (m/yr)
%
% Outputs:
% chi2_best - minimum chi squared value
% n_best - the value of n that yields the minimum root mean square error
% Klp_best - the Klp value that yields the minimum root mean square error
%
%-------------------------------------------------------------------------%
% tashi delek!
%-------------------------------------------------------------------------%
%
% parse inputs
	p = inputParser;         
	p.FunctionName = 'find_Klpn_beta';
    
	addRequired(p,'mar',@(x) isnumeric(x));
	addRequired(p,'E',@(x) isnumeric(x));
    addRequired(p,'E_err',@(x) isnumeric(x));
    addRequired(p,'ksn',@(x) isnumeric(x));
    addRequired(p,'ksn_err',@(x) isnumeric(x));
    addRequired(p,'theta',@(x) isscalar(x));
    addRequired(p,'n',@(x) isnumeric(x));
    addRequired(p,'nbins',@(x) isscalar(x));
    
	addParameter(p,'num',1e5,@(x) isscalar(x));
	addParameter(p,'R_max',max(mar),@(x) isscalar(x));
	addParameter(p,'R_min',min(mar),@(x) isscalar(x));
    addParameter(p,'R_step',0.1,@(x) isscalar(x));
    
	parse(p,mar,E,E_err,ksn,ksn_err,theta,n,nbins,varargin{:});
	mar = p.Results.mar;
	E = p.Results.E;
    E_err = p.Results.E_err;
    ksn = p.Results.ksn;
    ksn_err = p.Results.ksn_err;
    theta = p.Results.theta;
    n = p.Results.n;
    nbins = p.Results.nbins;
    
	num = p.Results.num;
	R_max = p.Results.R_max;
	R_min = p.Results.R_min;
    R_step = p.Results.R_step;

% initialize variables
    MAR_data = [];
    num_bins = fliplr(1:1:nbins);
    R_plot = 0:R_step:R_max;

% loop through n values
	for q = 1:length(n)
        K_data = [];
        K_std_data = [];
        R_data = [];
        size = [];
        % loop through each bin center
        for j = 1:nbins
                clear divisions centers
            % find rainfall division points (bin edges)
                R_step = (R_max - R_min)/num_bins(j);
                divisions = R_min:R_step:R_max;
            
            % divide data according to division points
                group(j,:) = discretize(mar,divisions);
            
            % calculate bin centers based on the divisions
                for k = 1:(length(divisions) - 1)
                    centers(k) = (divisions(k + 1) + divisions(k))/2; %#ok<*AGROW>
                end
            
            % concatonate data
                MAR_data = [MAR_data centers];
            
                for i = 1:nbins
                    ind = find(group(j,:) == i);
                    E_new = E(ind);
                    E_err_new = E_err(ind);
                    ksn_new = ksn(ind);
                    ksn_err_new = ksn_err(ind);
                    size = [size num_bins(j)];
                    
                    % call MC_York to regress ksn - E data
                        [C,~,~,~,C_std] = MC_York(n(q),E_new,ksn_new,E_err_new,ksn_err_new,'num',num);
                        K_LSE = (C.^-n(q))/1e6;
                        K_err = C_std/C*K_LSE;
                        
                    % if there are enough data in each bin then record
                    % output
                        if length(E_new) >= cutoff
                            K_data = [K_data K_LSE];
                            K_std_data = [K_std_data K_err];
                            R_data = [R_data centers(i)];
                        elseif length(E_new) < cutoff
                            K_data = [K_data NaN];
                            K_std_data = [K_std_data NaN];
                            R_data = [R_data NaN];
                        end
                end
        end
            K_all(q,:) = K_data;
            K_std_all(q,:) = K_std_data;
            Klp(q) = nanmean(K_data./(R_data.^(theta.*n(q))));
           
            % find the root mean square estimator for the n and Kl combination
                K = (Klp(q).*1e6.*mar).^(theta*n(q));
                ksn_predicted = (E./K).^(1./n(q));
                chi2_ksn(q) = sum(((ksn - ksn_predicted).^2)./(ksn_err.^2));
	end

% find best-fit stats
    chi2_best = min(chi2_ksn);
    Klp_best = Klp(chi2_ksn == chi2_best);
    n_best = n(chi2_ksn == chi2_best);

% plot stuff
    figure(1)
    hold on
    yyaxis left
    ax1 = plot(chi2_ksn,n,'-k','LineWidth',2);
    ax1.Parent.YColor = [0 0 0];
    ticks_y = 0:0.1:max(n);
    set(gca,'YTick',ticks_y)
    plot(chi2_ksn,n*theta,'-g','LineWidth',2)
    set(gca,'fontsize',14,'FontName','Arial')
    xlabel('\chi^2','FontSize',18,'FontWeight','bold','FontName','Arial')
    ylabel('n (black); m (green)','FontSize',18,'FontWeight','bold','FontName','Arial','Color','k')

    yyaxis right
    ax2 = plot(chi2_ksn,Klp,'-b','LineWidth',2);
    ax2.Parent.YColor = [0 0 0];
    set(gca, 'YScale', 'log')
    ylabel('K_l_p (blue)','FontSize',18,'FontWeight','bold','FontName','Arial','Color','k')
    
% plot regression data
    index = find(n == n_best);
	figure(2)
	hold on
	yyaxis left
    e = errorbar(fliplr(R_data),fliplr(K_all(index,:))*1e8,fliplr(K_std_all(index,:))*1e8*2,'.k');
    e.LineWidth = 1;
    e.CapSize = 0;
    ax1 = scatter(fliplr(R_data),fliplr(K_all(index,:))*1e8,size*15,fliplr(size),'filled','MarkerEdgeColor',[0 0 0]);
    ax1.Parent.YColor = [0 0 0];
	colormap(by_map)
	set(gca,'fontsize',14,'FontName','Arial')
	xlabel('Mean annual rainfall (MAR, m yr^-^1)','FontSize',18,'FontWeight','bold','FontName','Arial')
	ylabel('Erosion coefficient (K, x10^-^8)','FontSize',18,'FontWeight','bold','FontName','Arial')
	ylim([0 4])
	xlim([0 6])
   
% plot stream-power model
    yyaxis left
    K_SP = Klp(index).*R_plot.^(theta*n(index));
    line(R_plot,K_SP*1e8,'LineStyle','-','Color','k','LineWidth',2.5)
   
% plot histogram data      
    hist_bins = round(R_max,1)/0.1;
    edges = (1:1:hist_bins)*0.1;
    yyaxis right
    h = histogram(mar,edges);
    h.Parent.YColor = [0 0 0];
    set(gca,'fontsize',14,'FontName','Arial')
    set(gca,'Ydir','reverse')
    h.FaceColor = [1 1 1];
    h.EdgeColor = [0 0 0];
    h.LineWidth = 1;
    h.Parent.XMinorTick = 'on';
    h.Parent.YMinorTick = 'on';
    h.Parent.LineWidth = 1;
    ylabel('MAR bin sample counts','FontSize',18,'FontWeight','bold','FontName','Arial')
    box on