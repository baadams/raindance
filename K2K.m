function K2K(mar,E,E_err,ksn,ksn_err,m,n,Klp,nbins,varargin)
%-------------------------------------------------------------------------%
% Function written by Byron A. Adams - Updated: 3 Aug 2019
%-------------------------------------------------------------------------%
%
% Description:
% Function compares the erosional efficiencies predicted from regressions 
% of channel steepness and erosion rate data and erosional efficiencies
% predcited from the stream-power model for a given n and Klp pair. Before 
% regressions are carried out, the data are parsed based on the creation of 
% mean annual rainfall bins. The rainfall bins are calculated based on 
% R_max, R_min and nbins, where nbins is the maximum number of bins that 
% will be used to divded the R_max - R_min range. The function then 
% calculates further regressions where the rainfall range is divided by 
% nbins - 1, nbins - 2, and so on  until nbins = 1. This method creates a 
% large number of semi-random bins of varying size and center value.
%
% Usage:
% K2K(mar,E,E_err,ksn,ksn_err,ksn_q,ksn_q_err,m,n,Klp,nbins);
% K2K(mar,E,E_err,ksn,ksn_err,ksn_q,ksn_q_err,m,n,Klp,nbins,'num',1e6...);
%
% Required Inputs:
% mar - mean annual rainfall (m/yr)
% E - erosion rate (m/Myr)
% E_err - 1 sigma uncertainites on the erosion rate data (m/Myr)
% ksn - normalized channel steepness data
% ksn_err - 1 sigma uncertainites on the channel steepness data
% m - discharge exponent for stream-power law (unitless) (scalar)
% n - slope exponent for stream-power law (unitless) (scalar)
% Klp - partial coefficient of erosion (L^1-3m T^m-1) . Note that L is
%       (scalar) in meters and T is in years, and m is the discharge 
%       exponent
% nbins - maximum number or rainfall bins to created within the rainfall
%         range
%
% Optional Inputs:
% num - number of monte carlo simulations 
% R_max - maximum rainfall value to test (m/yr)
% R_min - minimum rainfall value to test (m/yr)
% R_step - rainfall step size for calculating (unitless)
% Klp_err - precent error on the partial coefficient of erosion (decimal
%           percentage, e.g. 0.2)
%
% Outputs:
% plotted comparision of erosional efficieneices from regression of
% observations and the stream-power model
%
%-------------------------------------------------------------------------%
% tashi delek!
%-------------------------------------------------------------------------%
%
% parse inputs
	p = inputParser;         
	p.FunctionName = 'K2K';
    %
	addRequired(p,'mar',@(x) isnumeric(x));
	addRequired(p,'E',@(x) isnumeric(x));
    addRequired(p,'E_err',@(x) isnumeric(x));
    addRequired(p,'ksn',@(x) isnumeric(x));
    addRequired(p,'ksn_err',@(x) isnumeric(x));
    addRequired(p,'m',@(x) isscalar(x));
    addRequired(p,'n',@(x) isnumeric(x));
    addRequired(p,'Klp',@(x) isscalar(x));
    addRequired(p,'nbins',@(x) isscalar(x));
    %
	addParameter(p,'num',1e6,@(x) isscalar(x));
	addParameter(p,'R_max',max(E),@(x) isscalar(x));
	addParameter(p,'R_min',min(ksn),@(x) isscalar(x));
    addParameter(p,'R_step',0.1,@(x) isscalar(x));
    addParameter(p,'Klp_err',0.2,@(x) isscalar(x));
    %
	parse(p,mar,E,E_err,ksn,ksn_err,m,n,Klp,nbins,varargin{:});
	mar = p.Results.mar;
	E = p.Results.E;
    E_err = p.Results.E_err;
    ksn = p.Results.ksn;
    ksn_err = p.Results.ksn_err;
    m = p.Results.m;
    n = p.Results.n;
    Klp = p.Results.Klp;
    nbins = p.Results.nbins;
    %
	num = p.Results.num;
	R_max = p.Results.R_max;
	R_min = p.Results.R_min;
    R_step = p.Results.R_step;
    Klp_err = p.Results.Klp_err;
%
% initialize variables
    MAR_data = [];
    K_data = [];
    size = [];
    Kl_array = [(Klp - Klp*Klp_err) Klp (Klp + Klp*Klp_err)];
    num_bins = fliplr(1:1:nbins);
    R_plot = R_min:R_step:R_max;
%
% calculate an envelope of predited erosion coefficients from the stream-
% power model
    for z = 1:length(Kl_array)
        for l = 1:length(R_plot)
                K_SP(z,l) = Kl_array(z)*R_plot(l)^m;
        end
    end
%    
for i = 1:length(num_bins)
    clear divisions K_LSE centers
	% find rainfall division points (bin edges)
        R_step(i) = (R_max - R_min)/num_bins(i);
        divisions = R_min:R_step(i):R_max;
        hold_bins = ones(1,num_bins(i))*R_step(i);
        size = [size hold_bins];
    %
    % divide data according to divsion points
    	group = discretize(mar,divisions);
    %
    % calculate bin centers based on the divisions
    	for k = 1:(length(divisions) - 1)
        	centers(k) = (divisions(k + 1) + divisions(k))/2; %#ok<*AGROW>
    	end
    %
    % for each bin, parse the data and regress these data to estimate
    % erosion coefficients
        for j = 1:length(centers)
            ind = find(group == j);
            %
            E_new = E(ind);
            E_err_new = E_err(ind);
            ksn_new = ksn(ind);
            ksn_err_new = ksn_err(ind);
            %
            % call MC_York to regress ksn - E data
                [C(j),~,~,~] = MC_York(n,E_new,ksn_new,E_err_new,ksn_err_new,num,'print','n');
                K_LSE(j) = (C(j).^-n)/1e6;
            %
            % concatonate data
                MAR_data = [MAR_data centers(j)];
                K_data = [K_data K_LSE(j)];
        end
end
%
% plot regression data
    figure(1)
    hold on
    yyaxis left
    scatter(MAR_data,K_data,size*50,size,'filled','MarkerEdgeColor',[0 0 0])
    colormap(flipud(baa_map))
% plot stream-power model
    yyaxis left
    line(R_plot,K_SP,'LineStyle','-','Color','k')
    xlabel('Mean annual rainfall (m/yr)')
    ylabel('Erosion coefficient')
% plot histogram data      
    num_bins = R_max/0.1;
    edges = (1:1:num_bins)*0.1;
    yyaxis right
    h = histogram(mar,edges);
    set(gca,'Ydir','reverse')
    h.FaceColor = [1 1 1];
    h.EdgeColor = [0 0 0];
    ylabel('Sample counts')
%