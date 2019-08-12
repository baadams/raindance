function [RMSE,min_RMSE,n_best,Klp_best] = find_Klp(mar,E,E_err,ksn,ksn_err,m,n,Klp,nbins,varargin)
%-------------------------------------------------------------------------%
% Function written by Byron A. Adams - Updated: 3 Aug 2019
%-------------------------------------------------------------------------%
%
% Description:
% Function finds a rainfall modulated erosional efficiency and n value that
% best fit the regressions of observed data, and the stream-power model
% for the same rainfall rate. From the rainfall modulated erosional
% efficiency value (K), an other erosional efficiency, this time a
% coeffcient lacking the influence of rainfall (Klp) can be calciulated.
% The script compares the erosional efficiencies predicted from regressions 
% of channel steepness and erosion rate data and erosional efficiencies
% predcited from the stream-power model. Before regressions are carried
% out, the data are parsed based on the creation of mean annual rainfall
% bins. The rainfall bins are calculated based on R_max, R_min and nbins,
% where nbins is the maximum number of bins that will be used to divded the
% R_max - R_min range. The function then calculates further regressions
% where the rainfall range is divided by nbins - 1, nbins - 2, and so on 
% until nbins = 1. This method creates a large number of semi-random bins 
% of varying size and center value. This process is repeated for a
% combintation of n and Klp values provided by the user. The root mean
% square error is calculated from the residuals of the regression- and
% stream-power model-based erosional efficiency values. The best-fit n and
% Klp values are found where the RMSE is minimized.
%
% Usage:
% find_Klp(mar,E,E_err,ksn,ksn_err,ksn_q,ksn_q_err,m,n,Klp,nbins);
% [RMSE,min_RMSE,n_best,Klp_best] = find_Klp(mar,E,E_err,ksn,ksn_err,ksn_q,ksn_q_err,m,n,Klp,nbins,'num',1e6...);
%
% Required Inputs:
% mar - mean annual rainfall (m/yr)
% E - erosion rate (m/Myr)
% E_err - 1 sigma uncertainites on the erosion rate data (m/Myr)
% ksn - normalized channel steepness data
% ksn_err - 1 sigma uncertainites on the channel steepness data
% m - discharge exponent for stream-power law (unitless)(scalar)
% n - an array of slope exponents for stream-power law (unitless)
% Klp - an array of partial coefficient of erosion (L^1-3m T^m-1). Note 
%       that L is in meters and T is in years, and m is the discharge 
%       exponent
% nbins - maximum number or rainfall bins to created within the rainfall
%         range
%
% Optional Inputs:
% num - number of monte carlo simulations used in regressions
% R_max - maximum rainfall value to test (m/yr)
% R_min - minimum rainfall value to test (m/yr)
%
% Outputs:
% RMSE - root mean square error
% min_RMSE - minimum root mean square error
% n_best - the value of n that yields the minimum root mean square error
% Klp_best - the Klp value that yields the minimum root mean square error
%
%-------------------------------------------------------------------------%
% tashi delek!
%-------------------------------------------------------------------------%
%
% parse inputs
	p = inputParser;         
	p.FunctionName = 'find_Klp';
    %
	addRequired(p,'mar',@(x) isnumeric(x));
	addRequired(p,'E',@(x) isnumeric(x));
    addRequired(p,'E_err',@(x) isnumeric(x));
    addRequired(p,'ksn',@(x) isnumeric(x));
    addRequired(p,'ksn_err',@(x) isnumeric(x));
    addRequired(p,'m',@(x) isscalar(x));
    addRequired(p,'n',@(x) isnumeric(x));
    addRequired(p,'Klp',@(x) isnumeric(x));
    addRequired(p,'nbins',@(x) isscalar(x));
    %
	addParameter(p,'num',1e6,@(x) isscalar(x));
	addParameter(p,'R_max',max(E),@(x) isscalar(x));
	addParameter(p,'R_min',min(ksn),@(x) isscalar(x));
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
%
% initialize variables
    MAR_data = [];
    num_bins = fliplr(1:1:nbins);
%  
% loop through each number of bin centers
	for i = 1:length(num_bins)
        clear divisions centers
        % find rainfall division points (bin edges)
        	R_step = (R_max - R_min)/num_bins(i);
            divisions = R_min:R_step:R_max;
        %
        % divide data according to divsion points
            group(i,:) = discretize(mar,divisions);
        %
        % calculate bin centers based on the divisions
            for k = 1:(length(divisions) - 1)
                centers(k) = (divisions(k + 1) + divisions(k))/2; %#ok<*AGROW>
            end
        %
        % concatonate data
            MAR_data = [MAR_data centers];
	end
%
% loop through n values
	for q = 1:length(n)
        K_data = [];
        % loop through each bin center
        for j = 1:nbins
                clear K_LSE
                for i = 1:nbins
                    ind = find(group(j,:) == i);
                    E_new = E(ind);
                    E_err_new = E_err(ind);
                    ksn_new = ksn(ind);
                    ksn_err_new = ksn_err(ind);
                    %
                    % call MC_York to regress ksn - E data
                        [C,~,~,~] = MC_York(n(q),E_new,ksn_new,E_err_new,ksn_err_new,num,'print','n');
                        K_LSE = (C.^-n(q))/1e6;
                    % save data
                        K_data = [K_data K_LSE];
                end
         end
            K_data(isnan(K_data)) = [];
            C_new(q,:) = K_data;
	end
%   
% loop through n values
	for q = 1:length(n)
    % loop through Kl values
        for l = 1:length(Klp)
            % find the root mean square estimator for the n and Kl combination
                res_sq = (C_new(q,:) - (1./(Klp(l).*1e6.*MAR_data).^m).^(1./n(q))).^2;
                RMSE(q,l) = sqrt((1/length(MAR_data))*sum(res_sq,'omitnan'));
        end
	end
%
% find the minimum RSME value
    min_RMSE = min(min(RMSE));
    [r,c] = find(RMSE == min_RMSE);
%
% find the n and Klp values that yield the minimum RMSE value
    n_best = n(r);
    Klp_best = Klp(c);
%
% plot RSME data
    figure(1)
    hold on
    tours = 0.1:0.1:20;
    [~,c] = contourf(RMSE,tours);
    set(gca,'Ydir','reverse')
    c.LineStyle = 'none';
    %
    colormap(baa_map)
    h = colorbar;
    ylabel(h,'RMSE');
    caxis([min(min(RMSE)) max(tours)])
    %
    ticks_x = linspace(0,(length(Klp)),length(Klp)/4);
    ticks_y = linspace(0,length(n),length(n));
    x_ticklabels = cellfun(@(v) sprintf('%2.1d',v),num2cell(Klp),'UniformOutput',false);
    y_ticklabels = cellfun(@(v) sprintf('%2.2g',v),num2cell(n),'UniformOutput',false);
    %
    set(gca,'XTick',ticks_x)
    set(gca,'XTickLabels',x_ticklabels)
    set(gca,'YTick',ticks_y)
    set(gca,'YTickLabels',y_ticklabels)
    ax =  gca;
    ax.Layer = 'top';
    box on
    xlabel('Klp')
    ylabel('n')
    %
    min_tour = ceil(min(min(RMSE)));
    [M,~] = contour(1:1:length(Klp),1:1:length(n),RMSE,[min_tour,min_tour]);
    plot(M(1,2:end),M(2,2:end),'-w','LineWidth',1.5)
%