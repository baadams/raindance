function gimmie_m(mar,E,E_err,ksn,ksn_err,n,P_min,P_max,nbins,varargin)
%-------------------------------------------------------------------------%
% Function written by Byron A. Adams - Updated: 3 Jun 2019
%-------------------------------------------------------------------------%
%
% Description:
% Function finds ...
%
% Usage:
% K2K(mar,E,E_err,ksn,ksn_err,ksn_q,ksn_q_err,m,n,Klp,divisions);
%
% Required Inputs:
% mar - mean annual rainfall/precipitation (m/yr)
% E - erosion rate (m/Myr)
% E_err - 1 sigma uncertainites on the erosion rate data (m/Myr)
% ksn - normalized channel steepness data
% ksn_err - 1 sigma uncertainites on the channel steepness data
% m - water discharge exponent for stream-power law (unitless)
% n - slope exponent for stream-power law (unitless)
% Klp - partial coefficient of erosion (L^1-3m T^m-1). Note that L is
%       in meters and T is in years
%
% Optional Inputs:
% num - number of monte carlo simulations 
% p_plot_step -
%
% Outputs:
% plotted modeled morphometric-E relationship data
%
%-------------------------------------------------------------------------%
% tashi delek
%-------------------------------------------------------------------------%
%
    num = 1e6;
%
% calculate the distriubtion of number of bins to test
    num_bins = fliplr(1:1:nbins);
%
% initialize variables
    MAR_data = [];
    K_data = [];
%  
for l = 1:legnth(n)
    for i = 1:length(num_bins)
        clear divisions K_LSE centers
        % precipitation division points (bin edges)
            P_step = (P_max - P_min)/num_bins(i);
            divisions = P_min:P_step:P_max;
        %
        % divide data according to divsion points
            group = discretize(mar,divisions);
        %
        % calculate bin centers based on the divisions
            for k = 1:(length(divisions) - 1)
                centers(k) = (divisions(k + 1) + divisions(k))/2; %#ok<*AGROW>
            end
        %
        for j = 1:length(centers)
            ind = find(group == j);
            %
            E_new = E(ind);
            E_err_new = E_err(ind);
            ksn_new = ksn(ind);
            ksn_err_new = ksn_err(ind);
            %
            % ***E vs ksn data processing and curve fitting***
                [C_LSE(j),~,~,~] = MC_York(1./n(l),E_new,ksn_new,E_err_new,ksn_err_new,num);
                close(1)
                K_LSE(j) = (C_LSE(j).^-n(l))/1e6;
                %
            % save data
                MAR_data = [MAR_data centers(j)];
                K_data = [K_data K_LSE(j)];
        end
    end
%
    K = transpose(K_data);
    R = transpose(MAR_data);
    trans_x = log(R);
    trans_y = log(K);
    pairs = [trans_x trans_y];
    meanX = mean(pairs,1);
    [coeff,score,~] = pca(pairs);
    Xfit = repmat(meanX,length(R),1) + score(:,1)*coeff(:,1)';
    p = polyfit(Xfit(:,1),Xfit(:,2),1);
    m(l) = p(1);
end
%
% plot
    figure(1)
    plot(n,m,'.k')
    xlabel('n')
    ylabel('m')
%