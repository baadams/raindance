function stream_power(m,n,Kl,R,E_max,fig_num,colour)
%-------------------------------------------------------------------------%
% Function written by Byron A. Adams - Updated: 11 Apr 2019
%-------------------------------------------------------------------------%
%
% Description:
% Function uses the stream-power model to predict channel steepness values
% for a given rainfall value. the form of the stream-power model used
% here comes from Royden and Perron (2013; JGR ES). This form of the
% equation assumes steady-state topography for a given rock uplift rate.
%
% Usage:
% [stream_power(m,n,Kl,P,E_max,fig_num,colour);
%
% Required Inputs:
% m - discharge exponent for stream-power model (unitless)
% n - slope exponent for stream-power model (unitless)
% Kl - lithology specific detatchment efficiency (L^1-3m T^m-1) where L is
%      in meters and T is in years
% R - mean annual rainfall rate (L T^-1)
% E_max - this is the maximum erosin/rock uplift rate to test (L T^-1). 
%         rates are expected to be in meters per million years.
% fig_num - this is the figure number where the output is to be plotted
% colour - this is the color of the plotted line needs to be in [r g b]
%          format where rgb values vary between 0 and 1
%
% Outputs:
% plotted modeled ksn-E relationship
%
%-------------------------------------------------------------------------%
% tashi delek
%-------------------------------------------------------------------------%
% 
% make an array of rock uplift rates to calculate channel steepness
    U = 0:10:E_max;
%
% intiatialize channel steepness array
    ksn = zeros(1,length(U));
%
% calculate steady-state channel steepness values at different rock uplift
% rates. adapted from eq (2) from Royden and Perron 2013 (JGR ES)
	for i = 1:length(U)
        ksn(i) = ((U(i)/1e6)/(Kl*(R^m)))^(1/n);
	end
%
% plot the modeled ksn values
	figure(fig_num)
	hold on
	line(U,ksn,'Color',colour,'LineWidth',1,'DisplayName','SPM')
%