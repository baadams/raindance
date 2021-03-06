function [E,ksn] = SPM(m,n,Klp,R,E_max)
%-------------------------------------------------------------------------%
% Function written by Byron A. Adams - Updated: 11 Apr 2019
%-------------------------------------------------------------------------%
%
% Description:
% Function uses the stream-power model to predict channel steepness values
% for a given rainfall value. The form of the stream-power model used
% here comes from Royden and Perron (2013; JGR ES). This form of the
% equation assumes steady-state topography for a given rock uplift rate.
%
% Usage:
% [E,ksn] = SPM(m,n,Klp,P,E_max);
%
% Required Inputs:
% m - discharge exponent for stream-power model (unitless)
% n - slope exponent for stream-power model (unitless)
% Klp - lithology specific detatchment efficiency (L^1-3m T^m-1) where L is
%       in meters and T is in years
% R - mean annual rainfall rate (L T^-1)
% E_max - this is the maximum erosin/rock uplift rate to test (L T^-1). 
%         rates are expected to be in meters per million years.
%
% Outputs:
% modeled ksn-E relationship
%
%-------------------------------------------------------------------------%
% tashi delek
%-------------------------------------------------------------------------%
% 
% parse inputs
    p = inputParser;         
	p.FunctionName = 'SPM';
    
    addRequired(p,'m',@(x) isscalar(x));
    addRequired(p,'n',@(x) isnumeric(x));
    addRequired(p,'Klp',@(x) isscalar(x));
	addRequired(p,'R',@(x) isscalar(x));
	addRequired(p,'E_max',@(x) isscalar(x));
    
	parse(p,m,n,Klp,R,E_max);
    m = p.Results.m;
    n = p.Results.n;
    Klp = p.Results.Klp;
	R = p.Results.R;
	E_max = p.Results.E_max;   

% make an array of rock uplift rates to calculate channel steepness
    E = 0:10:E_max;

% intiatialize channel steepness array
    ksn = zeros(1,length(E));

% calculate steady-state channel steepness values at different rock uplift
% rates. adapted from eq (2) from Royden and Perron 2013 (JGR ES)
	for i = 1:length(E)
        ksn(i) = ((E(i)/1e6)/(Klp*(R^m)))^(1/n);
	end
