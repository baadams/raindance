function [coeff,expo,MSWD,MSWD_std] = MC_York(b,x,y,x_err,y_err,num,varargin)
%-------------------------------------------------------------------------%
% Function written by Byron A. Adams - Updated: 11 Apr 2019
%-------------------------------------------------------------------------%
%
% Description:
% this code uses the equations of York 1967 (CND J. Phys.) to impliment a 
% least-squares linear regression. specifically this script implements the
% equation for the slope (b) of a line where the weights of the x and y
% data are equal. this precondition is selected such that distance
% minimized between a data point and the regression line is orthogonal to
% the regression line. therefore, the equation for the slope is found in
% case iii (a) on page 1083. this solution also assumes no correlation of 
% uncertainties. see York 1967 (EPSL) for more discussion.
%
% power-law related data must be linearized before conducting a linear
% regression. to achieve this, data are log-transformed.
%
% because the regressed xy data do not likely have equal uncertainties,
% this least-squares regression routine is wrapped in a monte carlo
% protocol where a random, but normally distributed value is selected
% from each x and y value based on the mean and 1 sigma uncertainties for
% each sample. these regression simulations are carried out N times and
% mean slope values are calculated from which intercept values are 
% calculated.
%
% Usage:
% [coeff,n,MSWD,MSWD_std] = MC_York(b,x,y,x_err,y_err,varargin);
% [coeff,n,MSWD,MSWD_std] = MC_York(b,x,y,x_err,y_err,'num',1e5,...);
%
% Required inputs:
% b - the exponent of the power law. this can be a scalar or an array
% x - data to be plotted on the x axis. assumed to be erosion rate (m/Myr)
% y - variable to be plotted on the y axis. assumed to be channel steepness
% x_err - 1 sigma uncertainites on the x data
% y_err - 1 sigma uncertainites on the y data
%
% Optional inputs:
% num - number of monte carlo simulations
% x_max - maximum value of x to use in calculations
% fig_num - this is the figure number to be used
% colour - color of plotted outputs. [r g b] format where rgb are 0-1
% draw_fig - true/false flag to draw (true, default) or suppress drawing 
% fix_b - true/false flag to fix the exponent in the power law when fitting        
%
% Outputs:
% n - the best-fit power law exponent
% coeff - the best-fit coefficient for the power law
% MSWD - the mean square weighted deviation
% MSWD_std - 1 sigma uncertainty on the MSWD
%
%-------------------------------------------------------------------------%
% tashi delek
%-------------------------------------------------------------------------%
%
% parse inputs
	p = inputParser;         
	p.FunctionName = 'MC_York';
    %
	addRequired(p,'b',@(x) isnumeric(x));
	addRequired(p,'x',@(x) isnumeric(x));
    addRequired(p,'y',@(x) isnumeric(x));
    addRequired(p,'x_err',@(x) isnumeric(x));
    addRequired(p,'y_err',@(x) isnumeric(x));
	addRequired(p,'num',@(x) isnumeric(x));
    %
	addParameter(p,'x_max',max(x),@(x) isscalar(x));
    addParameter(p,'fig_num',1,@(x) isscalar(x));
    addParameter(p,'colour',[0 0 0],@(x) isnumeric(x));
    addParameter(p,'draw_fig',true,@(x) isscalar(x) && islogical(x));
    addParameter(p,'fix_b',true,@(x) isscalar(x) && islogical(x));
    %
	parse(p,b,x,y,x_err,y_err,num,varargin{:});
	b = p.Results.b;
	x = p.Results.x;
    y = p.Results.y;
    x_err = p.Results.x_err;
    y_err = p.Results.y_err;
	num = p.Results.num;
    %
	x_max = p.Results.x_max;
    fig_num = p.Results.fig_num;
    colour = p.Results.colour;
    draw_fig = p.Results.draw_fig;
    fix_b = p.Results.fix_b;
%
% initialize matricies
    x_box = zeros(num,length(x));
    y_box = zeros(num,length(x)); 
%    
% make normal distributions of x and y based on 1 sigma uncertainties
    for k = 1:length(x)
        rng('shuffle');
        x_box(:,k) = x_err(k)*randn(1,num) + x(k);
        rng('shuffle');
        y_box(:,k) = y_err(k)*randn(1,num) + y(k);
    end
%
% reshape the data and remove negative values
	new_x = reshape(x_box,[num*length(x),1]);
	new_y = reshape(y_box,[num*length(x),1]);
	rmv_x = find(new_x < 0);
	rmv_y = find(new_y < 0);
	rmv = [rmv_x;rmv_y];
	new_x(rmv) = [];
	new_y(rmv) = [];
%
% if the exponent is is not known (b has length > 1), start by calculating
% the best-fit exponent, which is the slope of the data once log-
% transformed. if b is fixed (b has length = 1) then only calculate the
% intercept
	if ~fix_b
        % initialize variables
            b_all = zeros(1,num);
            r = 0;
        	tol = 1e-6;
        % set maximum number of iterations
        	runs = 1e6;
        %
        % monte carlo loop. calculate regression parameters for each
        % simulation
        w1=waitbar(0,'Running iterations...');
        for j = 1:num
            % intialize chi_sq with some large number
                chi_sq = 1e6;
            %
            % make an array of x and y points from the distributions
                xi = transpose(x_box(j,:));
                yi = transpose(y_box(j,:));
            %
            % create equal weighting values for x and y (e.g. make them 1)
            	omega_x = (ones(length(x),1));
            	omega_y = (ones(length(x),1));
            %
            % linearize the randomly selected data
            	xi = log(xi);
            	yi = log(yi);
            %
            % intialize b with an approximation for slope
                b = (yi(1)-yi(length(yi)))/(xi(1)-xi(length(xi)));
            %
            % find the slope (b)and intercept (a). loop iterates b until
            % chi_square < tol
            for i = 1:runs
            	Wi = (omega_x.*omega_y)./((omega_y.*(b^2)) + omega_x);
                X_bar = sum(Wi.*xi)/sum(Wi);
                Y_bar = sum(Wi.*yi)/sum(Wi);
                U = xi - X_bar;
                V = yi - Y_bar;
                ai = sqrt(omega_x.*omega_y);
                b = sum(((Wi.^2).*V).*(U./omega_y+(V.*b./omega_x)-(V.*r./ai)))./sum(((Wi.^2).*U).*(U./omega_y+((V.*b)./omega_x)-((U.*r*b)./ai)));
                a = Y_bar - X_bar.*b;
                chi_sq_new = sum(((yi - a - (xi.*b)).^2)); 
                if chi_sq - chi_sq_new <= tol
                    break
                else
                    chi_sq = chi_sq_new;
                end
            end
            %
            % remove simulations containing negative slopes. this can be
            % created by sampling distributions beyond 2 standard
            % deviations.
                if isreal(b) == 1
                    b_all(j) = b;
                else
                    b_all(j) = NaN;
                end
            waitbar(j/num);
        end
        close(w1);
        % find mean b
            b = mean(b_all,'omitnan');
        %
        % linearize the entire synthetic mote carlo dataset
        	trans_x = log(new_x);
        	trans_y = log(new_y);
        %    
        % create equal weighting the entire synthetic mote carlo dataset
        	omega_x = ones(length(trans_x),1);
        	omega_y = ones(length(trans_x),1);
        %
        % solve for a given the solution for b above
        	Wi = (omega_x.*omega_y)./((omega_y.*(b^2)) + omega_x);
        	X_bar = sum(Wi.*trans_x)./sum(Wi);
        	Y_bar = sum(Wi.*trans_y)./sum(Wi);
            a = Y_bar - X_bar*b;
    else
        % this else loop solves for a if b is fixed
        % linearize the data
        	trans_x = log(new_x);
        	trans_y = log(new_y);
        %
        % create equal weighting values for x and y
        	omega_x = ones(length(trans_x),1);
        	omega_y = ones(length(trans_x),1);
        %
        % because the value of b is a provided input in this case, use only
        % the essential equations from York to solve for a
            Wi = (omega_x.*omega_y)./((omega_y.*(b^2)) + omega_x);
            X_bar = sum(Wi.*trans_x)/sum(Wi);
            Y_bar = sum(Wi.*trans_y)/sum(Wi);
            a = Y_bar - X_bar*b;
	end
%
% based on the best b (slope) and a (intercept) values, create synthetic
% x and y data for plotting the best-fit curve
	x_plot = 0:1:x_max;
	y_plot = exp(a).*x_plot.^b;
%        
% write data for output
	coeff = exp(a);
	expo = b;
%
% calculate the MSWD of the best-fit parameters
	n = 1/expo;
	C = 1/coeff^n;
	chi_sq = sum(((x - C.*(y.^n)).^2)./((x_err.^2) + (C.*n.*y.^(n - 1)).^2.*(y_err.^2)));
	MSWD = chi_sq/(length(x) - 2);
	MSWD_std = sqrt(2/(length(x) - 2));
%
% plot best-fit curve and parameters
    if draw_fig
    	figure(fig_num)
    	hold on
        e1=errorbar(x,y,y_err,y_err,x_err,x_err,'.');
        e1.Color=[0.5 0.5 0.5];
        e1.CapSize=1;
        scatter(x,y,30,[0.5 0.5 0.5],'filled');
    	line(x_plot,y_plot,'Color',colour,'LineWidth',2.5,'DisplayName','LSE')
        hold off
    end
%
end