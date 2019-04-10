function CRN_splits(mar,E,E_err,ksn,ksn_err,ksn_q,ksn_q_err)
% written by byron adams (09/2018)
%
% this script was designed to explore the climate influences on
% topography and erosion rates. specficially, the
% goal is to see if there are meaningful relationships between basin
% averaged normalized channel steepness (ksn) with basin averaged cosmogenic nuclide erosion    
% rates (E). this processes is repeated for different climate regimes, here  
% defined by selected bins of mean annual precipitation/rainfall. the  
% programm seeks to find the best-fit power law between E and ksn within climate bins. the stream-power law is also used to predict the
% ksn-E relationship given the mean annual precipitation bin value.
%
%
% REQUIRED SCRIPTS:
% LSE_MC - this script uses a monte carlo framework to find best-fit curves
%          from power law or threshold relationships
% stream_power - this script uses the stream power law to predict channel 
%                steepness values for a given precipitation value
% ploterr - plots x and y errorbars
%
% OPTIONAL SCRIPTS:

%
%
% INPUTS:
% mar - mean annual rainfall (mm/yr)
% E - erosion rate (mm/yr)
% E_err - 1 sigma uncertainites on the erosion rate data
% ksn - normalized channel steepness data
% ksn_err - 1 sigma uncertainites on the channel steepness data
% ksn-q - 
% ksn-q_err - 
%
% ***Note that there are model and plot tuning parameters that may need to
% be changed below.
%
%
% OUTPUTS:
% plotted modeled morphometric-E relationship data
%
% tashi delek!
%-------------------------------------------------------------------------%
%
% ***BEGIN MODEL TUNING VARIABLES***
% water discharge exponent for stream-power law (unitless)
    m = 1;
% slope exponent for stream-power law (unitless)
    n = 2.22;
% lithology specific detatchment efficiency (L^1-3m T^m-1). Note that L is
% in meters and T is in years
    Kl = 1.9e-9;
% precipitation division points (bin edges)
    divisions = [0 1.5 2.5 3.5 6];
% power-law exponents to test for ksn vs E regressions
    %pksn = 0.3:0.01:1;
    pksn = 1/n;
% power-law exponents to test for R vs E regressions
    pR = 1/n;
% number of monte carlo simulations 
    N = 1e6;
% ***END MODEL TUNING VARIABLES***
%
% ***PLOT TUNING VARIABLES***
% set the x axis limit (m/Myr)
    E_max = 4000;
% set the relief axis limit (m)
    R_max = 2500;
% set the ksn axis limit
    ksn_max = 400;
% set the slope axis limit (m)
    S_max = 50;
% ***END PLOT TUNING VARIABLES***
%
%
% makes sure that all variables are in units of meters (assumes input data
% are output from Adam's TAK program
    E = E*1000;
    E_err = E_err*1000;
    mar = mar/1000;
%
% calcualte the size of each MAR bin
    for p = 1:length(divisions)-1
        bin_size(p) = divisions(p + 1) - divisions(p);
    end
% divide data according to divsion points
    group = discretize(mar,divisions);
%
% calculate bin centers based on the divisions
    for k = 1:(length(divisions) - 1)
        centers(k) = (divisions(k + 1) + divisions(k))/2; %#ok<*AGROW>
    end
%
% assign each data point a bin center based on the bin edges above
    for i = 1:length(group)
        mar_group(i) = centers(group(i));
    end
%
% calculate erosional effeciancy values of the stream-power law in (m^-1
% Myr^-1)
    K_SP = (1./(Kl.*1e6.*centers)).^(1/n);
%
% plot a histogram to see the distribution of data into selected bins
    figure(1)
    histogram(mar,divisions);
    ylabel('Counts')
    xlabel('Mean annual rainfall bins (m/yr)')
    xlim([0 max(divisions)])
    txt = ['Total counts = ' num2str(length(group))];
    legend(txt)
%
% based on the precipitation groupings split up the erosion and steepness 
% data and pass it to a regression protocol
    for j = 1:max(group)
        ind = find(group == j);
        %
        E_new = E(ind);
        E_err_new = E_err(ind);
        ksn_new = ksn(ind);
        ksn_err_new = ksn_err(ind);
        mar_grp_new = mar_group(ind); %#ok<*NASGU>
        mar_new = mar(ind);
        %
        figure(1)
        hold on
        num_bins = bin_size(j)/0.1;
        edges = [1:1:num_bins]*0.1 + divisions(j); %#ok<NBRAK>
        histogram(mar_new,edges)
        %
        % ***E vs ksn data processing and curve fitting***
            % call LSE_MC to calculate the best power-law realtionship
            % for the given division of ksn values
                [K_LSE(j),para(j)] = LSE_MC(pksn,E_new,ksn_new,E_err_new,ksn_err_new,N,E_max,ksn_max,j+2,[0.1300 0.1100 0.7750 0.8150],[0.5 0.5 0.5]);
            %
            % call model_curve to predict the power-law relationship for
            % the given division of channel steepness values and the 
            % precipitation 
                stream_power(m,n,Kl,mar_grp_new(1),E_max,j+2,[1 0 1]);
            %
            % calculate the MSWD of the stream-power law
                chi_sq = sum(((E_new - (1/K_SP(j)^n).*(ksn_new.^n)).^2)./((E_err_new.^2) + ((1/K_SP(j)^n).*n.*ksn_new.^(n - 1)).^2.*(ksn_err_new.^2)));
                MSWD = chi_sq/(length(E_new) - 2);
                MSWD_std = sqrt(2/(length(E_new) - 2));
                %
                txt = {'SPM',['C = ' num2str(round(K_SP(j),2,'significant'))],['MSWD = ' num2str(round(MSWD,2)) ' +/- ' num2str(round(MSWD_std*2,2))]};
                text(E_max*0.5,ksn_max*0.5,txt,'FontSize',8)
            %
            % make the plots look good
                figure(j+2)
                hold on
                ploterr(E_new,ksn_new,E_err_new*2,ksn_err_new*2,'.k','hhx',0.1,'hhy',0.1)
                scatter(E_new,ksn_new,30,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0 0 0])
                %
                xlabel('Erosion rate (m/Myr)')
                xlim([0 E_max])
                ylabel('k_s_n (m^{\0.9})')
                ylim([0 ksn_max]);
                legend({'MSWD power' 'LSE MC power' 'stream power'});
    end
    %
    % plot 
        figure(j+4)
        hold on
        subplot('Position',[0.12 0.6 0.78 0.35])
        h3 = gscatter(E,ksn,mar_group,[],'.',20);
        legend('off')
        xlabel('Erosion rate (m/Myr)')
        xlim([0 E_max])
        ylabel('k_s_n (m^{\0.9})')
        ylim([0 max(ksn_q)])
        %
        for l = 1:max(group)
            cgroup = h3(l);
            ind = find(group == l);
            E_new = E(ind);
            E_err_new = E_err(ind);
            ksn_new = ksn(ind);
            ksn_err_new = ksn_err(ind);
            mar_grp_new = mar_group(ind);
            new_mar = mar(ind);
            %
            % call LSE_MC to calculate the best power-law realtionship
            % for the given division of ksn values
                LSE_MC(pksn,E_new,ksn_new,E_err_new,ksn_err_new,N,E_max,ksn_max,j+4,[0.12 0.6 0.78 0.35],cgroup.Color);
                stream_power(m,n,Kl,mar_grp_new(1),E_max,j+4,cgroup.Color);
        end
        %
        % call LSE_MC to calculate the best power-law realtionship
        % for the given division of ksn values
        	LSE_MC(pksn,E,ksn_q,E_err,ksn_q_err,N,E_max,ksn_max+200,j+4,[0.12 0.12 0.78 0.35],[0.5 0.5 0.5]);
        % 
        % use Kl to predict ksn_q for range of E
            E_model = 0:10:E_max;
            ksn_q_model = (E_model./(Kl*1e6)).^(1/n);
        %
        % calculate the MSWD of the stream-power law
        	chi_sq = sum(((E - (Kl*1e6).*(ksn_q.^n)).^2)./((E_err.^2) + (Kl*1e6.*n.*ksn_q.^(n - 1)).^2.*(ksn_q_err.^2)));
        	MSWD = chi_sq/(length(E) - 2);
        	MSWD_std = sqrt(2/(length(E) - 2));
            %
            txt = {'SPM',['C = ' num2str(round((1/(Kl*1e6))^(1/n),2,'significant'))],['MSWD = ' num2str(round(MSWD,2)) ' +/- ' num2str(round(MSWD_std*2,2))]};
            text(E_max*0.5,ksn_max,txt,'FontSize',8)
        %
        subplot('Position',[0.12 0.12 0.78 0.35])
        h3 = gscatter(E,ksn_q,mar_group,[],'.',20);
        line(E_model,ksn_q_model,'Color',[0 0 0],'LineWidth',1)
        legend('off')
        xlabel('Erosion rate (m/Myr)')
        xlim([0 E_max])
        ylabel('k_s_n- q')
        ylim([0 600])
        %
        for l = 1:max(group)
            cgroup = h3(l);
            ind = find(group == l);
            E_new = E(ind);
            E_err_new = E_err(ind);
            ksn_new = ksn_q(ind);
            ksn_err_new = ksn_err(ind);
            %
            % call LSE_MC to calculate the best power-law realtionship
            % for the given division of ksn values
                LSE_MC(pksn,E_new,ksn_new,E_err_new,ksn_err_new,N,E_max,ksn_max,j+4,[0.12 0.12 0.78 0.35],cgroup.Color);
        end
    % plot
            % call LSE_MC to calculate the best power-law realtionship
            % for the given division of ksn values
                LSE_MC(pksn,E,ksn,E_err,ksn_err,N,E_max,ksn_max,j+5,[0.1300 0.1100 0.7750 0.8150],[0.5 0.5 0.5]);
            %
            % use Kl to predict ksn_q for range of E
                ksn_model = (E_model./(Kl*1e6*median(mar))).^(1/n);
            %
            % calculate the MSWD of the stream-power law
                chi_sq = sum(((E - (Kl*1e6*median(mar)).*(ksn.^n)).^2)./((E_err.^2) + ((Kl*1e6*median(mar)).*n.*ksn.^(n - 1)).^2.*(ksn_err.^2)));
                MSWD = chi_sq/(length(E) - 2);
                MSWD_std = sqrt(2/(length(E) - 2));
                %
                txt = {'SPM',['C = ' num2str(round((1/(Kl*1e6*median(mar)))^(1/n),2,'significant'))],['MSWD = ' num2str(round(MSWD,2)) ' +/- ' num2str(round(MSWD_std*2,2))]};
                text(E_max*0.5,ksn_max*0.9,txt,'FontSize',8)
        %
        figure(j+5)
        hold on
        ploterr(E,ksn,E_err,ksn_err,'.k','hhx',0.1,'hhy',0.1)
        scatter(E,ksn,30,mar,'filled');
        line(E_model,ksn_model,'Color',[0 0 0],'LineWidth',1)
        %
        xlabel('Erosion rate (m/Myr)')
        xlim([0 E_max])
        ylabel('k_s_n (m^{\0.9})')
        ylim([0 600])
        caxis([0 max(mar)])
        colormap(mar_map)
        colorbar
%