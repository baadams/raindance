% raindance
% Version 1.0  20-July-2020
%
% These are a series of MATLAB functions that were designed to explore  
% climate influences on topography and erosion rates. Specficially, the
% goal is to see if there are meaningful relationships between
% basin-averaged normalized channel steepness (ksn) with basin-averaged     
% cosmogenic nuclide erosion rates. This processes is repeated for   
% different climate regimes, here defined by selected bins of mean annual 
% rainfall. Each function contains a header with basic  functionality
% functionality info along with expected inputs and possible outputs. 
%
% If you encounter errors or have suggestions please contact:
% Byron A. Adams
% byron.adams@bristol.ac.uk
%
% Paper documenting the tools is now published in ___:
%
% URL: 
%
% If you use or modify these tools for use in a publication, please
% cite the above paper.
%
% Required MATLAB Toolbox list:
% None
%
% To fully utilize all of raindance, it is recommended that you have
% at least MATLAB 9.5 installed.
%
%
% Function List: 
% besties - finds a the partial erosional efficencny (Klp), n and, m values
%           that best fit the observed ksn and E data, and the stream-
%           power model
% bin_it - explores the sensitivity of fluvial relief and erosion rates to
%          rainfall rates
% MC_York - this script uses a monte carlo framework to find best-fit 
%           curves from power law or threshold relationships
% SPM - this script uses the stream power model to predict channel 
%       steepness values for a given precipitation value
% by_map - a colormap grading from blue to green to yellow
% mar_map - a colormap grading from blue to white to red
