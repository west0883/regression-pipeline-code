% Histogram.m
% Sarah West
% 4/19/22

% Is "histogram" function, but compatible with RunAnalysis

function [parameters] = Histogram(parameters)

    parameters.histogram = histogram(parameters.data, parameters.nbins);
    
end 