% Histogram.m
% Sarah West
% 4/19/22

% Is "histogram" function, but compatible with RunAnalysis

function [parameters] = Histogram(parameters)

    % If there's a "values" field from RunAnalysis, print updating message
    % for user. 
    if isfield(parameters, 'values')
        message = ['Histogram for '];
        for dispi = 1:numel(parameters.values)/2
           message = [message ', ' parameters.values{dispi}];
        end
        disp(message); 
    end
    
    % Find value ranges for bottom of histogram. 
    dx = (parameters.histogram_value_range(2) - parameters.histogram_value_range(1)) ./ parameters.histogram_nBins;
    x = (parameters.histogram_value_range(1) + dx): dx : parameters.histogram_value_range(2); 

    % Don't allow figures to pop up.
    set(groot,'defaultFigureVisible','off');

    % Calculate histogram
    h = histogram(parameters.data, [parameters.histogram_value_range(1) x]);

    % Re-allow figures to pop up.
    set(groot,'defaultFigureVisible','on');

    % Put value ranges and values into the output. 
    parameters.histogram = [x' h.Values'];

end 