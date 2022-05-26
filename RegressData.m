% RegressData.m
% Sarah West
% 4/20/22

% Abstracted regression using "regress"
% Inputs: parameters.explanator, parameters.response,
% parameters.predictorsDim

function [parameters] = RegressData(parameters)
    
    % If there's a "values" field from RunAnalysis, print updating message
    % for user. 
    if isfield(parameters, 'values')
        message = ['Regressing '];
        for dispi = 1:numel(parameters.values)/2
           message = [message ', ' parameters.values{dispi}];
        end
        disp(message); 
    end

    % Pull these out so you don't accidentally edit original values.
    response = parameters.response;
    explanatory = parameters.explanatory;

    % If the different predictors aren't already different columns, flip
    % matrix.
    if parameters.predictorsDim ~= 2
        explanatory = explanatory'; 
    end 

    % Make response match dimensions of explanatory.
    if size(response,1) ~= size(explanatory,1)
       response = response';
    end 

    % If uer said to & last column isn't all 1s, add a column of 1s for intercept.
    if isfield(parameters, 'useIntercept') && ~parameters.useIntercept
        % Don't do anything
    else
        if any(explanatory(:,end) ~= 1)

            explanatory = [explanatory ones(size(explanatory,1),1)]; 

        end
    end

    % Run regression.
    [b,~,~,~,stats] = regress(response, explanatory);

    % Put into output
    parameters.results_r2 = stats(1);
    parameters.results_betas = b; 
end