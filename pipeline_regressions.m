% pipeline_regression.m
% Sarah West
% 4/19/22

% Regressions(from locomotion paper):
% regression_against_treadmill_withtimelags_differentlags.m
% regression_against_treadmill_individuallags.m
% regression_against_paws.m
% regression_against_paws_individuallags.m

%% Initial Setup  
% Put all needed paramters in a structure called "parameters", which you
% can then easily feed into your functions. 
clear all; 

% Create the experiment name.
parameters.experiment_name='Random Motorized Treadmill';

% Output directory name bases
parameters.dir_base='Y:\Sarah\Analysis\Experiments\';
parameters.dir_exper=[parameters.dir_base parameters.experiment_name '\']; 

% Load mice_all, pass into parameters structure
load([parameters.dir_exper '\mice_all.mat']);
parameters.mice_all = mice_all;

% ****Change here if there are specific mice, days, and/or stacks you want to work with**** 
parameters.mice_all = parameters.mice_all(1);

periods_spontaneous = {'rest';'walk';'startwalk';'prewalk';'stopwalk';'postwalk';'full_onset';'full_offset'};

load([parameters.dir_exper 'periods_nametable.mat']);
periods_motorized = periods;

%% Spontaneous-- Look at histograms of response variable by itself 
% (each correlation, not dependent on velocity yet) Is important for
% picking the style of GLM.
number_of_sources = 29; 
indices = find(tril(ones(number_of_sources), -1));

period = 'walk';

% Always clear loop list first. 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterations.
parameters.loop_list.iterators = {'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
              'index', {'loop_variables.indices'}, 'index_iterator'};

parameters.loop_variables.mice_all = parameters.mice_all;
parameters.loop_variables.indices = indices;

% Number of bins to use in the histogram output.
parameters.histogram_nBins = 200;
parameters.histogram_value_range = [-1 1];

% Things to load
% Correlations 
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'fluorescence analysis\correlations\spontaneous\'] 'mouse', '\instances reshaped\'};
parameters.loop_list.things_to_load.data.filename= {'correlations_', period, '.mat'};
parameters.loop_list.things_to_load.data.variable= {'correlations_reshaped(', 'index', ', :)'}; 
parameters.loop_list.things_to_load.data.level = 'mouse';

% Things to save
% histogram values
parameters.loop_list.things_to_save.histogram.dir = {[parameters.dir_exper 'regression analysis\walk velocity\response distributions\spontaneous\'], 'mouse', '\'};
parameters.loop_list.things_to_save.histogram.filename= {'correlations_response_distribution_histograms_', period, '.mat'};
parameters.loop_list.things_to_save.histogram.variable= {'histograms{', 'index_iterator', ',1}'}; 
parameters.loop_list.things_to_save.histogram.level = 'mouse';

RunAnalysis({@Histogram}, parameters);

%% 

distributions = NaN(size(correlations_reshaped,2), 36); % 36 is the  default number of bins 
 a = squeeze(correlations(1,10,:));
[vels, places] = sort(average_velocity);
a = a(places);

%% Regress all spontaneous
number_of_sources = 29; 
indices = find(tril(ones(number_of_sources), -1));

period = 'walk';

% Always clear loop list first. 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterations.
parameters.loop_list.iterators = {'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
              'index', {'loop_variables.indices'}, 'index_iterator'};

parameters.loop_variables.mice_all = parameters.mice_all;
parameters.loop_variables.indices = indices;

% Dimension of different predictors (so you can orient variables correctly
% before regressing)
parameters.predictorsDim = 1; 

% Input 
% Correlations as response/ dependent variable.
parameters.loop_list.things_to_load.response.dir = {[parameters.dir_exper 'fluorescence analysis\correlations\spontaneous\'] 'mouse', '\instances reshaped\'};
parameters.loop_list.things_to_load.response.filename= {'correlations_', period, '.mat'};
parameters.loop_list.things_to_load.response.variable= {'correlations_reshaped(', 'index', ', :)'}; 
parameters.loop_list.things_to_load.response.level = 'mouse';

% Velocities as explanatory/independent variable.
parameters.loop_list.things_to_load.explanatory.dir = {[parameters.dir_exper 'behavior\spontaneous\concatenated velocity by behavior\'], 'mouse', '\'};
parameters.loop_list.things_to_load.explanatory.filename= {'average_velocity_byinstance_', period, '.mat'};
parameters.loop_list.things_to_load.explanatory.variable= {'average_velocity'}; 
parameters.loop_list.things_to_load.explanatory.level = 'mouse';

% Output
parameters.loop_list.things_to_save.results.dir = {[parameters.dir_exper '\regression analysis\walk velocity\spontaneous\'], 'mouse', '\'};
parameters.loop_list.things_to_save.results.filename= {'regression_results_', period, '.mat'};
parameters.loop_list.things_to_save.results.variable= {'regression_results{', 'index_iterator',', 1}'}; 
parameters.loop_list.things_to_save.results.level = 'mouse';

RunAnalysis({@RegressData}, parameters);


%% Make vector of all motorized velocities 
condition = 'c_walk';
periods = {'176'; '177'; '178'; '179'};

% Always clear loop list first. 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterations.
parameters.loop_list.iterators = {'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
              'period', {'loop_variables.periods(:)'}, 'period_iterator'};
parameters.loop_variables.mice_all = parameters.mice_all;
parameters.loop_variables.periods = periods;

% instructions for replicating.
parameters.toReplicate = {'parameters.period_velocities(', 'period_iterator', ')'};
parameters.replicateDims = {'[size(parameters.data, 2), 1]'};

% Velocities (cm/s) for each motor speed.
period_velocities = [2.25, 2.77, 3.33, 3.88]; 
parameters.period_velocities = period_velocities;

% Concatenation dimension
parameters.concatDim = 1;

% Inputs (reshaped correlations, for the number)
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'fluorescence analysis\correlations\motorized\'] 'mouse', '\instances reshaped\'};
parameters.loop_list.things_to_load.data.filename= {['correlations_' condition '_'], 'period', '.mat'};
parameters.loop_list.things_to_load.data.variable= {'correlations_reshaped'}; 
parameters.loop_list.things_to_load.data.level = 'period';

% Output
parameters.loop_list.things_to_save.concatenated_data.dir = {[parameters.dir_exper 'regression analysis\walk velocity\motorized\velocity vectors\'], 'mouse','\'};
parameters.loop_list.things_to_save.concatenated_data.filename= {'velocity_vector.mat'};
parameters.loop_list.things_to_save.concatenated_data.variable= {'velocity_vector'}; 
parameters.loop_list.things_to_save.concatenated_data.level = 'mouse';

% Items to rename between functions.
parameters.loop_list.things_to_rename = { {'data_replicated', 'data'}};

RunAnalysis({@ReplicateData, @ConcatenateData}, parameters);

%% Regress motorized

number_of_sources = 29; 
indices = find(tril(ones(number_of_sources), -1));

condition = 'c_walk';
periods = {'176', '177', '178', '179'};

% Always clear loop list first. 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterations.
parameters.loop_list.iterators = {'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
              'index', {'loop_variables.indices'}, 'index_iterator'};

parameters.loop_variables.mice_all = parameters.mice_all;
parameters.loop_variables.indices = indices;

% Dimension of different predictors (so you can orient variables correctly
% before regressing)
parameters.predictorsDim = 1; 

% Input 
% Correlations as response/ dependent variable.
parameters.loop_list.things_to_load.response.dir = {[parameters.dir_exper 'fluorescence analysis\correlations\spontaneous\'] 'mouse', '\instances reshaped\'};
parameters.loop_list.things_to_load.response.filename= {'correlations_', period, '.mat'};
parameters.loop_list.things_to_load.response.variable= {'correlations_reshaped(', 'index', ', :)'}; 
parameters.loop_list.things_to_load.response.level = 'mouse';

% Velocities as explanatory/independent variable.
parameters.loop_list.things_to_load.explanatory.dir = {[parameters.dir_exper 'behavior\spontaneous\concatenated velocity by behavior\'], 'mouse', '\'};
parameters.loop_list.things_to_load.explanatory.filename= {'average_velocity_byinstance_', period, '.mat'};
parameters.loop_list.things_to_load.explanatory.variable= {'average_velocity'}; 
parameters.loop_list.things_to_load.explanatory.level = 'mouse';

% Output
parameters.loop_list.things_to_save.results.dir = {[parameters.dir_exper '\regression analysis\walk velocity\spontaneous\'], 'mouse', '\'};
parameters.loop_list.things_to_save.results.filename= {'regression_results_', period, '.mat'};
parameters.loop_list.things_to_save.results.variable= {'regression_results{', 'index_iterator',', 1}'}; 
parameters.loop_list.things_to_save.results.level = 'mouse';

RunAnalysis({@RegressData}, parameters);
