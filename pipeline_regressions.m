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
indices = find(tril(ones(number_of_sources)));

period = 'walk';

% Always clear loop list first. 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterations.
parameters.loop_list.iterators = {'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
              'index', {'loop_list.indices'}, 'index_iterator'};

parameters.loop_variables.mice_all = parameters.mice_all;
parameters.loop_variables.indices = indices;

% Number of bins to use in the histogram output.
parameters.histogram_nBins

% Things to load
% Correlations 
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'fluorescence analysis\correlations\spontaneous\instances\'], 'mouse', '\'};
parameters.loop_list.things_to_load.data.filename= {'correlations_', 'period', '.mat'};
parameters.loop_list.things_to_load.data.variable= {'correlations'}; 
parameters.loop_list.things_to_load.data.level = 'mouse';

% Things to save
% histogram values
parameters.loop_list.things_to_save.histogram.dir = {[parameters.dir_exper 'regression analysis\walk velocity\response distributions\spontaneous'], 'mouse', '\'};
parameters.loop_list.things_to_save.histogram.filename= {'correlations_response_distribution_', 'period', '.mat'};
parameters.loop_list.things_to_save.histogram.variable= {'distribution{', 'index_iterator', '}'}; 
parameters.loop_list.things_to_save.histogram.level = 'mouse';


%% 

distributions = NaN(size(correlations_reshaped,2), 36); % 36 is the  default number of bins 
 a = squeeze(correlations(1,10,:));
[vels, places] = sort(average_velocity);
a = a(places);



