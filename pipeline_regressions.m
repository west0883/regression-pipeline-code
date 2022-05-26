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
parameters.loop_list.things_to_save.results_r2.dir = {[parameters.dir_exper '\regression analysis\walk velocity\spontaneous\results\'], 'mouse', '\'};
parameters.loop_list.things_to_save.results_r2.filename= {'regression_results_r2s.mat'};
parameters.loop_list.things_to_save.results_r2.variable= {'r2s(', 'index_iterator', ', 1)'}; 
parameters.loop_list.things_to_save.results_r2.level = 'mouse';

parameters.loop_list.things_to_save.results_betas.dir = {[parameters.dir_exper '\regression analysis\walk velocity\spontaneous\results\'], 'mouse', '\'};
parameters.loop_list.things_to_save.results_betas.filename= {'regression_results_betas.mat'};
parameters.loop_list.things_to_save.results_betas.variable= {'betas(', 'index_iterator', ', :)'}; 
parameters.loop_list.things_to_save.results_betas.level = 'mouse';

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

%% Motorized -- Concatenate all continued walk correlations into one matrix
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

% Concatenation dimension
parameters.concatDim = 2;

% Inputs (reshaped correlations)
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'fluorescence analysis\correlations\motorized\'] 'mouse', '\instances reshaped\'};
parameters.loop_list.things_to_load.data.filename= {['correlations_' condition '_'], 'period', '.mat'};
parameters.loop_list.things_to_load.data.variable= {'correlations_reshaped'}; 
parameters.loop_list.things_to_load.data.level = 'period';

% Output
parameters.loop_list.things_to_save.concatenated_data.dir = {[parameters.dir_exper 'fluorescence analysis\correlations\motorized\'] 'mouse', '\all concatenated\'};
parameters.loop_list.things_to_save.concatenated_data.filename= {'correlations_all_c_walk.mat'};
parameters.loop_list.things_to_save.concatenated_data.variable= {'correlations'}; 
parameters.loop_list.things_to_save.concatenated_data.level = 'mouse';

RunAnalysis({@ConcatenateData}, parameters);

%% Motorized -- regress

number_of_sources = 29; 
indices = find(tril(ones(number_of_sources), -1));

condition = 'c_walk';

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
parameters.predictorsDim = 2; 

% Input 
% Correlations as response/ dependent variable.
parameters.loop_list.things_to_load.response.dir = {[parameters.dir_exper 'fluorescence analysis\correlations\motorized\'] 'mouse', '\all concatenated\'};
parameters.loop_list.things_to_load.response.filename= {'correlations_all_c_walk.mat'};
parameters.loop_list.things_to_load.response.variable= {'correlations(', 'index', ',:)'}; 
parameters.loop_list.things_to_load.response.level = 'mouse';

% Velocities as explanatory/independent variable.
parameters.loop_list.things_to_load.explanatory.dir = {[parameters.dir_exper 'regression analysis\walk velocity\motorized\velocity vectors\'], 'mouse', '\'};
parameters.loop_list.things_to_load.explanatory.filename= {'velocity_vector.mat'};
parameters.loop_list.things_to_load.explanatory.variable= {'velocity_vector'}; 
parameters.loop_list.things_to_load.explanatory.level = 'mouse';

% Output
parameters.loop_list.things_to_save.results_r2.dir = {[parameters.dir_exper '\regression analysis\walk velocity\motorized\results\'], 'mouse', '\'};
parameters.loop_list.things_to_save.results_r2.filename= {'regression_results_r2s.mat'};
parameters.loop_list.things_to_save.results_r2.variable= {'r2s(', 'index_iterator', ', 1)'}; 
parameters.loop_list.things_to_save.results_r2.level = 'mouse';

parameters.loop_list.things_to_save.results_betas.dir = {[parameters.dir_exper '\regression analysis\walk velocity\motorized\results\'], 'mouse', '\'};
parameters.loop_list.things_to_save.results_betas.filename= {'regression_results_betas.mat'};
parameters.loop_list.things_to_save.results_betas.variable= {'betas(', 'index_iterator', ', :)'}; 
parameters.loop_list.things_to_save.results_betas.level = 'mouse';

RunAnalysis({@RegressData}, parameters);

%% Reshape regression results 
% (you end up using them all in reshaped form later anyway)
number_of_sources = 29; 
indices = find(tril(ones(number_of_sources), -1));

parameters.loop_variables.conditions = {'motorized'; 'spontaneous'};
parameters.loop_variables.result_types = {'betas'; 'r2s'};
parameters.loop_variables.mice_all = parameters.mice_all;

% Always clear loop list first. 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator';
                                  'condition', {'loop_variables.conditions(:)'}, 'condition_iterator';
                                  'result_type', {'loop_variables.result_types(:)'}, 'result_type_iterator'};
parameters.pixels = [number_of_sources, number_of_sources];
parameters.indices_of_mask = indices;
parameters.pixelDim = 2; 

parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper '\regression analysis\walk velocity\'], 'condition', '\results\', 'mouse', '\'};
parameters.loop_list.things_to_load.data.filename= {'regression_results_', 'result_type', '.mat'};
parameters.loop_list.things_to_load.data.variable= {'result_type'}; 
parameters.loop_list.things_to_load.data.level = 'result_type';

parameters.loop_list.things_to_save.data_matrix_filled.dir = {[parameters.dir_exper '\regression analysis\walk velocity\'], 'condition', '\results\', 'mouse', '\'};
parameters.loop_list.things_to_save.data_matrix_filled.filename= {'regression_results_reshaped_', 'result_type', '.mat'};
parameters.loop_list.things_to_save.data_matrix_filled.variable= {'result_type'}; 
parameters.loop_list.things_to_save.data_matrix_filled.level = 'result_type';

RunAnalysis({@FillMasks_forRunAnalysis}, parameters);

%% Compare differences in betas between conditions
% Load, subtract, save
number_of_sources = 29; 
indices = find(tril(ones(number_of_sources), -1));

parameters.loop_variables.conditions = {'motorized'; 'spontaneous'};
parameters.loop_variables.mice_all = parameters.mice_all;

% Always clear loop list first. 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'};
                                  
parameters.pixels = [number_of_sources, number_of_sources];
parameters.indices_of_mask = indices;

% Motorized
parameters.loop_list.things_to_load.subtract_this.dir = {[parameters.dir_exper '\regression analysis\walk velocity\spontaneous\results\'], 'mouse', '\'};
parameters.loop_list.things_to_load.subtract_this.filename= {'regression_results_reshaped_betas.mat'};
parameters.loop_list.things_to_load.subtract_this.variable= {'betas'}; 
parameters.loop_list.things_to_load.subtract_this.level = 'mouse';

% Spontaneous 
parameters.loop_list.things_to_load.subtract_from_this.dir = {[parameters.dir_exper 'regression analysis\walk velocity\motorized\results\'], 'mouse', '\'};
parameters.loop_list.things_to_load.subtract_from_this.filename= {'regression_results_reshaped_betas.mat'};
parameters.loop_list.things_to_load.subtract_from_this.variable= {'betas'}; 
parameters.loop_list.things_to_load.subtract_from_this.level = 'mouse';

% Output
parameters.loop_list.things_to_save.data_subtracted.dir = {[parameters.dir_exper '\regression analysis\walk velocity\difference from spontaneous to motorized\'], 'mouse', '\'};
parameters.loop_list.things_to_save.data_subtracted.filename= {'betas_difference.mat'};
parameters.loop_list.things_to_save.data_subtracted.variable= {'betas_difference'}; 
parameters.loop_list.things_to_save.data_subtracted.level = 'mouse';

%parameters.loop_list.things_to_rename = {{'data_subtracted', 'data'}};

RunAnalysis({@SubtractData}, parameters);

%% Plot the difference in betas 
cmap_diffs = flipud(cbrewer('div', 'RdBu', 256, 'nearest'));
mouse = '1087'; 

figure; 

% plot spontaneous;
load('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\regression analysis\walk velocity\spontaneous\results\1087\regression_results_reshaped_betas.mat'); 
subplot(3,3,1); imagesc(betas(:,:,1)); axis square;
title('b1 spontaneous'); colorbar;
caxis([-0.08 0.01]);

subplot(3,3,2); imagesc(betas(:,:,2)); axis square;
title('b0 spontaneous'); colorbar;
caxis([0 1]);

load('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\regression analysis\walk velocity\spontaneous\results\1087\regression_results_reshaped_r2s.mat'); 
subplot(3,3,3); imagesc(r2s); axis square;
title('r2s spontaneous'); colorbar;
caxis([0 0.02]);

% plot motorized
load('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\regression analysis\walk velocity\motorized\results\1087\regression_results_reshaped_betas.mat'); 
subplot(3,3,4); imagesc(betas(:,:,1)); axis square;
title('b1 motorized'); colorbar;
caxis([-0.08 0.01]);

subplot(3,3,5); imagesc(betas(:,:,2)); axis square;
title('b0 motorized'); colorbar;
caxis([0 1]);

load('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\regression analysis\walk velocity\motorized\results\1087\regression_results_reshaped_r2s.mat'); 
subplot(3,3,6); imagesc(r2s); axis square;
title('r2s motorized'); colorbar;
caxis([0 0.02]);

% plot differences
load('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\regression analysis\walk velocity\difference from spontaneous to motorized\1087\betas_difference.mat'); 
subplot(3,3,7); imagesc(betas_difference(:,:,1)); axis square;
title('b1 difference'); colorbar; colormap(gca, cmap_diffs); 
caxis([-0.09 0.09]); 

subplot(3,3,8); imagesc(betas_difference(:,:,2)); axis square;
title('b0 difference'); colorbar; colormap(gca, cmap_diffs);
caxis([-0.5 0.5]);

sgtitle([mouse ', diff spontaneous to motorized']);

%% Concatenate all spontaneous & motorized velocities together

%% Concatenate all spontaneous & motorized corr together