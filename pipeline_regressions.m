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

load([parameters.dir_exper 'periods_nametable.mat']);
periods_motorized = periods;

load([parameters.dir_exper 'periods_nametable_spontaneous.mat']);
periods_spontaneous = periods;

% Create a shared motorized & spontaneous table.
periods_bothConditions = [periods_motorized; periods_spontaneous]; 

% Make list of transformation types for iterating later.
transformations = {'not transformed'; 'Fisher transformed'};

% Number of sources (spatial sources across mice)
number_of_sources = 32; 

% Number of principal components calculated on correlations
number_of_PCs = 100;

% Lower triangle-only for correlations.
indices = find(tril(ones(number_of_sources), -1));

% Correlations vs scores 
parameters.loop_variables.data_type = {'correlations', 'PCA scores individual mouse'};
value_indices{1} = 1:((number_of_sources^2 - number_of_sources)/2);
value_indices{2} = 1:number_of_PCs; 
parameters.loop_variables.value_indices = value_indices;

% Loop variables.
parameters.loop_variables.mice_all = parameters.mice_all;
parameters.loop_variables.indices = indices;
parameters.loop_variables.transformations = transformations;
parameters.loop_variables.conditions = {'motorized'; 'spontaneous'};
parameters.loop_variables.result_types = {'betas'; 'r2s'};


%% Spontaneous-- Look at histograms of response variable by itself 
% (each correlation, not dependent on velocity yet) Is important for
% picking the style of GLM.
% 
% period = 'walk';
% 
% % Always clear loop list first. 
% if isfield(parameters, 'loop_list')
% parameters = rmfield(parameters,'loop_list');
% end
% 
% % Iterations.
% parameters.loop_list.iterators = {'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
%               'index', {'loop_variables.indices'}, 'index_iterator'};
% 
% % Number of bins to use in the histogram output.
% parameters.histogram_nBins = 200;
% parameters.histogram_value_range = [-1 1];
% 
% % Things to load
% % Correlations 
% parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'fluorescence analysis\correlations\spontaneous\'] 'mouse', '\instances reshaped\'};
% parameters.loop_list.things_to_load.data.filename= {'correlations_', period, '.mat'};
% parameters.loop_list.things_to_load.data.variable= {'correlations_reshaped(', 'index', ', :)'}; 
% parameters.loop_list.things_to_load.data.level = 'mouse';
% 
% % Things to save
% % histogram values
% parameters.loop_list.things_to_save.histogram.dir = {[parameters.dir_exper 'regression analysis\walk velocity\response distributions\spontaneous\'], 'mouse', '\'};
% parameters.loop_list.things_to_save.histogram.filename= {'correlations_response_distribution_histograms_', period, '.mat'};
% parameters.loop_list.things_to_save.histogram.variable= {'histograms{', 'index_iterator', ',1}'}; 
% parameters.loop_list.things_to_save.histogram.level = 'mouse';
% 
% RunAnalysis({@Histogram}, parameters);

%% 

distributions = NaN(size(correlations_reshaped,2), 36); % 36 is the  default number of bins 
 a = squeeze(correlations(1,10,:));
[vels, places] = sort(average_velocity);
a = a(places);

%% Regress all spontaneous
period = 'walk';
period_index = '190';

% Always clear loop list first. 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterations.
parameters.loop_list.iterators = {
              'data_type', {'loop_variables.data_type'}, 'data_type_iterator';
              'transformation', {'loop_variables.transformations'}, 'transformation_iterator';
              'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
              'index', {'loop_variables.value_indices{', 'data_type_iterator', '}'}, 'index_iterator'};

% Dimension of different predictors (so you can orient variables correctly
% before regressing)
parameters.predictorsDim = 1; 

% Squeeze data first
parameters.evaluation_instructions = {'data_evaluated = squeeze(parameters.response);'};

% Input 
% Correlations as response/ dependent variable.
parameters.loop_list.things_to_load.response.dir = {[parameters.dir_exper 'fluorescence analysis\'], 'data_type', '\', 'transformation', '\', 'mouse', '\instances reshaped\'};
parameters.loop_list.things_to_load.response.filename= {'values.mat'};
parameters.loop_list.things_to_load.response.variable= {['values{' period_index '}('], 'index_iterator' ,', 1,:)'}; 
parameters.loop_list.things_to_load.response.level = 'mouse';

% Velocities as explanatory/independent variable.
parameters.loop_list.things_to_load.explanatory.dir = {[parameters.dir_exper 'behavior\spontaneous\concatenated velocity by behavior\'], 'mouse', '\'};
parameters.loop_list.things_to_load.explanatory.filename= {'average_velocity_byinstance_', period, '.mat'};
parameters.loop_list.things_to_load.explanatory.variable= {'average_velocity'}; 
parameters.loop_list.things_to_load.explanatory.level = 'mouse';

% Output
parameters.loop_list.things_to_save.results_r2.dir = {[parameters.dir_exper '\regression analysis\walk velocity\'],  'data_type', '\', 'transformation',  '\spontaneous\results\', 'mouse', '\'};
parameters.loop_list.things_to_save.results_r2.filename= {'regression_results_r2s.mat'};
parameters.loop_list.things_to_save.results_r2.variable= {'r2s(', 'index_iterator', ', 1)'}; 
parameters.loop_list.things_to_save.results_r2.level = 'mouse';

parameters.loop_list.things_to_save.results_betas.dir = {[parameters.dir_exper '\regression analysis\walk velocity\'], 'data_type', '\','transformation','\spontaneous\results\', 'mouse', '\'};
parameters.loop_list.things_to_save.results_betas.filename= {'regression_results_betas.mat'};
parameters.loop_list.things_to_save.results_betas.variable= {'betas(', 'index_iterator', ', :)'}; 
parameters.loop_list.things_to_save.results_betas.level = 'mouse';

parameters.loop_list.things_to_rename = {{'data_evaluated', 'response'}};

RunAnalysis({@EvaluateOnData, @RegressData}, parameters);


%% Make vector of all motorized velocities 
condition = 'c_walk';
periods = {'176'; '177'; '178'; '179'};

% Always clear loop list first. 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterations.
parameters.loop_list.iterators = {
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
              'period', {'loop_variables.periods(:)'}, 'period_iterator'};

parameters.loop_variables.periods = periods;

% instructions for replicating.
parameters.toReplicate = {'parameters.period_velocities(', 'period_iterator', ')'};
parameters.replicateDims = {'[size(parameters.data, 3), 1]'};

% Velocities (cm/s) for each motor speed.
period_velocities = [2.25, 2.77, 3.33, 3.88]; 
parameters.period_velocities = period_velocities;

% Concatenation dimension
parameters.concatDim = 1;

% Inputs (reshaped correlations, for the number)
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'fluorescence analysis\correlations\not transformed\'], 'mouse', '\instances reshaped\'};
parameters.loop_list.things_to_load.data.filename= {'values.mat'};
parameters.loop_list.things_to_load.data.variable= {'values{', 'period', '}'}; 
parameters.loop_list.things_to_load.data.level = 'mouse';

% Output
parameters.loop_list.things_to_save.concatenated_data.dir = {[parameters.dir_exper 'regression analysis\walk velocity\velocity vectors\motorized\'], 'mouse','\'};
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
parameters.loop_variables.data_type = {'PCA scores individual mouse'};
% Iterations.
parameters.loop_list.iterators = {
               'data_type', {'loop_variables.data_type'}, 'data_type_iterator';
               'transformation', {'loop_variables.transformations'}, 'transformation_iterator';
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
              'period', {'loop_variables.periods(:)'}, 'period_iterator'};
parameters.loop_variables.periods = periods;

% Concatenation dimension
parameters.concatDim = 3;

% Inputs (reshaped correlations)
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'fluorescence analysis\'], 'data_type', '\', 'transformation', '\' 'mouse', '\instances reshaped\'};
parameters.loop_list.things_to_load.data.filename= {'values.mat'};
parameters.loop_list.things_to_load.data.variable= {'values{', 'period', '}'}; 
parameters.loop_list.things_to_load.data.level = 'mouse';

% Output
parameters.loop_list.things_to_save.concatenated_data.dir = {[parameters.dir_exper 'fluorescence analysis\'], 'data_type', '\', 'transformation', '\', 'mouse', '\all concatenated\'};
parameters.loop_list.things_to_save.concatenated_data.filename= {'values_all_c_walk.mat'};
parameters.loop_list.things_to_save.concatenated_data.variable= {'values_all'}; 
parameters.loop_list.things_to_save.concatenated_data.level = 'mouse';

RunAnalysis({@ConcatenateData}, parameters);

%% Motorized -- regress

condition = 'c_walk';

% Always clear loop list first. 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterations.
parameters.loop_list.iterators = {
              'data_type', {'loop_variables.data_type'}, 'data_type_iterator';
              'transformation', {'loop_variables.transformations'}, 'transformation_iterator';
              'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
              'index', {'loop_variables.value_indices{', 'data_type_iterator', '}'}, 'index_iterator'};

% Dimension of different predictors (so you can orient variables correctly
% before regressing)
parameters.predictorsDim = 2; 

% Input 
% Correlations as response/ dependent variable.
parameters.loop_list.things_to_load.response.dir = {[parameters.dir_exper 'fluorescence analysis\'], 'data_type', '\','transformation', '\', 'mouse', '\all concatenated\'};
parameters.loop_list.things_to_load.response.filename= {'values_all_c_walk.mat'};
parameters.loop_list.things_to_load.response.variable= {'values_all(', 'index', ',:)'}; 
parameters.loop_list.things_to_load.response.level = 'mouse';

% Velocities as explanatory/independent variable.
parameters.loop_list.things_to_load.explanatory.dir = {[parameters.dir_exper 'regression analysis\walk velocity\velocity vectors\motorized\'], 'mouse', '\'};
parameters.loop_list.things_to_load.explanatory.filename= {'velocity_vector.mat'};
parameters.loop_list.things_to_load.explanatory.variable= {'velocity_vector'}; 
parameters.loop_list.things_to_load.explanatory.level = 'mouse';

% Output
parameters.loop_list.things_to_save.results_r2.dir = {[parameters.dir_exper '\regression analysis\walk velocity\'], 'data_type', '\', 'transformation', '\motorized\results\', 'mouse', '\'};
parameters.loop_list.things_to_save.results_r2.filename= {'regression_results_r2s.mat'};
parameters.loop_list.things_to_save.results_r2.variable= {'r2s(', 'index_iterator', ', 1)'}; 
parameters.loop_list.things_to_save.results_r2.level = 'mouse';

parameters.loop_list.things_to_save.results_betas.dir = {[parameters.dir_exper '\regression analysis\walk velocity\'], 'data_type', '\', 'transformation', '\motorized\results\', 'mouse', '\'};
parameters.loop_list.things_to_save.results_betas.filename= {'regression_results_betas.mat'};
parameters.loop_list.things_to_save.results_betas.variable= {'betas(', 'index_iterator', ', :)'}; 
parameters.loop_list.things_to_save.results_betas.level = 'mouse';

RunAnalysis({@RegressData}, parameters);

%% Reshape regression results (Correlations only)
% (you end up using them all in reshaped form later anyway)

% Always clear loop list first. 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {
              'transformation', {'loop_variables.transformations'}, 'transformation_iterator';
              'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
              'condition', {'loop_variables.conditions(:)'}, 'condition_iterator';
              'result_type', {'loop_variables.result_types(:)'}, 'result_type_iterator'};
     
parameters.pixels = [number_of_sources, number_of_sources];
parameters.indices_of_mask = indices;
parameters.pixelDim = 2; 

parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper '\regression analysis\walk velocity\correlations\'],'transformation', '\', 'condition', '\results\', 'mouse', '\'};
parameters.loop_list.things_to_load.data.filename= {'regression_results_', 'result_type', '.mat'};
parameters.loop_list.things_to_load.data.variable= {'result_type'}; 
parameters.loop_list.things_to_load.data.level = 'result_type';

parameters.loop_list.things_to_save.data_matrix_filled.dir = {[parameters.dir_exper '\regression analysis\walk velocity\correlations\'], 'transformation', '\', 'condition', '\results\', 'mouse', '\'};
parameters.loop_list.things_to_save.data_matrix_filled.filename= {'regression_results_reshaped_', 'result_type', '.mat'};
parameters.loop_list.things_to_save.data_matrix_filled.variable= {'result_type'}; 
parameters.loop_list.things_to_save.data_matrix_filled.level = 'result_type';

RunAnalysis({@FillMasks_forRunAnalysis}, parameters);

%% Compare differences in betas between conditions
% Load, subtract, save 

% Always clear loop list first. 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {
                    'data_type', {'loop_variables.data_type'}, 'data_type_iterator';
                    'transformation', {'loop_variables.transformations'}, 'transformation_iterator';
                    'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'};
                                  
parameters.pixels = [number_of_sources, number_of_sources];
parameters.indices_of_mask = indices;

% Motorized
parameters.loop_list.things_to_load.subtract_this.dir = {[parameters.dir_exper '\regression analysis\walk velocity\'],'data_type', '\', 'transformation', '\spontaneous\results\', 'mouse', '\'};
parameters.loop_list.things_to_load.subtract_this.filename= {'regression_results_betas.mat'};
parameters.loop_list.things_to_load.subtract_this.variable= {'betas'}; 
parameters.loop_list.things_to_load.subtract_this.level = 'mouse';

% Spontaneous 
parameters.loop_list.things_to_load.subtract_from_this.dir = {[parameters.dir_exper 'regression analysis\walk velocity\'], 'data_type', '\', 'transformation', '\motorized\results\', 'mouse', '\'};
parameters.loop_list.things_to_load.subtract_from_this.filename= {'regression_results_betas.mat'};
parameters.loop_list.things_to_load.subtract_from_this.variable= {'betas'}; 
parameters.loop_list.things_to_load.subtract_from_this.level = 'mouse';

% Output
parameters.loop_list.things_to_save.data_subtracted.dir = {[parameters.dir_exper '\regression analysis\walk velocity\'], 'data_type', '\', 'transformation', '\difference from spontaneous to motorized\', 'mouse', '\'};
parameters.loop_list.things_to_save.data_subtracted.filename= {'betas_difference.mat'};
parameters.loop_list.things_to_save.data_subtracted.variable= {'betas_difference'}; 
parameters.loop_list.things_to_save.data_subtracted.level = 'mouse';

%parameters.loop_list.things_to_rename = {{'data_subtracted', 'data'}};

RunAnalysis({@SubtractData}, parameters);

%% Concatenate all spontaneous & motorized velocities

% Always clear loop list first. 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterations.
parameters.loop_list.iterators = {'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
              'condition', {'loop_variables.conditions(:)'}, 'condition_iterator'};

parameters.concatDim = 1; 

% Input
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'regression analysis\walk velocity\'], '\velocity vectors\', 'condition', '\', 'mouse','\'};
parameters.loop_list.things_to_load.data.filename= {'velocity_vector.mat'};
parameters.loop_list.things_to_load.data.variable= {'velocity_vector'}; 
parameters.loop_list.things_to_load.data.level = 'condition';

% Output
parameters.loop_list.things_to_save.concatenated_data.dir = {[parameters.dir_exper 'regression analysis\walk velocity\velocity vectors\motorized & spontaneous together\'], 'mouse','\'};
parameters.loop_list.things_to_save.concatenated_data.filename= {'velocity_vector.mat'};
parameters.loop_list.things_to_save.concatenated_data.variable= {'velocity_vector'}; 
parameters.loop_list.things_to_save.concatenated_data.level = 'mouse';

RunAnalysis({@ConcatenateData}, parameters);

%% Concatenate all spontaneous & motorized corr together
% (Use a trick with the way Concatenate data works to load in from 2
% separate locations, works only when concatenating 2 things): load spontaneous in
% as data (needs to be second to match with above), motorized as
% concatenated_data.

% Always clear loop list first. 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

parameters.loop_variables.periods = {'176'; '177'; '178'; '179'; '190'};

% Iterations.
parameters.loop_list.iterators = {
    'data_type', {'loop_variables.data_type'}, 'data_type_iterator';
    'transformation', {'loop_variables.transformations'}, 'transformation_iterator';
    'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator';
    'period', {'loop_variables.periods'}, 'period_iterator';
     };

parameters.concatDim = 3; 

% Input
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'fluorescence analysis\'], 'data_type', '\', 'transformation', '\', 'mouse', '\instances reshaped\'};
parameters.loop_list.things_to_load.data.filename= {'values.mat'};
parameters.loop_list.things_to_load.data.variable= {'values{', 'period', '}'}; 
parameters.loop_list.things_to_load.data.level = 'mouse';

% Output
parameters.loop_list.things_to_save.concatenated_data.dir = {[parameters.dir_exper 'regression analysis\walk velocity\'],  'data_type', '\', 'transformation', '\motorized & spontaneous together\values together\', 'mouse','\'};
parameters.loop_list.things_to_save.concatenated_data.filename= {'values_all_walk.mat'};
parameters.loop_list.things_to_save.concatenated_data.variable= {'values_all_walk'}; 
parameters.loop_list.things_to_save.concatenated_data.level = 'mouse';

RunAnalysis({@ConcatenateData}, parameters);

%% Regress spontaneous & motorized together.

% Always clear loop list first. 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterations.
parameters.loop_list.iterators = {
              'data_type', {'loop_variables.data_type'}, 'data_type_iterator';
              'transformation', {'loop_variables.transformations'}, 'transformation_iterator';
              'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
              'index', {'loop_variables.value_indices{', 'data_type_iterator', '}'}, 'index_iterator'};

% Squeeze data first
parameters.evaluation_instructions = {'data_evaluated = squeeze(parameters.response);'};

% Dimension of different predictors (so you can orient variables correctly
% before regressing)
parameters.predictorsDim = 2; 

% Input 
% (correlations as response/dependent variable)
parameters.loop_list.things_to_load.response.dir = {[parameters.dir_exper 'regression analysis\walk velocity\'], 'data_type', '\', 'transformation', '\motorized & spontaneous together\values together\', 'mouse','\'};
parameters.loop_list.things_to_load.response.filename= {'values_all_walk.mat'};
parameters.loop_list.things_to_load.response.variable= {'values_all_walk(', 'index', ',:)'}; 
parameters.loop_list.things_to_load.response.level = 'mouse';

% (Velocities as explanatory/independent variable.)
parameters.loop_list.things_to_load.explanatory.dir = {[parameters.dir_exper 'regression analysis\walk velocity\velocity vectors\motorized & spontaneous together\'], 'mouse', '\'};
parameters.loop_list.things_to_load.explanatory.filename= {'velocity_vector.mat'};
parameters.loop_list.things_to_load.explanatory.variable= {'velocity_vector'}; 
parameters.loop_list.things_to_load.explanatory.level = 'mouse';

% Output
parameters.loop_list.things_to_save.results_r2.dir = {[parameters.dir_exper '\regression analysis\walk velocity\'], 'data_type', '\',  'transformation', '\motorized & spontaneous together\results\', 'mouse', '\'};
parameters.loop_list.things_to_save.results_r2.filename= {'regression_results_r2s.mat'};
parameters.loop_list.things_to_save.results_r2.variable= {'r2s(', 'index_iterator', ', 1)'}; 
parameters.loop_list.things_to_save.results_r2.level = 'mouse';

parameters.loop_list.things_to_save.results_betas.dir = {[parameters.dir_exper '\regression analysis\walk velocity\'], 'data_type', '\',  'transformation', '\motorized & spontaneous together\results\', 'mouse', '\'};
parameters.loop_list.things_to_save.results_betas.filename= {'regression_results_betas.mat'};
parameters.loop_list.things_to_save.results_betas.variable= {'betas(', 'index_iterator', ', :)'}; 
parameters.loop_list.things_to_save.results_betas.level = 'mouse';

parameters.loop_list.things_to_rename = {{'data_evaluated', 'response'}};

RunAnalysis({@EvaluateOnData, @RegressData}, parameters);

%% Reshape motorized & spontaneous regression results 
% Always clear loop list first. 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {'transformation', {'loop_variables.transformations'}, 'transformation_iterator';
                                  'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator';
                                  'result_type', {'loop_variables.result_types(:)'}, 'result_type_iterator'};
parameters.pixels = [number_of_sources, number_of_sources];
parameters.indices_of_mask = indices;
parameters.pixelDim = 2; 

parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'regression analysis\walk velocity\'], 'transformation', '\motorized & spontaneous together\results\', 'mouse', '\'};
parameters.loop_list.things_to_load.data.filename= {'regression_results_', 'result_type', '.mat'};
parameters.loop_list.things_to_load.data.variable= {'result_type'}; 
parameters.loop_list.things_to_load.data.level = 'result_type';

parameters.loop_list.things_to_save.data_matrix_filled.dir = {[parameters.dir_exper 'regression analysis\walk velocity\'], 'transformation', '\motorized & spontaneous together\results\', 'mouse', '\'};
parameters.loop_list.things_to_save.data_matrix_filled.filename= {'regression_results_reshaped_', 'result_type', '.mat'};
parameters.loop_list.things_to_save.data_matrix_filled.variable= {'result_type'}; 
parameters.loop_list.things_to_save.data_matrix_filled.level = 'result_type';

RunAnalysis({@FillMasks_forRunAnalysis}, parameters);

%% Plot correlation regression results
cmap_diffs = flipud(cbrewer('div', 'RdBu', 256, 'nearest'));
mouse = '1087'; 

for i = 1:numel(transformations)
    transformation = transformations{i};
    figure; 
    
    % plot spontaneous;
    load(['Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\regression analysis\walk velocity\' transformation '\spontaneous\results\1087\regression_results_reshaped_betas.mat']); 
    subplot(3, 3, 1); imagesc(betas(:,:,1)); axis square;
    title('b1 spontaneous'); colorbar; colormap(gca, cmap_diffs);
    caxis([-0.1 0.1]);
    
    subplot(3, 3, 2); imagesc(betas(:,:,2)); axis square;
    title('b0 spontaneous'); colorbar;
    caxis([0 2]);
    
    load(['Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\regression analysis\walk velocity\' transformation '\spontaneous\results\1087\regression_results_reshaped_r2s.mat']); 
    subplot(3, 3,3); imagesc(r2s); axis square;
    title('r2s spontaneous'); colorbar; 
    caxis([0 0.05]);
    
    % plot motorized
    load(['Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\regression analysis\walk velocity\' transformation '\motorized\results\1087\regression_results_reshaped_betas.mat']); 
    subplot(3, 3,4); imagesc(betas(:,:,1)); axis square;
    title('b1 motorized'); colorbar; colormap(gca, cmap_diffs);
    caxis([-0.1 0.1]);
    
    subplot(3, 3, 5); imagesc(betas(:,:,2)); axis square;
    title('b0 motorized'); colorbar;
    caxis([0 2]);
    
    load(['Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\regression analysis\walk velocity\' transformation '\motorized\results\1087\regression_results_reshaped_r2s.mat']); 
    subplot(3, 3, 6); imagesc(r2s); axis square;
    title('r2s motorized'); colorbar;
    caxis([0 0.05]);
    
    % % plot differences
    % load('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\regression analysis\walk velocity\difference from spontaneous to motorized\1087\betas_difference.mat'); 
    % subplot(4, 3, 7); imagesc(betas_difference(:,:,1)); axis square;
    % title('b1 difference'); colorbar; colormap(gca, cmap_diffs); 
    % caxis([-0.09 0.09]); 
    % 
    % subplot(4, 3, 8); imagesc(betas_difference(:,:,2)); axis square;
    % title('b0 difference'); colorbar; colormap(gca, cmap_diffs);
    % caxis([-0.5 0.5]);
    
    % results from motorized & spontaneous together.
    load(['Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\regression analysis\walk velocity\' transformation '\motorized & spontaneous together\results\1087\regression_results_reshaped_betas.mat']); 
    subplot(3, 3, 7); imagesc(betas(:,:,1)); axis square;
    title('b1 motorized & spontaneous'); colorbar; colormap(gca, cmap_diffs);
    caxis([-0.1 0.1]);
    
    subplot(3, 3, 8); imagesc(betas(:,:,2)); axis square;
    title('b0 motorized & spontaneous'); colorbar;
    caxis([0 2]);
    
    load(['Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\regression analysis\walk velocity\' transformation '\motorized & spontaneous together\results\1087\regression_results_reshaped_r2s.mat']); 
    subplot(3, 3, 9); imagesc(r2s); axis square;
    title('r2s motorized & spontaneous'); colorbar;
    caxis([0 0.05]);
    
    sgtitle([mouse ', ' transformation]);

end 

%% Plot results of PC score regressions

mouse = '1087'; 
conditions = {'spontaneous', 'motorized', 'motorized & spontaneous together'};
for i = 1:numel(transformations)
    transformation = transformations{i};
    
    figure; 
    for condi = 1:numel(conditions)
        condition = conditions{condi};
        load(['Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\regression analysis\walk velocity\PCA scores individual mouse\' transformation '\' condition '\results\' mouse '\regression_results_betas.mat'])
        load(['Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\regression analysis\walk velocity\PCA scores individual mouse\' transformation '\' condition '\results\' mouse '\regression_results_r2s.mat'])    
        subplot(2,3,condi); 
        imagesc(betas(1:20, :) ); colorbar; caxis([-1 1])
        xticks([1:2]); 
        xticklabels({'b', 'intercept'});
        title(condition);

        subplot(2,3,condi + 3); 
        imagesc(r2s(1:20)); colorbar; caxis([0  0.05])
        xlabel('r2');
    
    end
    sgtitle([mouse ', ' transformation]);
end 

%% Try fitting exponential
% [I don't think I feel comfortable *really* trying this or another shape 
% until I get more mice analyzed to this point]

x = velocity_vector; 
y = all_correlations(indices(15),:)';
g = fittype('a-b*exp(-c*x)');
f0 = fit(x,y,g,'StartPoint',[[ones(size(x)), -exp(-x)]\y; 1]);

plotx = [0 : 0.1: 18];
ploty_main = f0.a -  f0.b * exp(-f0.c * plotx); 

intervals =  confint(f0);

ploty_low = intervals(1,1) -  intervals(1,2) * exp(- intervals(1,3) *plotx); 
ploty_high = intervals(2,1) -  intervals(2,2) * exp(- intervals(2,3) *plotx); 

figure; hold on;
plot(velocity_vector, y, 'ok')
plot(plotx, ploty_main, 'b');
plot(plotx, ploty_low, 'r');
plot(plotx, ploty_high, 'r');
%ylim([-1 1]);


%% Regress std deviation values? Would work only for motorized.
