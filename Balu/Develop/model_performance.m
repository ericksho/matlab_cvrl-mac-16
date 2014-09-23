function svm_model = model_performance( train_data, test_data, sf, best_thres, svm_opts)
% select the best classifier using performance 's F-score in the testing
% set

% include toolboxes
run('../cp_toolbox/cp_setup.m');
addpath('../../thesis_toolboxes/libsvm-3.1/matlab');

% load data
train_features = train_data.features * sf.X_load;
train_targets = train_data.targets;

test_features = test_data.features * sf.X_load;
test_targets = test_data.targets;

% train model using all train instances
tr_params = get_svm_parameters( svm_opts );
% te_params =  '-b 0';

% libsvm train
fprintf('\nTrainning ... \n');
tic

model = svmtrain( train_targets, train_features, tr_params);

tElapsed = toc;
cp_toc2time(tElapsed);
fprintf( 'done.\n' );

fprintf('\n--- Evaluate in train set ---\n');
[~, ~, train_scores] = svmpredict( ...
    train_targets, ...
    train_features, ...
    model );

fprintf('\n--- Evaluate in test set ---\n');
[~, ~, test_scores] = svmpredict( ...
    test_targets, ...
    test_features, ...
    model );


fprintf('\n--- Computing precision-recall curves in train set ---\n');
[train_recall, train_precision, train_thres, train_AUC] = perfcurve( ...
    train_targets, ...
    train_scores, ...
    1, ...
    'xCrit', 'reca', ...
    'yCrit', 'prec');

fprintf('\n--- Computing precision-recall curves in test set ---\n');
[test_recall, test_precision, test_thres, test_AUC] = perfcurve( ...
    test_targets, ...
    test_scores, ...
    1, ...
    'xCrit', 'reca', ...
    'yCrit', 'prec');

% statistics
f_score = 2 * (test_recall .* test_precision) ./ (eps + test_recall + test_precision);
mean_recall = mean( test_recall( ~isnan(test_recall) ) );
mean_precision = mean( test_precision( ~isnan(test_precision) ) );
mean_AUC = mean( test_AUC( ~isnan(test_AUC) ) );
mean_fscore = mean( f_score( ~isnan( f_score ) ) );

fprintf('Best Threshold: %.4f, avg_recall: %.4f\tavg_precision: %.4f\tavg_AUC: %.4f\t avg_F-measure: %.4f\n', ...
        best_thres, mean_recall, mean_precision, mean_AUC, mean_fscore);
    
figure, plot( train_recall, train_precision, 'r'),
hold on,
plot( test_recall, test_precision, 'b'),
hold off
title('Comparison presion-recall curves')

figure, plot( test_thres, test_recall, 'r'),
title('threshold vs recall')

figure, plot( test_thres, test_precision, 'r'),
title('threshold vs precision')


% build model structure
svm_model.name = 'svm';
svm_model.model = model;
svm_model.thres = best_thres;
svm_model.opts = svm_opts;
end

function parameters = get_svm_parameters( svm_opts )

% build string for svm parameters
str_param = '';
switch svm_opts.type
    case 'c-svc'
        str_param = '-s 0';
    case 'nu-svc'
        str_param = '-s 1';
    case 'epsilon-SVR'
        str_param = '-s 3';
end

switch svm_opts.kernel
    case 'linear'
        str_param = sprintf('%s -t 0', str_param);
    case 'polynomial'
        str_param = sprintf('%s -t 1', str_param);
    case 'rbf'
        str_param = sprintf('%s -t 2', str_param);
end

switch svm_opts.prob_est
    case 0
        str_param = sprintf('%s -b 0', str_param);
    case 1
        str_param = sprintf('%s -b 1', str_param);
end

parameters = str_param;
end
