function best_thres = train_classifier_svm(features, train_targets, num_folders, sf, svm_opts)

% include toolboxes
run('../cp_toolbox/cp_setup.m');
run('../../thesis_toolboxes/ToolboxBalu3/balu_setup.m');
addpath('../../thesis_toolboxes/libsvm-3.1/matlab');

% feature selection
train_features = features * sf.X_load;

% build arguments
tr_params = get_svm_parameters( svm_opts );

% tr_params = '-s 0 -t 1 -d 2 -r 1 -b 1';
% model = svmtrain(Y_tr, X_tr, tr_params);
te_params =  '-b 0';

% we build and train the SVM classifier
indices = crossvalind('KFold', train_targets, num_folders );
input_data = train_features;
input_targets = train_targets;

labels_cell = cell(1,num_folders);
scores_cell = cell(1,num_folders);

for fold = 1:num_folders
    fprintf('\nTrainning Fold %d/%d ...\n', fold, num_folders);
    
    test = indices == fold;
    train = ~test;

    % libsvm train
    tic
    
    svm_model = svmtrain( input_targets(train,1), input_data(train,:), tr_params);
    
    tElapsed = toc;
    cp_toc2time(tElapsed);
    
    fprintf( 'done.\n' );
    
    % evaluate model in test set
    fprintf( '--- evaluate model in test set ... ' );
    [classes, acc_test, prob_test] = svmpredict( input_targets(test,:), ...
        input_data(test,:), ...
        svm_model, ...
        te_params);
    
    fprintf('done.');
  
    % save performance in testing set
    labels_cell{fold} = input_targets(test,:);

    % save temporally the model
    scores_cell{fold} = prob_test;

end


% --- EVALUATE TEST FOLDS

scores = cell2mat( scores_cell' );
targets = cell2mat( labels_cell' );

% compute precision-recall curve
fprintf('--- computing precision-recall curves...\n');
[recall, precision, thres, AUC] = perfcurve( ...
    targets, ...
    scores, ...
    1, ...
    'xCrit', 'reca', ...
    'yCrit', 'prec');

figure, plot( recall, precision, 'r'),
title('Comparison presion-recall curves')

figure, plot( thres, recall, 'r'),
title('threshold vs recall')

figure, plot( thres, precision, 'r'),
title('threshold vs precision')

% find best threshold
best_thres = find_best_thres( precision, recall, thres, 'n-p', 95);
        
% --- EVALUATION IN VALIDATION DATA ---
%     fprintf('evaluation in testing data:\n'); 
    
%     % evaluation over real test dataset
%     test_targets(test_targets == -1) = 0;
%     cp_validation = classperf( test_targets);
   
%     F_score = (2 * PRC * TPR) / (eps + PRC + TPR);

%     predicted_classes = ones( numel(st), 1);
%     predicted_classes( st < best_thres ) = 0;
%     classperf( cp_validation, predicted_classes);
%     
%     % F-score evaluation
%     TP = cp_validation.DiagnosticTable(1);
%     FP = cp_validation.DiagnosticTable(3);
%     FN = cp_validation.DiagnosticTable(2);
%     TN = cp_validation.DiagnosticTable(4);
%     
%     TPR = TP / (TP + FN);   %true positive rate or recall
%     PRC = TP / (TP + FP);   % precision or positive predicted value
%     F_score = (2 * PRC * TPR) / (eps + PRC + TPR);
%     
%     fprintf('Best Threshold: %.4f, Pf: %.5f\tSn: %.5f\t1-Sp: %.5f\t F-measure: %.5f\n', ...
%         best_thres, cp_validation.CorrectRate, cp_validation.Sensitivity, ...
%         1-cp_validation.Specificity, F_score);

%     fprintf('Best Threshold: %.4f, Pf: %.5f\tSn: %.5f\t1-Sp: %.5f\t F-measure: %.5f\n', ...
%         best_thres, cp_validation.CorrectRate, cp_validation.Sensitivity, ...
%         1-cp_validation.Specificity, F_score);



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

