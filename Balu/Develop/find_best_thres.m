function best_thres = find_best_thres( precision, recall, thres, criterion, margin )

switch lower(criterion)
    case 'n-p'
        % find the best threshold using Neyman-Person strategy
        % margin must be in [0,1]
        lambda = 1e6;
        J = (lambda * abs( recall - margin ) ) + (1 - precision);
        best_thres = thres( find( J == min(J), 1, 'first') );
        
    case 'f-score'
        % find the best threshold using F-score measure
        J = 2 * (recall .* precision) ./ (eps + recall + precision);
        best_thres = thres( find( J == max(J), 1, 'first') );
end

end
