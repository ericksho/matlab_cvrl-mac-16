function out_struct = feature_selection_plsr(varargin)
%feature_selection( method , varargin )
%   this function select a feature subset from a set of predictive features
%
%   inputs:
%   

% manage input parameters
input_features = varargin{1};
input_labels = varargin{2};

if iscell( input_features )
    features = cell2mat( cellfun(@(x) x', input_features', 'Uni', 0) );
    labels = cell2mat( cellfun(@(x) x', input_labels', 'Uni', 0) );
    labels( labels == -1 ) = 0;
    
    num_components = size( features, 2);
    if numel( varargin ) == 3
        num_components = varargin{3};
    end
    
    [XL,YL,XS,YS,dummy,PCTVAR] = plsregress( features , labels , num_components);
    
else
    num_components = size( input_features, 2);
    if numel( varargin ) == 3
        num_components = varargin{3};
    end
    
    [XL,YL,XS,YS,dummy,PCTVAR] = plsregress( input_features , input_labels , num_components);
    
end

out_struct.X_load = XL;
out_struct.Y_load = YL;
out_struct.X_scored = XS;
out_struct.Y_scored = YS;
out_struct.PCTVAR = PCTVAR;
out_struct.num_components = num_components;

% figure
figure, plot(1:num_components , cumsum( 100 * PCTVAR(2,:) ) , '-b');
xlabel('Number of PLS components');
ylabel('Percent Variance Explained in y');

end
