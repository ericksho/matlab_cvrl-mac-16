%seteamos vlfeat
run('vlfeat-0.9.19/toolbox/vl_setup')
%seteamos spams
run('spams-matlab/start_spams')
%agregamos balu al path
addpath(genpath('Balu'));
%agregamos bsif al path
addpath(genpath('bsif_code_and_data'));