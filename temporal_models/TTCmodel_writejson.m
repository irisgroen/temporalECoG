function opts_json = TTCmodel_writejson()
% Writes out a json file with parameters for model fitting.
%
% TTC model parameters: 
%          1. weight - relative weight on transient channel ([0 1]) 
%          2. shift -- time between stimulus onset and when the signal reaches
%               the cortex (seconds)
%          3. scale -- response gain

%% define starting points and bounds on parameters

%params:    [weight, shift, gain]

% starting point
opts.x0   = [0.5,    0.06,  2];
% lower bounds
opts.lb   = [0,      0,     0.01];
% upper bounds
opts.ub   = [1,      0.1,   200];

% plausible lower and upper bounds (required for bads search algorithm)
opts.plb  = [0.1,    0.01,  0.5];
opts.pub  = [0.9,    0.08,  100];       
        
%% write out json
fname = fullfile(tdeRootPath, 'temporal_models', 'TTCmodel.json');

%json_opts.indent = '    '; % this just makes the json file look prettier 
%jsonwrite(json_fname,opts,json_opts);

opts_json = savejson('',opts,fname);