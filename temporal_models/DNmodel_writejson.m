function opts_json = DNmodel_writejson()
% Writes out a json file with parameters for model fitting.
%
% DN model parameters:
%          1. tau1 -- time to peak for positive IRF (seconds)
%          2. weight -- ratio of negative to positive IRFs (set to 0
%               for uniphasic irf)
%          3. tau2 -- time window of adaptation (seconds)
%          4. n -- exponent
%          5. sigma -- semi-saturation constant
%          6. shift -- time between stimulus onset and when the signal reaches
%               the cortex (seconds)
%          7. scale -- response gain.

%% define starting points and bounds on parameters

% params:   [t1,   w, t2,   n,   sigma, shift, scale]

% starting point
opts.x0   = [0.03, 0, 0.07, 1.5, 0.15, 0.06, 2];  
% lower bounds
opts.lb   = [0.01, 0, 0.01, 1,   0,    0,    0.01];
% upper bounds
opts.ub   = [1,    1, 2,    5,   1,    0.1,  200];

% plausible lower and upper bounds (required for bads search algorithm)
opts.plb  = [0.1, 0,   0.1, 1.5, 0.01, 0.01, 0.5];
opts.pub  = [0.9, 0.5, 1,   3,   0.5,  0.08, 100];


%% write out json
fname = fullfile(tdeRootPath, 'temporal_models', 'DNmodel.json');

%json_opts.indent = '    '; % this just makes the json file look prettier 
%jsonwrite(json_fname,opts,json_opts);

opts_json = savejson('',opts,fname);