function opts_json = tde_writeParamsToJson(modelName)
% Writes out a json file with parameters for model fitting.%

%% define starting points and bounds on parameters

switch modelName
    
    case {'DN', 'DNCASCADE'}
        
        opts.params = 't1,w,t2,n,sigma,shift,scale';
        opts.x0     = [0.03, 0, 0.07, 1.5, 0.15, 0.06, 2];    % starting point
        opts.lb     = [0.01, 0, 0.01, 1,   0,    0,    0.01]; % lower bounds
        opts.ub     = [1,    1, 2,    5,   1,    0.1,  200];  % upper bounds
        opts.plb    = [0.1, 0,   0.1, 1.5, 0.01, 0.01, 0.5];  % plausible lower bound (required for bads search algorithm)
        opts.pub    = [0.9, 0.5, 1,   3,   0.5,  0.08, 100];  % plausible upper bound (required for bads search algorithm)
    
    case {'TTC', 'TTCSTIG17'}
        
        opts.params = 'weight,shift,gain';
        opts.x0   = [0.5,    0.06,  2];    % starting point
        opts.lb   = [0,      0,     0.01]; % lower bounds
        opts.ub   = [1,      0.1,   200];  % upper bounds
        opts.plb  = [0.1,    0.01,  0.5];  % plausible lower bound (required for bads search algorithm)
        opts.pub  = [0.9,    0.08,  100];  % plausible upper bound (required for bads search algorithm)
        
    case {'TTCSTIG19'}
        
        opts.params = 'weight,shift,scale,tau,k_on,k_off,lambda,alpha';
        opts.x0   = [0.5,    0.06,  2,    4.93, 3,      3,     0.1,     1];  % starting point
        opts.lb   = [0,      0,     0.01, 0,    0.01,   0.01,  0.001,   1];      % lower bounds
        opts.ub   = [1,      0.1,   200,  100,  10,     10,    100,     100000]; % upper bounds
        opts.plb  = [0.1,    0.01,  0.5,  1,    0.01,   0.01,  0.001,   1];      % plausible lower bound (required for bads search algorithm)
        opts.pub  = [0.9,    0.08,  100,  100,   10,     10,   100,    100000]; % plausible upper bound (required for bads search algorithm)
        
     case {'HEEGER92'}
        
        opts.params = 'tau1,shift,sigma,alpha,n';
        opts.x0   = [0.5,   0.06, 0.1, 0.01,  1    ];  % starting point
        opts.lb   = [0.01,  0,    0,   0,     0    ];  % lower bounds
        opts.ub   = [1,     0.1,  2,   1,     10   ];  % upper bounds
        opts.plb  = [0.1,   0.01, 0.01,0,     0.01 ];  % plausible lower bound (required for bads search algorithm)
        opts.pub  = [0.9,   0.08, 2,   1,     2    ];  % plausible upper bound (required for bads search algorithm)
       
end

%% write out json
fname = fullfile(tdeRootPath, 'temporal_models', sprintf('%s.json', modelName));
opts_json = savejson('',opts,fname);