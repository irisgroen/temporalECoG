function opts_json = tde_writeParamsToJson(modelName)
% Writes out a json file with parameters for model fitting.
% 2020 Iris Groen

%% define starting points and bounds on parameters

switch modelName
    
    case {'DN', 'DNCASCADE'}
        
        opts.params = 't1,w,t2,n,sigma,shift,scale';
        opts.x0     = [0.03, 0,   0.07, 1.5, 0.15, 0.06, 2];    % starting point
        opts.lb     = [0.01, 0,   0.01, 1,   0,    0,    0.01]; % lower bounds
        opts.ub     = [1,    1,   2,    5,   1,    0.1,  200];  % upper bounds
        opts.plb    = [0.1,  0,   0.1,  1.5, 0.01, 0.01, 0.5];  % plausible lower bound (required for bads search algorithm)
        opts.pub    = [0.9,  0.5, 1,    3,   0.5,  0.08, 100];  % plausible upper bound (required for bads search algorithm)
    
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
        
     case {'HEEGER92', 'HEEGER93'}
        
        opts.params = 'tau1,shift,sigma,alpha,rmax';
        opts.x0   = [0.02,  0.06, 0.1, 0.01,  10];  % starting point
        opts.lb   = [0.01,  0,    0,   0,     0.01]; % lower bounds
        opts.ub   = [1,     0.1,  2,   1,     200];  % upper bounds
        opts.plb  = [0.1,   0.01, 0.01,0,     0.1];  % plausible lower bound (required for bads search algorithm)
        opts.pub  = [0.9,   0.08, 2,   1,     50];   % plausible upper bound (required for bads search algorithm)
     
    case {'FIR'} 
        
        n = 300;
        opts.params = '';
        opts.x0   = zeros(1,n); % starting point
        opts.lb   = opts.x0-100;
        opts.ub   = opts.x0+100;
        opts.plb   = opts.x0-100;
        opts.pub   = opts.x0+100;
        
	case {'LINEAR', 'LINEAR_RECTH', 'LINEAR_RECTF'} 
         
        opts.params = 'tau1, tau2, n_irf, weight, shift, scale';
        opts.x0     = [0.01, 0.03, 2,     0,      0.01,  2];    % starting point
        opts.lb     = [0.01, 0.01, 2,     0,      0,     0.01]; % lower bounds
        opts.ub     = [10,   10,   20,    1,      1,     200];  % upper bounds
        opts.plb    = [0.1,  0.1,  2,     0.01,   0.01,  0.5];  % plausible lower bound (required for bads search algorithm)
        opts.pub    = [0.9,  0.9,  10,    0.5,    0.08,  100];  % plausible upper bound (required for bads search algorithm)
    
	case {'LINEAR_RECTH_EXP', 'LINEAR_RECTF_EXP'} 
        
        opts.params = 'tau1, tau2, n_irf, weight, shift, scale, n';
        opts.x0     = [0.01, 0.03, 2,     0,      0.01,   2,    1.5];   % starting point
        opts.lb     = [0.01, 0.01, 2,     0,      0,      0.01, 1];     % lower bounds
        opts.ub     = [10,   10,   20,    1,      1,      200,  5];     % upper bounds
        opts.plb    = [0.1,  0.1,  2,     0.01,   0.01,   0.5,  1.5];   % plausible lower bound (required for bads search algorithm)
        opts.pub    = [0.9,  0.9,  10,    0.5,    0.08,   100,  3];     % plausible upper bound (required for bads search algorithm)
    
    case {'LINEAR_RECTH_EXP_NORM', 'LINEAR_RECTF_EXP_NORM'} 
        
        opts.params = 'tau1, tau2, n_irf, weight, shift, scale, n,   sigma';
        opts.x0     = [0.01, 0.03, 2,     0,      0.01,   2,    1.5, 0.15];    % starting point
        opts.lb     = [0.01, 0.01, 2,     0,      0,      0.01, 1,   0];       % lower bounds
        opts.ub     = [10,   10,   20,    1,      1,      200,  5,   1];       % upper bounds
        opts.plb    = [0.1,  0.1,  2,     0.01,   0.01,   0.5,  1.5  0.01];    % plausible lower bound (required for bads search algorithm)
        opts.pub    = [0.9,  0.9,  10,    0.5,    0.08,   100,  3,   0.5];     % plausible upper bound (required for bads search algorithm)
    
	case {'LINEAR_RECTH_EXP_NORM_DELAY', 'LINEAR_RECTF_EXP_NORM_DELAY'} 
        
        opts.params = 'tau1, tau2, n_irf,  weight, shift, scale, n,   sigma, tau_a';
        opts.x0     = [0.01, 0.03, 2,      0,     0.01,   2,     1.5, 0.15,  0.07];    % starting point
        opts.lb     = [0.01, 0.01, 2,      0,     0,      0.01,  1,   0,     0.01];    % lower bounds
        opts.ub     = [10,   10,   20,     1,     1,      200,   5,   1,     2];       % upper bounds
        opts.plb    = [0.1,  0.1,  2,      0.01,  0.01,   0.5,   1.5  0.01,  0.1];     % plausible lower bound (required for bads search algorithm)
        opts.pub    = [0.9,  0.9,  10,     0.5,   0.08,   100,   3,   0.5,   1];       % plausible upper bound (required for bads search algorithm)
   
end

%% write out json
fname = fullfile(tdeRootPath, 'temporal_models', sprintf('%s.json', modelName));
opts_json = savejson('',opts,fname);