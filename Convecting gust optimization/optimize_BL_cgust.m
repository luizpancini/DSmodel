function E = optimize_BL_cgust(p,INPUTS,varargin)

% Handle inputs
plotting = 0;           % TF to plot  
if ~isempty(varargin)
    plotting = varargin{1};
end

%% Loop over cases
% Pre-process
if ~isrow(INPUTS.cases), INPUTS.cases = INPUTS.cases'; end
INPUTS.base_name = "Gust test case ";

% Initialize error sum
E = 0;

% Loop
j = 0;
for case_now = INPUTS.cases
    
    j = j+1;
    
    %% Initialization
    [data,params,x0,y0,ODE_options] = initialize_model(case_now,INPUTS,0,1E-4,0);
    % Set convection gust parameters
    params = set_cgust_params(params,p);
    
    %% ODEs solver
    [t,x,y] = BLgust_RKF45(params.tspan,x0,y0,params,ODE_options);
    
    %% Get output variables
    outputs = BLgust_output_vars(x,y,t);
    
    %% Interpolate data and find normalized errors
    % Interpolate reference data
    N = length(t); tau = params.U/params.b*t;
    cl_interp = interp1(data.time_exp_cl,data.clt_exp,tau,'pchip');
    % Sum of squares of the errors
    eps_cl = sum((outputs.c_l(:)-cl_interp(:)).^2);    
    % Normalized cumulative error sum
    E = E + sqrt(eps_cl/N)/abs(max(cl_interp)-min(cl_interp));
    % Plot
%     figure;plot(tau,outputs.c_l,'r-',tau,cl_interp,'k-o');grid
    
    %% Plots
    if plotting
        % Unpack
        authors = data.authors;
        model = INPUTS.model;
        source = INPUTS.source;
        base_name = INPUTS.base_name;
        t = outputs.t;
        c_l = outputs.c_l;
        model_name = "Present model";
        % Initialize
        [h,plot_opt] = init_plotter(struct(),base_name,case_now,source,params);
        % Plot
        time_plotter('Lift coefficient','cl x t',t,c_l,t,cl_interp,data.time_exp_cl,data.clt_exp,data.time_mod_cl,data.clt_mod,data.time_cycles,data.cl_cycles,model,outputs,params,[],authors,source,h.tabgp,model_name,plot_opt);
    end  
    
end
   
end