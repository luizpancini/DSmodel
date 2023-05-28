function [wg,vg] = gust_kin_comp(t,U,U_0,lambda_g0,wg_fun)

% Gust-induced chordwise velocity component
lambda_g = U/(U+U_0*(1-lambda_g0)/lambda_g0); % Current gust speed ratio
vg = U/lambda_g-U;

% Gust-induced normal velocity component
if ~isempty(wg_fun)
    wg = wg_fun(t);  
else
    wg = 0; 
end

end