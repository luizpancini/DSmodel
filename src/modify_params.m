function params = modify_params(params,data,source,model,gust_ind,case_now)

switch source
    case {"gust_test","hargen_test"}
        if isfield(data,'case_now') && data.source == "gust_test", case_now = data.case_now; end
        if isnumeric(case_now) && ismember(case_now,[14:20])
            % Better c_n_alpha fit
            params.c_n_alpha = 2*pi/params.beta*1.07;
        end
end

end