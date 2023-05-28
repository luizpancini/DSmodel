function s = variables2struct(s,varargin)

% Add fields given by input variables to structure s
for i=1:nargin-1
    s.(inputname(i+1)) = varargin{i}; % Overwrite field if existent
end

end