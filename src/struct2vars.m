function varargout = struct2vars(s)

% Get structure's field's names
fields = fieldnames(s);

% Output structure's fields
varargout = cell(length(fields),1);
for i=1:length(fields)
    varargout{i} = s.(fields{i});
end

end