function bool = structContains(structIn, fieldIn)
bool = isa(structIn, 'struct') && isfield(structIn, fieldIn) && ...
    ~isempty(structIn.(fieldIn));
end