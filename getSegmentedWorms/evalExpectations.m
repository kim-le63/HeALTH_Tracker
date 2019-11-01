function score = evalExpectations(results, expectations)
% Used as an input to getThreshold
% Pulls WormSizes out of the output of findWorm and compares them to
%    expectations, so that a threshold can be determined based on
%    the size of the worms at a certain threshold
assert((isstruct(expectations) && isstruct(results)) && all(...
    size(expectations) == size(results)))
expectationFields = fieldnames(expectations);
score = zeros(size(results));
for i = 1:numel(expectationFields) % Need to make sure no missing results
    expField = expectationFields{i};
    if ~contains(expField,'Std') && ~contains(expField, 'Report') && ...
            ~contains(expField, 'Confidence')
        % It's a valid field that should have a result
        if ~isfield(results, expField)
            error('WormHotel:MissingExpectationResult', ...
                ['No result reported for expectation ' expField])
        else
            if isfield(expectations, [expField 'Std'])
                score = score + arrayfun(@(expEl, resEl) ... 
                    abs(expEl.(expField) - resEl.(expField)) ./ ...
                    expEl.([expField 'Std']), expectations, results);
            elseif isfield(expectations, [expField 'Confidence'])
                score = score + arrayfun(@(expEl, resEl) ...
                    abs(expEl.(expField) - resEl.(expField)) .* ...
                    expEl.([expField 'Confidence']) ./ expEl.(expField), ...
                    expectations, results);
            else
                warning(['field ' expField ' had no method of comparison'])
            end
        end
    end
end
end
%{
results = orderfields(results);
expectations = orderfields(expectations);
scoreSet = fieldnames(results);
for i = numel(scoreSet):-1:1
    if ~any(cellfun(@(cell) strcmp(cell, [scoreSet{i} 'Std']), ...
            fieldnames(expectations))) || ... % no 'Std' field in expect.
            ~any(cellfun(@(cell) strcmp(cell, scoreSet{i}), ...
            fieldnames(expectations))) % no expectation to compare against
        if ~contains(scoreSet{i}, 'Report')
            warning(['WormHotel:FieldNotReported', ...
                'No way to analyze result ' scoreSet{i}])
        end
        try results = rmfield(results, scoreSet{i}); catch, ''; end
        scoreSet(i) = [];
    end
end
results = squeeze(struct2cell(orderfields(results)));
% Keep only expectations that can be evaluated
expectations = rmfield(expectations, setdiff(fieldnames(expectations), scoreSet));
expectations = squeeze(struct2cell(orderfields(expectations)));
score = cellfun(@minus, results, expectations, 'UniformOutput', false);
score = cellfun(@(cell) score(:).(cell) ./ expectations(:).(...
    [cell 'Std']), scoreSet);
if isempty(score)
    score = 0;
else
    score = sum(abs([score{:}]));
end
                %}