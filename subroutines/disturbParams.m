% JN Kather 2017
% randomly disturbs selected parameters by a given percentage
function allParams = disturbParams(allParams,targetvars,magnitude)

for i=1:numel(targetvars)
    crand = 1 + magnitude * 2 * (rand() - 0.5); % varies between 1-magnitude and 1+magnitude
    allParams.(targetvars{i}) = allParams.(targetvars{i}) * crand;
end

end