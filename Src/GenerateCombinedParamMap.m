%
% This function returns a map of three parameters, the combinations of
% which are arranged along one dimension. This allows a triple-nested for 
% loop to be converted to a single for loop where its parameters can be
% indexed by a single integer iterator, as required by parfor loops.
%
% Use case: when using parfor, this method maximises the number of cores
% which can be used for processing.
%
% Inputs: param_1-3 should have the shape (1xN), but don't need to match
% e.g. param_1 = [1 2 3], param_2 = [-2 -4 -6 -8], param_3 = [10 20].
%
% Usage example:
% 
% combined_param_map = GenerateCombinedParamMap(a_params, b_params, c_params);
%
% parfor combined_index = 1:size(combined_param_map, 2)
%     a = combined_param_map(1, combined_index);
%     b = combined_param_map(2, combined_index);
%     c = combined_param_map(3, combined_index);
% 
%     Foo(a, b, c);
% end

function combined_param_map = GenerateCombinedParamMap(param_1, param_2, param_3)
    combined_param_map = zeros(3, size(param_1, 2) * size(param_2, 2) * size(param_3, 2));
    combined_index = 1;
    
    for param_1_index = param_1
        for param_2_index = param_2
            for param_3_index = param_3
                combined_param_map(:,combined_index) = [param_1_index param_2_index param_3_index];
                combined_index = combined_index + 1;
            end
        end
    end
end