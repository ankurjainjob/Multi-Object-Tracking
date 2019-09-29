function root_mean_square_error = RMSE(state_sequence1,state_sequence2)
    %RMSE calculates the root mean square error between state_sequence1 and
    %state_sequence2 
    %INPUT: state_sequence1: a cell array of size (total tracking time x 1),
    %       each cell contains an object state mean vector of size
    %       (state dimension x 1)
    %       state_sequence2: a cell array of size (total tracking time x 1),
    %       each cell contains an object state mean vector of size
    %       (state dimension x 1)
    %OUTPUT:root_mean_square_error: root mean square error --- scalar
    
    root_mean_square_error = mean( (cell2mat(state_sequence1)-cell2mat(state_sequence2)).^2 ) ^0.5
end



%% test function

% N = 20;
% x = cell(N,1);
% y = cell(N,1);
% 
% for i = 1:N
%     x{i} = randn(4,1);
%     y{i} = randn(4,1);
% end
% 
% mse(cell2mat(x) - cell2mat(y))
% 
% root_mean_square_error = RMSE(x,y);