N = 20;
x = cell(N,1);
y = cell(N,1);

for i = 1:N
    x{i} = randn(4,1);
    y{i} = randn(4,1);
end


% Run learner solution.
root_mean_square_error = RMSE(x,y);

% Run reference solution.
root_mean_square_error_ref = reference.RMSE(x,y);

% Compare.
assert(abs(root_mean_square_error-root_mean_square_error_ref) < 1e-3, 'Expected %f, got %f. Your implementation is incorrect, please modify your code and submit again.', root_mean_square_error_ref, root_mean_square_error);
% assessVariableEqual('root_mean_square_error', root_mean_square_error_ref,'RelativeTolerance',0.03,'Feedback','Your implementation is incorrect, please modify your code and submit again.');