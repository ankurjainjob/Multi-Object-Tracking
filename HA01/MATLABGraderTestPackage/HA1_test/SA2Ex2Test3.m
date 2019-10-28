numGaussians = 10;

%Moment matching
numGaussians = 5;
nDim = 4;
states = [];
for i = 1:numGaussians
    states(i,1).x = 5*randn(nDim,1);
    S = randn(nDim);
    states(i,1).P = S*S';
end

w = abs(randn(numGaussians,1));
w = log(w/sum(w));

% Run learner solution.
state = GaussianDensity.momentMatching(w, states);

% Run reference solution.
state_ref = reference.GaussianDensity.momentMatching(w, states);

% Compare.
tol = 1e-4;
assert(norm(state.x-state_ref.x) < tol, 'The merged mean is not correctly computed. The Euclidean norm between your solution and the reference solution is %f. Your implementation is incorrect, please modify your code and submit again.', norm(state.x-state_ref.x));
assert(norm(state.P-state_ref.P) < tol, 'The merged covariance is not correctly computed. The Euclidean norm between your solution and the reference solution is %f. Your implementation is incorrect, please modify your code and submit again.', norm(state.P-state_ref.P));
% assessVariableEqual('state', state_ref,'RelativeTolerance',0.03,'Feedback','Your implementation is incorrect, please modify your code and submit again.');