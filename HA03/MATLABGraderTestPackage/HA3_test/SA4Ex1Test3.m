numGaussian = 20;

% Get generator settings.
s = rng;

% Run learner solution
PPP = PHDfilter();
PPP.density = feval(@GaussianDensity);
PPP.paras.w = log(rand(numGaussian,1));
for i = 1:numGaussian
    PPP.paras.states(i,1).x = rand(2,1);
end

estimates = PHD_estimator(PPP);

[~,I] = sort(estimates(1,:));
estimates = estimates(:,I);
num_objects = size(estimates,2);

% Restore previous generator settings.
rng(s);

% Run reference solution.
PPP_ref = reference.PHDfilter();
PPP_ref.density = feval(@GaussianDensity);
PPP_ref.paras.w = log(rand(numGaussian,1));
for i = 1:numGaussian
    PPP_ref.paras.states(i,1).x = rand(2,1);
end

estimates_ref = PHD_estimator(PPP_ref);


[~,I] = sort(estimates_ref(1,:));
estimates_ref = estimates_ref(:,I);
num_objects_ref = size(estimates_ref,2);

% Compare.
assert(num_objects == num_objects_ref,'The number of estimated objects: Expected %d, got %d.', num_objects_ref, num_objects);
assessVariableEqual('estimates', estimates_ref,'RelativeTolerance',0.001,'Feedback','Your implementation is incorrect, please modify your code and submit again.');