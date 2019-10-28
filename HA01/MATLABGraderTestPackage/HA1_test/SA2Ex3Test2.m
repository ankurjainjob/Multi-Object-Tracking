threshold = 50;

hypothesis.x = 0;
multiHypotheses = repmat(hypothesis,100,1);
hypothesesWeight = rand(100,1);
hypothesesWeight = log(hypothesesWeight/sum(hypothesesWeight));

% Run learner solution.
[hypothesesWeight_hat, multiHypotheses_hat] = hypothesisReduction.cap(hypothesesWeight, multiHypotheses, threshold); 

hypothesesWeight_hat = sort(hypothesesWeight_hat);

% Run reference solution.
[hypothesesWeight_hat_ref, multiHypotheses_hat_ref] = reference.hypothesisReduction.cap(hypothesesWeight, multiHypotheses, threshold); 

hypothesesWeight_hat_ref = sort(hypothesesWeight_hat_ref);

% Compare.
assessVariableEqual('hypothesesWeight_hat', hypothesesWeight_hat_ref,'RelativeTolerance',0.03,'Feedback','Your implementation is incorrect, please modify your code and submit again.');