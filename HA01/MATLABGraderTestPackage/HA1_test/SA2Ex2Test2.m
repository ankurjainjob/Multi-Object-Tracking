nbirths = 1;
K = 10;
initial_state.x = [0; 0; 5; 5];
ground_truth = modelgen.groundtruth(nbirths,initial_state.x,1,K+1,K);

T = 1;
sigma_q = 2;
motion_model = motionmodel.cvmodel(T,sigma_q);

P_D = 1;
lambda_c = 10;
range_c = [-50 50;-50 50];
sensor_model = modelgen.sensormodel(P_D,lambda_c,range_c);

sigma_r = 5;
meas_model = measmodel.cvmeasmodel(sigma_r);

objectdata = objectdatagen(ground_truth,motion_model,0);
     
measdata = measdatagen(objectdata,sensor_model,meas_model);


for i = 1:K
     state.x = objectdata.X{i};
     state.P = eye(motion_model.d);
     % Run learner solution.
     predict_likelihood = GaussianDensity.predictedLikelihood(state,measdata{i},meas_model);
     % Run reference solution.
     predict_likelihood_ref = reference.GaussianDensity.predictedLikelihood(state,measdata{i},meas_model);
     % Compare.
     assert(norm(predict_likelihood-predict_likelihood_ref) < 1e-3, 'The Euclidean norm between your solution and the reference solution is %f. Your implementation is incorrect, please modify your code and submit again.', norm(predict_likelihood-predict_likelihood_ref));
%      assessVariableEqual('predict_likelihood', predict_likelihood_ref,'RelativeTolerance',0.01,'Feedback','Your implementation is incorrect, please modify your code and submit again.');
end
