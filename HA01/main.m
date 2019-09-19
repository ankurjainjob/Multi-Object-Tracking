clear all; clc; close all;

%Call ellipsoidal gating
%Choose object state
state_pred.x = [0;0;0;0];
state_pred.P = eye(4);
z = rand(2,100);
%Create measurement model
sigma_r = 1;
meas_model = measmodel.cvmeasmodel(sigma_r);
%Choose gating size
gating_size = chi2inv(0.99,meas_model.d);
%Ellipsoidal gating
[z_ingate, meas_in_gate] = GaussianDensity.ellipsoidalGating(state_pred, z, meas_model, gating_size);

%Calculate predicted likelihood
predict_likelihood = GaussianDensity.predictedLikelihood(state_pred,z,meas_model);

%Moment matching
numGaussians = 5;
for i = 1:numGaussians
    states(i,1).x = 10*randn(4,1);
    states(i,1).P = 10*eye(4);
end
w = rand(numGaussians,1);
w = log(w/sum(w));
state = GaussianDensity.momentMatching(w, states);

%Mixture reduction
threshold = 2;
[w_hat,states_hat] = GaussianDensity.mixtureReduction(w,states,threshold);

