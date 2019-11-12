nbirths = 1;
K = 10;
initial_state.x = [0; 0; 5; 5];
ground_truth = modelgen.groundtruth(nbirths,initial_state.x,1,K+1,K);

T = 1;
sigma_q = 2;
motion_model = motionmodel.cvmodel(T,sigma_q);

P_D = 1;
lambda_c = 60;
range_c = [-200 200;-200 200];
sensor_model = modelgen.sensormodel(P_D,lambda_c,range_c);

sigma_r = 10;
meas_model = measmodel.cvmeasmodel(sigma_r);

objectdata = objectdatagen(ground_truth,motion_model,0);
     
measdata = measdatagen(objectdata,sensor_model,meas_model);

gating_size = chi2inv(0.9999,meas_model.d);

for i = 1:K-1
     state.x = objectdata.X{i};
     state.P = eye(motion_model.d);
     state_pred = GaussianDensity.predict(state, motion_model);
     % Run learner solution.
     [z_ingate, meas_in_gate] = GaussianDensity.ellipsoidalGating(state_pred, measdata{i+1}, meas_model, gating_size);
     m_ingate = size(z_ingate,2);
     if m_ingate > 1
         z_ingate = sort(z_ingate,2);
     end
     
     % Run reference solution.
     [z_ingate_ref, meas_in_gate_ref] = reference.GaussianDensity.ellipsoidalGating(state_pred, measdata{i+1}, meas_model, gating_size);
     m_ingate_ref = size(z_ingate_ref,2);
     if m_ingate_ref > 1
         z_ingate_ref = sort(z_ingate_ref,2);
     end
     
     % Compare.
     assessVariableEqual('meas_in_gate', meas_in_gate_ref,'Feedback','The output boolean vector does not match the reference solution. Your implementation is incorrect, please modify your code and submit again.');
     assert(abs(m_ingate-m_ingate_ref) == 0, 'Number of measurements inside the gate: Expected %f, got %f. Your implementation is incorrect, please modify your code and submit again.', m_ingate_ref, m_ingate);
     assessVariableEqual('z_ingate', z_ingate_ref,'Feedback','The output measurements do not match the reference solution. Your implementation is incorrect, please modify your code and submit again.');
end
