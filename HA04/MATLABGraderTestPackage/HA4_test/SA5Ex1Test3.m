%Set detection probability
P_D = rand;
%Create hypothesis tree
Bern.r = rand;
Bern.state.x = rand(5,1);
Bern.state.P = eye(5);
tt_entry = [1,1];
z = rand(2,3);
%Create nonlinear measurement model (range/bearing)
sigma_r = 5;
sigma_b = pi/180;
s = [300;400];
meas_model = measmodel.rangebearingmeasmodel(sigma_r, sigma_b, s);

%Run learner solution.
PMBM_ref = reference.PMBMfilter();
PMBM_ref.density = feval(@GaussianDensity);
PMBM_ref.paras.MBM.tt{1} = Bern;
lik_detected_ref = Bern_detected_update_lik(PMBM_ref,tt_entry,z,meas_model,P_D);

% Run reference solution.
PMBM = PMBMfilter();
PMBM.density = feval(@GaussianDensity);
PMBM.paras.MBM.tt{1} = Bern;
lik_detected = Bern_detected_update_lik(PMBM,tt_entry,z,meas_model,P_D); 

% Compare.
assessVariableEqual('lik_detected', lik_detected_ref,'RelativeTolerance',0.0001,'Feedback','Parameter lik_detected is incorrect, please modify your code and submit again.');