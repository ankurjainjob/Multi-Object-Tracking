%Create hypothesis tree
Bern.r = rand;
Bern.state.x = rand(5,1);
Bern.state.P = eye(5);
tt_entry = [1,1];
z = rand(2,1);
%Create nonlinear measurement model (range/bearing)
sigma_r = 5;
sigma_b = pi/180;
s = [300;400];
meas_model = measmodel.rangebearingmeasmodel(sigma_r, sigma_b, s);

%Run learner solution.
PMBM_ref = reference.PMBMfilter();
PMBM_ref.density = feval(@GaussianDensity);
PMBM_ref.paras.MBM.tt{1} = Bern;
Bern_ref = Bern_detected_update_state(PMBM_ref,tt_entry,z,meas_model);
r_ref = Bern_ref.r;
state_ref = Bern_ref.state;

% Run reference solution.
PMBM = PMBMfilter();
PMBM.density = feval(@GaussianDensity);
PMBM.paras.MBM.tt{1} = Bern;
Bern = Bern_detected_update_state(PMBM,tt_entry,z,meas_model);
r = Bern.r;
state = Bern.state;

% Compare.
assessVariableEqual('r', r_ref,'Feedback','Parameter Bern.r is incorrect, please modify your code and submit again.');
assessVariableEqual('state', state_ref,'Feedback','Parameter Bern.state is incorrect, please modify your code and submit again.');
assessVariableEqual('Bern', Bern_ref,'RelativeTolerance',0.01,'Feedback','Parameter Bern is incorrect, please modify your code and submit again.');