P_D = rand;
clutter_intensity = rand/100;
%Create nonlinear measurement model (range/bearing)
sigma_r = 5;
sigma_b = pi/180;
s = [300;400];
meas_model = measmodel.rangebearingmeasmodel(sigma_r, sigma_b, s);
%Set Poisson RFS
PPP.w = log(rand(3,1));
PPP.states = repmat(struct('x',rand(5,1),'P',eye(5)),[3,1]);

PPP.states(1).x = [0; 0; 5; 0; pi/180];
PPP.states(2).x = [20; 20; -20; 0; pi/90];
PPP.states(3).x = [-20; 10; -10; 0; pi/360];

indices = [true;false;true];
z = rand(2,1);

% Run reference solution.
PMBM_ref = reference.PMBMfilter();
PMBM_ref.density = feval(@GaussianDensity);
PMBM_ref.paras.PPP.w = PPP.w;
PMBM_ref.paras.PPP.states = PPP.states;
[Bern_ref, lik_new_ref] = PPP_detected_update(PMBM_ref,indices,z,meas_model,P_D,clutter_intensity);
r_ref = Bern_ref.r;
state_ref = Bern_ref.state;
x_ref = state_ref.x;
P_ref = state_ref.P;

% Run learner solution.
PMBM = PMBMfilter();
%Create density handle
PMBM.density = feval(@GaussianDensity);
PMBM.paras.PPP.w = PPP.w;
PMBM.paras.PPP.states = PPP.states;
[Bern, lik_new] = PPP_detected_update(PMBM,indices,z,meas_model,P_D,clutter_intensity);
r = Bern.r;
state = Bern.state;
x = state.x;
P = state.P;

% Compare.
assessVariableEqual('r', r_ref,'RelativeTolerance',0.0001,'Feedback','Parameter Bern.r is incorrect, please modify your code and submit again.');
assessVariableEqual('x', x_ref,'AbsoluteTolerance',0.001,'Feedback','Parameter Bern.state.x is incorrect, please modify your code and submit again.');
assessVariableEqual('P', P_ref,'AbsoluteTolerance',0.001,'Feedback','Parameter Bern.state.P is incorrect, please modify your code and submit again.');
assessVariableEqual('lik_new', lik_new_ref,'RelativeTolerance',0.0001,'Feedback','Parameter lik_new is incorrect, please modify your code and submit again.');