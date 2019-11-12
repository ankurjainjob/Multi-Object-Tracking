%Create nonlinear motion model (coordinate turn)
T = 1;
sigmaV = 1;
sigmaOmega = pi/180;
motion_model = motionmodel.ctmodel(T,sigmaV,sigmaOmega);
%Set probability of existence
P_S = rand;
%Set Bernoulli RFS
Bern.r = rand;
Bern.state.x = rand(5,1);
Bern.state.P = eye(5);

% Run reference solution.
PMBM_ref = reference.PMBMfilter();
PMBM_ref.density = feval(@GaussianDensity);
Bern_ref = Bern_predict(PMBM_ref,Bern,motion_model,P_S); 
r_ref = Bern_ref.r;
state_ref = Bern_ref.state;

% Run learner solution.
PMBM = PMBMfilter();
PMBM.density = feval(@GaussianDensity);
Bern = Bern_predict(PMBM,Bern,motion_model,P_S);
r = Bern.r;
state = Bern.state;

% Compare.
assessVariableEqual('r', r_ref,'Feedback','Parameter Bern.r is incorrect, please modify your code and submit again.');
assessVariableEqual('state', state_ref,'Feedback','Parameter Bern.state is incorrect, please modify your code and submit again.');