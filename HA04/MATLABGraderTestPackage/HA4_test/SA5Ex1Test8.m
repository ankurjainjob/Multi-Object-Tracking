%Create nonlinear motion model (coordinate turn)
T = 1;
sigmaV = 1;
sigmaOmega = pi/180;
motion_model = motionmodel.ctmodel(T,sigmaV,sigmaOmega);
%Set birth model
birth_model = repmat(struct('w',log(rand),'x',rand(5,1),'P',diag([1 1 1 1*pi/90 1*pi/90].^2)),[1,4]);
%Set probability of existence
P_S = rand;
%Set Poisson RFS
PPP.w = log(rand(3,1));
PPP.states = repmat(struct('x',rand(5,1),'P',eye(5)),[3,1]);
PPP.states(1).x = [0; 0; 5; 0; pi/180];
PPP.states(2).x = [20; 20; -20; 0; pi/90];
PPP.states(3).x = [-20; 10; -10; 0; pi/360];
%Set Bernoulli RFS
tt = cell(2,1);
tt{1} = repmat(struct('r',rand,'state',struct('x',rand(5,1),'P',eye(5))),[3,1]);
tt{2} = repmat(struct('r',rand,'state',struct('x',rand(5,1),'P',eye(5))),[2,1]);

% Run reference solution.
PMBM_ref = reference.PMBMfilter();
PMBM_ref.density = feval(@GaussianDensity);
PMBM_ref.paras.PPP.w = PPP.w;
PMBM_ref.paras.PPP.states = PPP.states;
PMBM_ref.paras.MBM.tt = tt;
PMBM_ref = PMBM_predict(PMBM_ref,P_S,motion_model,birth_model);

tt_ref = PMBM_ref.paras.MBM.tt;
[PPP_w_ref,I] = sort(PMBM_ref.paras.PPP.w);
PPP_states_ref = PMBM_ref.paras.PPP.states(I);

% Run learner solution.
PMBM = PMBMfilter();
%Create density handle
PMBM.density = feval(@GaussianDensity);
PMBM.paras.PPP.w = PPP.w;
PMBM.paras.PPP.states = PPP.states;
PMBM.paras.MBM.tt = tt;
PMBM = PMBM_predict(PMBM,P_S,motion_model,birth_model);

tt = PMBM.paras.MBM.tt;
[PPP_w,I] = sort(PMBM.paras.PPP.w);
PPP_states = PMBM.paras.PPP.states(I);

% Compare.
assessVariableEqual('PPP_w', PPP_w_ref,'RelativeTolerance',0.0001,'Feedback','Parameter obj.paras.PPP.w is incorrect, please modify your code and submit again.');
assessVariableEqual('PPP_states', PPP_states_ref,'RelativeTolerance',0.001,'Feedback','Parameter obj.paras.PPP.states is incorrect, please modify your code and submit again.');
assessVariableEqual('tt', tt_ref,'Feedback','Parameter obj.paras.MBM.tt is incorrect, please modify your code and submit again.');