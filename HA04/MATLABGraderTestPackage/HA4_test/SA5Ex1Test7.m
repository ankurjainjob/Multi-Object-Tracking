P_D = rand;
PPP.w = log(rand(3,1));

% Run reference solution.
PMBM_ref = reference.PMBMfilter();
PMBM_ref.paras.PPP.w = PPP.w;
PMBM_ref = PPP_undetected_update(PMBM_ref,P_D);

PPP_w_ref = PMBM_ref.paras.PPP.w;

% Run learner solution.
PMBM = PMBMfilter();
PMBM.paras.PPP.w = PPP.w;
PMBM = PPP_undetected_update(PMBM,P_D);

PPP_w = PMBM.paras.PPP.w;

% Compare.
assessVariableEqual('PPP_w', PPP_w_ref,'Feedback','Parameter PPP.w is incorrect, please modify your code and submit again.');