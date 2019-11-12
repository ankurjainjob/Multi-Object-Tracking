%Set detection probability
P_D = rand;
%Create hypothesis tree
Bern.r = rand;
tt_entry = [1,1];

% Run reference solution.
PMBM_ref = reference.PMBMfilter();
PMBM_ref.paras.MBM.tt{1} = Bern;
[Bern_ref, lik_undetected_ref] = Bern_undetected_update(PMBM_ref,tt_entry,P_D);
r_ref = Bern_ref.r;

% Run learner solution.
PMBM = PMBMfilter();
PMBM.paras.MBM.tt{1} = Bern;
[Bern, lik_undetected] = Bern_undetected_update(PMBM,tt_entry,P_D);
r = Bern.r;

% Compare.
assessVariableEqual('r', r_ref,'Feedback','Parameter Bern.r is incorrect, please modify your code and submit again.');
assessVariableEqual('lik_undetected', lik_undetected_ref,'Feedback','Parameter lik_undetected is incorrect, please modify your code and submit again.');