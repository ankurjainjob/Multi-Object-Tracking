w = sort(rand(3,1));
ht = [1 0 1;
      1 1 1;
      0 1 1];
r = rand(3,1);
threshold = mean(r);
tt = cell(3,1);
tt{1}.r = r(1);
tt{1}.state.x = rand(4,1);
tt{2}.r = r(2);
tt{2}.state.x = rand(4,1);
tt{3}.r = r(3);
tt{3}.state.x = rand(4,1);

% Run reference solution.
PMBM_ref = reference.PMBMfilter();
PMBM_ref.density = feval(@GaussianDensity);
PMBM_ref.paras.MBM.w = w;
PMBM_ref.paras.MBM.ht = ht;
PMBM_ref.paras.MBM.tt = tt;
estimates_ref = PMBM_estimator(PMBM_ref,threshold);

% Run learner solution.
PMBM = PMBMfilter();
PMBM.density = feval(@GaussianDensity);
PMBM.paras.MBM.w = w;
PMBM.paras.MBM.ht = ht;
PMBM.paras.MBM.tt = tt;
estimates = PMBM_estimator(PMBM,threshold);

% Compare.
assessVariableEqual('estimates', estimates_ref,'Feedback','Output estimates is incorrect, please modify your code and submit again.');