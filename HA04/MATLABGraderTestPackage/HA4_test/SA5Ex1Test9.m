%Choose object detection probability
P_D = 0.9;
%Choose clutter rate
lambda_c = 10;
%Choose object survival probability
P_S = 0.99;
%Create sensor model
%Range/bearing measurement range
range_c = [-1000 1000;-pi pi];
sensor_model = modelgen.sensormodel(P_D,lambda_c,range_c);
        
%Create nonlinear motion model (coordinate turn)
T = 1;
sigmaV = 1;
sigmaOmega = pi/180;
motion_model = motionmodel.ctmodel(T,sigmaV,sigmaOmega);
        
%Create nonlinear measurement model (range/bearing)
sigma_r = 5;
sigma_b = pi/180;
s = [300;400];
meas_model = measmodel.rangebearingmeasmodel(sigma_r, sigma_b, s);
        
%Creat ground truth model
nbirths = 4;
K = 10;
tbirth = zeros(nbirths,1);
tdeath = zeros(nbirths,1);
        
initial_state(1).x = [0; 0; 5; 0; pi/180];       tbirth(1) = 1;   tdeath(1) = 50;
initial_state(2).x = [20; 20; -20; 0; pi/90];    tbirth(2) = 20;  tdeath(2) = 70;
initial_state(3).x = [-20; 10; -10; 0; pi/360];  tbirth(3) = 40;  tdeath(3) = 90;
initial_state(4).x = [-10; -10; 8; 0; pi/270];   tbirth(4) = 60;  tdeath(4) = K;

birth_model = repmat(struct('w',log(0.03),'x',[],'P',diag([1 1 1 1*pi/90 1*pi/90].^2)),[1,4]);
birth_model(1).x = [0; 0; 5; 0; pi/180];
birth_model(2).x = [20; 20; -20; 0; pi/90];
birth_model(3).x = [-20; 10; -10; 0; pi/360];
birth_model(4).x = [-10; -10; 8; 0; pi/270];

%Generate true object data (noisy or noiseless) and measurement data
ground_truth = modelgen.groundtruth(nbirths,[initial_state.x],tbirth,tdeath,K);
ifnoisy = 0;
objectdata = objectdatagen(ground_truth,motion_model,ifnoisy);
measdata = measdatagen(objectdata,sensor_model,meas_model);

%Object tracker parameter setting
P_G = 0.999;            %gating size in percentage
w_min = 1e-3;           %hypothesis pruning threshold
merging_threshold = 2;  %hypothesis merging threshold
M = 100;                %maximum number of hypotheses kept
r_min = 1e-3;           %Bernoulli component pruning threshold
r_recycle = 0.1;        %Bernoulli component recycling threshold
r_estimate = 0.4;       %Threshold used to extract estimates from Bernoullis
density_class_handle = feval(@GaussianDensity);    %density class handle
tracker = multiobjectracker();
tracker = tracker.initialize(density_class_handle,P_S,P_G,meas_model.d,w_min,merging_threshold,M,r_min,r_recycle,r_estimate);

%Create a class instance
PMBM_ref = reference.PMBMfilter();
PMBM_ref = PMBM_ref.initialize(tracker.density,birth_model);

%Create a class instance
PMBM = PMBMfilter();
PMBM = PMBM.initialize(tracker.density,birth_model);

for k = 1:K
    % Run learner solution.
%     PMBM = PMBM_ref;
    PMBM = PMBM_update(PMBM,measdata{k},meas_model,sensor_model,tracker.gating,tracker.reduction.w_min,tracker.reduction.M);
    % Run reference solution.
    %PMBM update
    PMBM_ref = PMBM_ref.PMBM_update(measdata{k},meas_model,sensor_model,tracker.gating,tracker.reduction.w_min,tracker.reduction.M);
    
    PPP_ref.w = PMBM_ref.paras.PPP.w;
    [PPP_ref_w,I] = sort(PPP_ref.w);
    PPP_ref.states = PMBM_ref.paras.PPP.states;
    PPP_ref_states = PPP_ref.states(I);
    
    MBM.w_ref = PMBM_ref.paras.MBM.w;
    [MBM_w_ref,Iw] = sort(MBM.w_ref);
    
    tt_ref = cellfun(@(x) x(1).state.x(1), PMBM_ref.paras.MBM.tt);
    [~,I] = sort(tt_ref);
    temp = PMBM_ref.paras.MBM.tt(I);
    MBM_tt_ref = cell(length(temp),1);
    for i = 1:length(temp)
        state_ref = arrayfun(@(x) x.state.x(1), temp{i});
        [~,I] = sort(state_ref);
        MBM_tt_ref{i} = temp{i}(I);
    end
    
    PPP.w = PMBM.paras.PPP.w;
    [PPP_w,I] = sort(PPP.w);
    PPP.states = PMBM.paras.PPP.states;
    PPP_states = PPP.states(I);
    
    MBM.w = PMBM.paras.MBM.w;
    [MBM_w,Iw] = sort(MBM.w);
    
    tt = cellfun(@(x) x(1).state.x(1), PMBM.paras.MBM.tt);
    [~,I] = sort(tt);
    temp = PMBM.paras.MBM.tt(I);
    MBM_tt = cell(length(temp),1);
    for i = 1:length(temp)
        state = arrayfun(@(x) x.state.x(1), temp{i});
        [~,I] = sort(state);
        MBM_tt{i} = temp{i}(I);
    end
    
    %PMBM prediction
    PMBM_ref = PMBM_ref.PMBM_predict(P_S,motion_model,birth_model);
    PMBM = PMBM.PMBM_predict(P_S,motion_model,birth_model);
    
    % Compare.
    assessVariableEqual('PPP_ref_w', PPP_w,'RelativeTolerance',0.0001,'Feedback','Parameter PPP.w is incorrect, please modify your code and submit again.');
    assessVariableEqual('PPP_ref_states', PPP_states,'RelativeTolerance',0.001,'Feedback','Parameter PPP.states is incorrect, please modify your code and submit again.');
    assessVariableEqual('MBM_w_ref', MBM_w,'RelativeTolerance',0.0001,'Feedback','Parameter MBM.w is incorrect, please modify your code and submit again.');
    assessVariableEqual('MBM_tt_ref', MBM_tt,'RelativeTolerance',0.001,'Feedback','Parameter MBM.tt is incorrect, please modify your code and submit again.');
end


