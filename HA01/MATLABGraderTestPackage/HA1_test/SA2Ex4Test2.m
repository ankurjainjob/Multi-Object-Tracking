%Choose object detection probability
P_D = 0.7;
%Choose clutter rate
lambda_c = 60;

%Create sensor model
%Range/bearing measurement range
range_c = [-1000 1000;-pi pi];
sensor_model = modelgen.sensormodel(P_D,lambda_c,range_c);
        
%Creat ground truth model
nbirths = 1;
K = 50;
initial_state.x = [0; 0; 10; 0; pi/180];
initial_state.P = diag([1 1 1 1*pi/180 1*pi/180].^2);
ground_truth = modelgen.groundtruth(nbirths,initial_state.x,1,K+1,K);
        
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

%% Generate true object data (noisy or noiseless) and measurement data
ifnoisy = 0;
objectdata = objectdatagen(ground_truth,motion_model,ifnoisy);
measdata = measdatagen(objectdata,sensor_model,meas_model);

%% Single object tracker parameter setting
P_G = 0.999;            %gating size in percentage
wmin = 1e-3;            %hypothesis pruning threshold
merging_threshold = 4;  %hypothesis merging threshold
M = 50;                %maximum number of hypotheses kept in Gaussian sum filter
density_class_handle = feval(@GaussianDensity);    %density class handle


% Run learner solution.
tracker = singleobjectracker();
tracker = tracker.initialize(density_class_handle,P_G,meas_model.d,wmin,merging_threshold,M);

probDataAssocEstimates = probDataAssocFilter(tracker, initial_state, measdata, sensor_model, motion_model, meas_model);
probDataAssocRMSE = RMSE(probDataAssocEstimates,objectdata.X);

% Run reference solution.
tracker_ref = reference.singleobjectracker();
tracker_ref = tracker_ref.initialize(density_class_handle,P_G,meas_model.d,wmin,merging_threshold,M);

probDataAssocEstimates_ref = probDataAssocFilter(tracker_ref, initial_state, measdata, sensor_model, motion_model, meas_model);
probDataAssocRMSE_ref = RMSE(probDataAssocEstimates_ref,objectdata.X);

% Compare.
assert(abs(probDataAssocRMSE-probDataAssocRMSE_ref) < 1e-3, 'Root Mean Square Error: Expected %f, got %f.', probDataAssocRMSE_ref, probDataAssocRMSE);
assessVariableEqual('probDataAssocEstimates', probDataAssocEstimates_ref,'RelativeTolerance',0.03,'Feedback','Your Probabilistic Data Association Filter implementation is incorrect, please modify your code and submit again.');