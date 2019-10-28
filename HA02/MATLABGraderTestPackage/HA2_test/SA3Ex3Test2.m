%Choose object detection probability
P_D = 0.5;
%Choose clutter rate
lambda_c = 30;
     
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
K = 20;
        
initial_state = repmat(struct('x',[],'P',diag([1 1 1 1*pi/90 1*pi/90].^2)),[1,nbirths]);
        
initial_state(1).x = [0; 0; 5; 0; pi/180];       tbirth(1) = 1;   tdeath(1) = K;
initial_state(2).x = [20; 20; -20; 0; pi/90];    tbirth(2) = 1;   tdeath(2) = K;
initial_state(3).x = [-20; 10; -10; 0; pi/360];  tbirth(3) = 1;   tdeath(3) = K;
initial_state(4).x = [-10; -10; 8; 0; pi/270];   tbirth(4) = 1;   tdeath(4) = K;

%% Generate true object data (noisy or noiseless) and measurement data
ground_truth = modelgen.groundtruth(nbirths,[initial_state.x],tbirth,tdeath,K);
ifnoisy = 0;
objectdata = objectdatagen(ground_truth,motion_model,ifnoisy);
measdata = measdatagen(objectdata,sensor_model,meas_model);

%% N-object tracker parameter setting
P_G = 0.999;            %gating size in percentage
w_min = 1e-3;           %hypothesis pruning threshold
merging_threshold = 2;  %hypothesis merging threshold
M = 100;                %maximum number of hypotheses kept in MHT
density_class_handle = feval(@GaussianDensity);    %density class handle
tracker = n_objectracker();
tracker = tracker.initialize(density_class_handle,P_G,meas_model.d,w_min,merging_threshold,M);

% Run learner solution.
TOMHTestimates = TOMHT(tracker, initial_state, measdata, sensor_model, motion_model, meas_model);
TOMHT_RMSE = RMSE_n_objects(objectdata.X,TOMHTestimates);

for k = 1:K
    [~,I] = sort(TOMHTestimates{k}(1,:));
    TOMHTestimates{k} = TOMHTestimates{k}(:,I);
end


tracker_ref = reference.n_objectracker();
tracker_ref = tracker_ref.initialize(density_class_handle,P_G,meas_model.d,w_min,merging_threshold,M);

% Run reference solution.
TOMHTestimates_ref = TOMHT(tracker_ref, initial_state, measdata, sensor_model, motion_model, meas_model);
TOMHT_RMSE_ref = RMSE_n_objects(objectdata.X,TOMHTestimates_ref);

for k = 1:K
    [~,I] = sort(TOMHTestimates_ref{k}(1,:));
    TOMHTestimates_ref{k} = TOMHTestimates_ref{k}(:,I);
end

% Compare.
assert(abs(TOMHT_RMSE-TOMHT_RMSE_ref) < 1e-3, 'Root Mean Square Error: Expected %f, got %f. Your TOMHT implementation is incorrect, please modify your code and submit again.', TOMHT_RMSE_ref, TOMHT_RMSE);
assessVariableEqual('TOMHTestimates', TOMHTestimates_ref,'RelativeTolerance',0.03,'Feedback','Your implementation is incorrect, please modify your code and submit again.');