%Choose object detection probability
P_D = 0.6;
%Choose clutter rate
lambda_c = 10 *5;
%Create sensor model
range_c = [-1000 1000;-1000 1000];
sensor_model = modelgen.sensormodel(P_D,lambda_c,range_c);
        
%Create ground truth model
nbirths = 1;
K = 100;
initial_state.x = [0; 0; 10; 10];
ground_truth = modelgen.groundtruth(nbirths,initial_state.x,1,K+1,K);
   
% initial guess
initial_state.x = [10; 10; 0; 0];
initial_state.P = eye(4) * 300^2;


%Create linear motion model
T = 1;
sigma_q = 5;
motion_model = motionmodel.cvmodel(T,sigma_q);
        
%Create linear measurement model
sigma_r = 10;
meas_model = measmodel.cvmeasmodel(sigma_r);

%Generate true object data (noisy or noiseless) and measurement data
ifnoisy = 0;
objectdata = objectdatagen(ground_truth,motion_model,ifnoisy);
measdata = measdatagen(objectdata,sensor_model,meas_model);

%Single object tracker parameter setting
P_G = 0.999;            %gating size in percentage
w_min = 1e-3;           %hypothesis pruning threshold
merging_threshold = 2;  %hypothesis merging threshold
M = 100;                %maximum number of hypotheses kept in Gaussian sum filter
density_class_handle = feval(@GaussianDensity);    %density class handle
tracker = singleobjectracker();
tracker = tracker.initialize(density_class_handle,P_G,meas_model.d,w_min,merging_threshold,M);

%Nearest neighbour filter
[x_NN, P_NN] = nearestNeighbourFilter(tracker, initial_state, measdata, sensor_model, motion_model, meas_model);

%Probabilistic data association filter
[x_PDA, P_PDA] = probDataAssocFilter(tracker, initial_state, measdata, sensor_model, motion_model, meas_model);

%Gaussian sum filter
[x_GSF, P_GSF] = GaussianSumFilter(tracker, initial_state, measdata, sensor_model, motion_model, meas_model);


est = struct('x',x_GSF, 'P', P_GSF)
animate = Animate_2D_tracking();
animate.animate(est, initial_state, measdata, meas_model, range_c);



true_state = cell2mat(objectdata.X');
NN_estimated_state = cell2mat(x_NN');
PDA_estimated_state = cell2mat(x_PDA');
GS_estimated_state = cell2mat(x_GSF');

figure
hold on
grid on

plot(true_state(1,:), true_state(2,:), 'g','Linewidth', 2)
plot(NN_estimated_state(1,:), NN_estimated_state(2,:), 'r-s' , 'Linewidth', 1)
plot(PDA_estimated_state(1,:), PDA_estimated_state(2,:), 'm-o' , 'Linewidth', 1)
plot(GS_estimated_state(1,:), GS_estimated_state(2,:), 'b-d' , 'Linewidth', 1)

xlabel('x (m)')
ylabel('y (m)')
legend('Ground Truth','Nearest Neighbour', 'Probalistic Data Association', 'Gaussian Sum', 'Location', 'best')

set(gca,'FontSize',12)