classdef n_objectracker
    %N_OBJECTRACKER is a class containing functions to track n object in
    %clutter. 
    %Model structures need to be called:
    %sensormodel: a structure specifies the sensor parameters
    %           P_D: object detection probability --- scalar
    %           lambda_c: average number of clutter measurements per time
    %           scan, Poisson distributed --- scalar 
    %           pdf_c: clutter (Poisson) intensity --- scalar
    %           intensity_c: clutter (Poisson) intensity --- scalar
    %motionmodel: a structure specifies the motion model parameters
    %           d: object state dimension --- scalar
    %           F: function handle return transition/Jacobian matrix
    %           f: function handle return predicted object state
    %           Q: motion noise covariance matrix
    %measmodel: a structure specifies the measurement model parameters
    %           d: measurement dimension --- scalar
    %           H: function handle return transition/Jacobian matrix
    %           h: function handle return the observation of the object
    %           state 
    %           R: measurement noise covariance matrix
    
    properties
        gating      %specify gating parameter
        reduction   %specify hypothesis reduction parameter
        density     %density class handle
    end
    
    methods
        
        function obj = initialize(obj,density_class_handle,P_G,m_d,w_min,merging_threshold,M)
            %INITIATOR initializes n_objectracker class
            %INPUT: density_class_handle: density class handle
            %       P_D: object detection probability
            %       P_G: gating size in decimal --- scalar
            %       m_d: measurement dimension --- scalar
            %       wmin: allowed minimum hypothesis weight --- scalar
            %       merging_threshold: merging threshold --- scalar
            %       M: allowed maximum number of hypotheses --- scalar
            %OUTPUT:  obj.density: density class handle
            %         obj.gating.P_G: gating size in decimal --- scalar
            %         obj.gating.size: gating size --- scalar
            %         obj.reduction.w_min: allowed minimum hypothesis
            %         weight in logarithmic scale --- scalar 
            %         obj.reduction.merging_threshold: merging threshold
            %         --- scalar 
            %         obj.reduction.M: allowed maximum number of hypotheses
            %         used in TOMHT --- scalar 
            obj.density = density_class_handle;
            obj.gating.P_G = P_G;
            obj.gating.size = chi2inv(obj.gating.P_G,m_d);
            obj.reduction.w_min = log(w_min);
            obj.reduction.merging_threshold = merging_threshold;
            obj.reduction.M = M;
        end
        
        function estimates = GNNfilter(obj, states, Z, sensormodel, motionmodel, measmodel)
            %GNNFILTER tracks n object using global nearest neighbor
            %association 
            %INPUT: obj: an instantiation of n_objectracker class
            %       states: structure array of size (1, number of objects)
            %       with two fields: 
            %                x: object initial state mean --- (object state
            %                dimension) x 1 vector 
            %                P: object initial state covariance --- (object
            %                state dimension) x (object state dimension)
            %                matrix  
            %       Z: cell array of size (total tracking time, 1), each
            %       cell stores measurements of size (measurement
            %       dimension) x (number of measurements at corresponding
            %       time step)  
            %OUTPUT:estimates: cell array of size (total tracking time, 1),
            %       each cell stores estimated object state of size (object
            %       state dimension) x (number of objects)

        end
        

        function estimates = JPDAfilter(obj, states, Z, sensormodel, motionmodel, measmodel)
            %JPDAFILTER tracks n object using joint probabilistic data
            %association
            %INPUT: obj: an instantiation of n_objectracker class
            %       states: structure array of size (1, number of objects)
            %       with two fields: 
            %                x: object initial state mean --- (object state
            %                dimension) x 1 vector 
            %                P: object initial state covariance --- (object
            %                state dimension) x (object state dimension)
            %                matrix  
            %       Z: cell array of size (total tracking time, 1), each
            %       cell stores measurements of size (measurement
            %       dimension) x (number of measurements at corresponding
            %       time step)  
            %OUTPUT:estimates: cell array of size (total tracking time, 1),
            %       each cell stores estimated object state of size (object
            %       state dimension) x (number of objects)
            
        end 
        
        
        function estimates = TOMHT(obj, states, Z, sensormodel, motionmodel, measmodel)
            %TOMHT tracks n object using track-oriented multi-hypothesis tracking
            %INPUT: obj: an instantiation of n_objectracker class
            %       states: structure array of size (1, number of objects)
            %       with two fields: 
            %                x: object initial state mean --- (object state
            %                dimension) x 1 vector 
            %                P: object initial state covariance --- (object
            %                state dimension) x (object state dimension)
            %                matrix  
            %       Z: cell array of size (total tracking time, 1), each
            %       cell stores measurements of size (measurement
            %       dimension) x (number of measurements at corresponding
            %       time step)  
            %OUTPUT:estimates: cell array of size (total tracking time, 1),
            %       each cell stores estimated object state of size (object
            %       state dimension) x (number of objects)
            
        end
        
        
    end
end






% %Choose object detection probability
% P_D = 0.9;
% %Choose clutter rate
% lambda_c = 10;
% 
% %Choose linear or nonlinear scenario
% scenario_type = 'linear';
% 
% %Creat sensor model
% range_c = [-1000 1000;-1000 1000];
% sensor_model = modelgen.sensormodel(P_D,lambda_c,range_c);
%         
% %Creat linear motion model
% T = 1;
% sigma_q = 5;
% motion_model = motionmodel.cvmodel(T,sigma_q);
%         
% %Create linear measurement model
% sigma_r = 10;
% meas_model = measmodel.cvmeasmodel(sigma_r);
%         
% %Creat ground truth model
% nbirths = 5;
% K = 100;
% tbirth = zeros(nbirths,1);
% tdeath = zeros(nbirths,1);
%         
% initial_state = repmat(struct('x',[],'P',eye(motion_model.d)),[1,nbirths]);
%         
% initial_state(1).x = [0; 0; 0; -10];        tbirth(1) = 1;   tdeath(1) = K;
% initial_state(2).x = [400; -600; -10; 5];   tbirth(2) = 1;   tdeath(2) = K;
% initial_state(3).x = [-800; -200; 20; -5];  tbirth(3) = 1;   tdeath(3) = K;
% initial_state(4).x = [0; 0; 7.5; -5];       tbirth(4) = 1;   tdeath(4) = K;
% initial_state(5).x = [-200; 800; -3; -15];  tbirth(5) = 1;   tdeath(5) = K;
% 
% %% Generate true object data (noisy or noiseless) and measurement data
% ground_truth = modelgen.groundtruth(nbirths,[initial_state.x],tbirth,tdeath,K);
% ifnoisy = 0;
% objectdata = objectdatagen(ground_truth,motion_model,ifnoisy);
% measdata = measdatagen(objectdata,sensor_model,meas_model);
% 
% %% N-object tracker parameter setting
% P_G = 0.999;            %gating size in percentage
% w_min = 1e-3;           %hypothesis pruning threshold
% merging_threshold = 2;  %hypothesis merging threshold
% M = 100;                %maximum number of hypotheses kept in TOMHT
% density_class_handle = feval(@GaussianDensity);    %density class handle
% tracker = n_objectracker();
% tracker = tracker.initialize(density_class_handle,P_G,meas_model.d,w_min,merging_threshold,M);
% 
% GNNestimates = GNNfilter(tracker, initial_state, measdata, sensor_model, motion_model, meas_model);
% JPDAestimates = JPDAfilter(tracker, initial_state, measdata, sensor_model, motion_model, meas_model);
% TOMHTestimates = TOMHT(tracker, initial_state, measdata, sensor_model, motion_model, meas_model);
%
% figure
% hold on
% grid on
% 
% for i = 1:nbirths
%     h1 = plot(cell2mat(cellfun(@(x) x(1,i), objectdata.X, 'UniformOutput', false)), ...
%         cell2mat(cellfun(@(x) x(2,i), objectdata.X, 'UniformOutput', false)), 'g', 'Linewidth', 2);
%     h2 = plot(cell2mat(cellfun(@(x) x(1,i), GNNestimates, 'UniformOutput', false)), ...
%         cell2mat(cellfun(@(x) x(2,i), GNNestimates, 'UniformOutput', false)), 'r-s', 'Linewidth', 1);
% end
% 
% xlabel('x'); ylabel('y')
% 
% xlim([-1000 1000])
% ylim([-1000 1000])
% 
% legend([h1 h2],'Ground Truth','GNN Estimates', 'Location', 'best')
% 
% set(gca,'FontSize',12) 