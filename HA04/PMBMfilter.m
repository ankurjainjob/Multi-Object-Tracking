classdef PMBMfilter
    %PMBMFILTER is a class containing necessary functions to implement the
    %PMBM filter
    %Model structures need to be called:
    %    sensormodel: a structure specifies the sensor parameters
    %           P_D: object detection probability --- scalar
    %           lambda_c: average number of clutter measurements per time scan, 
    %                     Poisson distributed --- scalar
    %           pdf_c: value of clutter pdf --- scalar
    %           intensity_c: Poisson clutter intensity --- scalar
    %       motionmodel: a structure specifies the motion model parameters
    %           d: object state dimension --- scalar
    %           F: function handle return transition/Jacobian matrix
    %           f: function handle return predicted object state
    %           Q: motion noise covariance matrix
    %       measmodel: a structure specifies the measurement model parameters
    %           d: measurement dimension --- scalar
    %           H: function handle return transition/Jacobian matrix
    %           h: function handle return the observation of the target state
    %           R: measurement noise covariance matrix
    %       birthmodel: a structure array specifies the birth model (Gaussian
    %       mixture density) parameters --- (1 x number of birth components)
    %           w: weights of mixture components (in logarithm domain)
    %           x: mean of mixture components
    %           P: covariance of mixture components
    properties
        density %density class handle
        paras   %%parameters specify a PMBM
    end
    
    methods
        function obj = initialize(obj,density_class_handle,birthmodel)
            %INITIATOR initializes PMBMfilter class
            %INPUT: density_class_handle: density class handle
            %       birthmodel: a struct specifying the intensity (mixture)
            %       of a PPP birth model
            %OUTPUT:obj.density: density class handle
            %       obj.paras.PPP.w: weights of mixture components in PPP
            %       intensity --- vector of size (number of mixture
            %       components x 1) in logarithmic scale
            %       obj.paras.PPP.states: parameters of mixture components
            %       in PPP intensity struct array of size (number of
            %       mixture components x 1)
            %       obj.paras.MBM.w: weights of MBs --- vector of size
            %       (number of MBs (global hypotheses) x 1) in logarithmic 
            %       scale
            %       obj.paras.MBM.ht: hypothesis table --- matrix of size
            %       (number of global hypotheses x number of hypothesis
            %       trees). Entry (h,i) indicates that the (h,i)th local
            %       hypothesis in the ith hypothesis tree is included in
            %       the hth global hypothesis. If entry (h,i) is zero, then
            %       no local hypothesis from the ith hypothesis tree is
            %       included in the hth global hypothesis.
            %       obj.paras.MBM.tt: local hypotheses --- cell of size
            %       (number of hypothesis trees x 1). The ith cell contains
            %       local hypotheses in struct form of size (number of
            %       local hypotheses in the ith hypothesis tree x 1). Each
            %       struct has two fields: r: probability of existence;
            %       state: parameters specifying the object density
            
            obj.density = density_class_handle;
            obj.paras.PPP.w = [birthmodel.w]';
            obj.paras.PPP.states = rmfield(birthmodel,'w')';
            obj.paras.MBM.w = [];
            obj.paras.MBM.ht = [];
            obj.paras.MBM.tt = {};
        end
        
        function Bern = Bern_predict(obj,Bern,motionmodel,P_S)
            %BERN_PREDICT performs prediction step for a Bernoulli component
            %INPUT: Bern: a struct that specifies a Bernoulli component,
            %             with fields: r: probability of existence ---
            %                          scalar;
            %                          state: a struct contains parameters
            %                          describing the object pdf
            %       P_S: object survival probability
            
        end
        
        function [Bern, lik_undetected] = Bern_undetected_update(obj,tt_entry,P_D)
            %BERN_UNDETECTED_UPDATE calculates the likelihood of missed
            %detection, and creates new local hypotheses due to missed
            %detection.
            %INPUT: tt_entry: a (2 x 1) array that specifies the index of
            %       local hypotheses. (i,j) indicates the jth local
            %       hypothesis in the ith hypothesis tree. 
            %       P_D: object detection probability --- scalar
            %OUTPUT:Bern: a struct that specifies a Bernoulli component,
            %       with fields: r: probability of existence --- scalar;
            %                    state: a struct contains parameters
            %                    describing the object pdf
            %       lik_undetected: missed detection likelihood --- scalar
            %       in logorithmic scale
            
        end
        
        function lik_detected = Bern_detected_update_lik(obj,tt_entry,z,measmodel,P_D)
            %BERN_DETECTED_UPDATE_LIK calculates the predicted likelihood
            %for a given local hypothesis. 
            %INPUT: tt_entry: a (2 x 1) array that specifies the index of
            %       local hypotheses. (i,j) indicates the jth
            %       local hypothesis in the ith hypothesis tree.
            %       z: measurement array --- (measurement dimension x
            %       number of measurements)
            %       P_D: object detection probability --- scalar
            %OUTPUT:lik_detected: predicted likelihood --- (number of
            %measurements x 1) array in logarithmic scale 
            
        end
        
        function Bern = Bern_detected_update_state(obj,tt_entry,z,measmodel)
            %BERN_DETECTED_UPDATE_STATE creates the new local hypothesis
            %due to measurement update. 
            %INPUT: tt_entry: a (2 x 1) array that specifies the index of
            %                 local hypotheses. (i,j) indicates the jth
            %                 local hypothesis in the ith hypothesis tree.
            %       z: measurement vector --- (measurement dimension x 1)
            %OUTPUT:Bern: a struct that specifies a Bernoulli component,
            %             with fields: r: probability of existence ---
            %                          scalar; 
            %                          state: a struct contains parameters
            %                          describing the object pdf 
            
        end
        
        function obj = PPP_predict(obj,motionmodel,birthmodel,P_S)
            %PPP_PREDICT performs predicion step for PPP components
            %hypothesising undetected objects.
            %INPUT: P_S: object survival probability --- scalar          

        end
        
        function [Bern, lik_new] = PPP_detected_update(obj,indices,z,measmodel,P_D,clutter_intensity)
            %PPP_DETECTED_UPDATE creates a new local hypothesis by
            %updating the PPP with measurement and calculates the
            %corresponding likelihood.
            %INPUT: z: measurement vector --- (measurement dimension x 1)
            %       P_D: object detection probability --- scalar
            %       clutter_intensity: Poisson clutter intensity --- scalar
            %       indices: boolean vector, if measurement z is inside the
            %       gate of mixture component i, then indices(i) = true
            %OUTPUT:Bern: a struct that specifies a Bernoulli component,
            %             with fields: r: probability of existence ---
            %             scalar;
            %             state: a struct contains parameters describing
            %             the object pdf
            %       lik_new: predicted likelihood of PPP --- scalar in
            %       logarithmic scale 

        end
        
        function obj = PPP_undetected_update(obj,P_D)
            %PPP_UNDETECTED_UPDATE performs PPP update for missed detection.
            %INPUT: P_D: object detection probability --- scalar
            
        end
        
        function obj = PPP_reduction(obj,prune_threshold,merging_threshold)
            %PPP_REDUCTION truncates mixture components in the PPP
            %intensity by pruning and merging
            %INPUT: prune_threshold: pruning threshold --- scalar in
            %       logarithmic scale
            %       merging_threshold: merging threshold --- scalar
            [obj.paras.PPP.w, obj.paras.PPP.states] = hypothesisReduction.prune(obj.paras.PPP.w, obj.paras.PPP.states, prune_threshold);
            if ~isempty(obj.paras.PPP.w)
                [obj.paras.PPP.w, obj.paras.PPP.states] = hypothesisReduction.merge(obj.paras.PPP.w, obj.paras.PPP.states, merging_threshold, obj.density);
            end
        end
        
        function obj = Bern_recycle(obj,prune_threshold,recycle_threshold)
            %BERN_RECYCLE recycles Bernoulli components with small
            %probability of existence, adds them to the PPP component, and
            %re-index the hypothesis table. If a hypothesis tree contains no
            %local hypothesis after pruning, this tree is removed. After
            %recycling, merge similar Gaussian components in the PPP
            %intensity
            %INPUT: prune_threshold: Bernoulli components with probability
            %       of existence smaller than this threshold are pruned ---
            %       scalar
            %       recycle_threshold: Bernoulli components with probability
            %       of existence smaller than this threshold needed to be
            %       recycled --- scalar
            
            n_tt = length(obj.paras.MBM.tt);
            for i = 1:n_tt
                idx = arrayfun(@(x) x.r<recycle_threshold & x.r>=prune_threshold, obj.paras.MBM.tt{i});
                if any(idx)
                    %Here, we should also consider the weights of different MBs
                    idx_t = find(idx);
                    n_h = length(idx_t);
                    w_h = zeros(n_h,1);
                    for j = 1:n_h
                        idx_h = obj.paras.MBM.ht(:,i) == idx_t(j);
                        [~,w_h(j)] = normalizeLogWeights(obj.paras.MBM.w(idx_h));
                    end
                    %Recycle
                    temp = obj.paras.MBM.tt{i}(idx);
                    obj.paras.PPP.w = [obj.paras.PPP.w;log([temp.r]')+w_h];
                    obj.paras.PPP.states = [obj.paras.PPP.states;[temp.state]'];
                end
                idx = arrayfun(@(x) x.r<recycle_threshold, obj.paras.MBM.tt{i});
                if any(idx)
                    %Remove Bernoullis
                    obj.paras.MBM.tt{i} = obj.paras.MBM.tt{i}(~idx);
                    %Update hypothesis table, if a Bernoulli component is
                    %pruned, set its corresponding entry to zero
                    idx = find(idx);
                    for j = 1:length(idx)
                        temp = obj.paras.MBM.ht(:,i);
                        temp(temp==idx(j)) = 0;
                        obj.paras.MBM.ht(:,i) = temp;
                    end
                end
            end
            
            %Remove tracks that contains no valid local hypotheses
            idx = sum(obj.paras.MBM.ht,1)~=0;
            obj.paras.MBM.ht = obj.paras.MBM.ht(:,idx);
            obj.paras.MBM.tt = obj.paras.MBM.tt(idx);
            if isempty(obj.paras.MBM.ht)
                %Ensure the algorithm still works when all Bernoullis are
                %recycled
                obj.paras.MBM.w = [];
            end
            
            %Re-index hypothesis table
            n_tt = length(obj.paras.MBM.tt);
            for i = 1:n_tt
                idx = obj.paras.MBM.ht(:,i) > 0;
                [~,~,obj.paras.MBM.ht(idx,i)] = unique(obj.paras.MBM.ht(idx,i),'rows','stable');
            end
            
            %Merge duplicate hypothesis table rows
            if ~isempty(obj.paras.MBM.ht)
                [ht,~,ic] = unique(obj.paras.MBM.ht,'rows','stable');
                if(size(ht,1)~=size(obj.paras.MBM.ht,1))
                    %There are duplicate entries
                    w = zeros(size(ht,1),1);
                    for i = 1:size(ht,1)
                        indices_dupli = (ic==i);
                        [~,w(i)] = normalizeLogWeights(obj.paras.MBM.w(indices_dupli));
                    end
                    obj.paras.MBM.ht = ht;
                    obj.paras.MBM.w = w;
                end
            end
            
        end
        
        function obj = PMBM_predict(obj,P_S,motionmodel,birthmodel)
            %PMBM_PREDICT performs PMBM prediction step.

        end
        
        function obj = PMBM_update(obj,z,measmodel,sensormodel,gating,w_min,M)
            %PMBM_UPDATE performs PMBM update step.
            %INPUT: z: measurements --- array of size (measurement
            %       dimension x number of measurements)
            %       gating: a struct with two fields that specifies gating
            %       parameters: P_G: gating size in decimal --- scalar;
            %                   size: gating size --- scalar.
            %       wmin: hypothesis weight pruning threshold --- scalar in
            %       logarithmic scale
            %       M: maximum global hypotheses kept
            
        end
        
        function estimates = PMBM_estimator(obj,threshold)
            %PMBM_ESTIMATOR performs object state estimation in the PMBM
            %filter
            %INPUT: threshold (if exist): object states are extracted from
            %       Bernoulli components with probability of existence no
            %       less than this threhold in Estimator 1. Given the
            %       probabilities of detection and survival, this threshold
            %       determines the number of consecutive misdetections
            %OUTPUT:estimates: estimated object states in matrix form of
            %       size (object state dimension) x (number of objects)
            %%%
            %First, select the multi-Bernoulli with the highest weight.
            %Second, report the mean of the Bernoulli components whose
            %existence probability is above a threshold. 
            
        end
    
    end
end





% 
% %A simple sanity check is to compare the GOSPA error of the PMBM filter and the PHD filter. 
% 
% %Choose object detection probability
% P_D = 0.9;
% %Choose clutter rate
% lambda_c = 10;
% %Choose object survival probability
% P_S = 0.99;
% 
% %Choose linear or nonlinear scenario
% scenario_type = 'linear';
% 
% %% Create tracking scenario
% switch(scenario_type)
%     case 'linear'
%         %Create sensor model
%         range_c = [-1000 1000;-1000 1000];
%         sensor_model = modelgen.sensormodel(P_D,lambda_c,range_c);
%         
%         %Create linear motion model
%         T = 1;
%         sigma_q = 5;
%         motion_model = motionmodel.cvmodel(T,sigma_q);
%         
%         %Create linear measurement model
%         sigma_r = 10;
%         meas_model = measmodel.cvmeasmodel(sigma_r);
%         
%         %Create ground truth model
%         nbirths = 12;
%         K = 100;
%         tbirth = zeros(nbirths,1);
%         tdeath = zeros(nbirths,1);
%         
%         initial_state(1).x  = [ 0; 0; 0; -10 ];            tbirth(1)  = 1;     tdeath(1)  = 70;
%         initial_state(2).x  = [ 400; -600; -10; 5 ];       tbirth(2)  = 1;     tdeath(2)  = K+1;
%         initial_state(3).x  = [ -800; -200; 20; -5 ];      tbirth(3)  = 1;     tdeath(3)  = 70;
%         initial_state(4).x  = [ 400; -600; -7; -4 ];       tbirth(4)  = 20;    tdeath(4)  = K+1;
%         initial_state(5).x  = [ 400; -600; -2.5; 10 ];     tbirth(5)  = 20;    tdeath(5)  = K+1;
%         initial_state(6).x  = [ 0; 0; 7.5; -5 ];           tbirth(6)  = 20;    tdeath(6)  = K+1;
%         initial_state(7).x  = [ -800; -200; 12; 7 ];       tbirth(7)  = 40;    tdeath(7)  = K+1;
%         initial_state(8).x  = [ -200; 800; 15; -10 ];      tbirth(8)  = 40;    tdeath(8)  = K+1;
%         initial_state(9).x  = [ -800; -200; 3; 15 ];       tbirth(9)   = 60;   tdeath(9)  = K+1;
%         initial_state(10).x  = [ -200; 800; -3; -15 ];     tbirth(10)  = 60;   tdeath(10) = K+1;
%         initial_state(11).x  = [ 0; 0; -20; -15 ];         tbirth(11)  = 80;   tdeath(11) = K+1;
%         initial_state(12).x  = [ -200; 800; 15; -5 ];      tbirth(12)  = 80;   tdeath(12) = K+1;
%         
%         birth_model = repmat(struct('w',log(0.03),'x',[],'P',400*eye(motion_model.d)),[1,4]);
%         birth_model(1).x = [ 0; 0; 0; 0];
%         birth_model(2).x = [ 400; -600; 0; 0];
%         birth_model(3).x = [ -800; -200; 0; 0 ];
%         birth_model(4).x = [ -200; 800; 0; 0 ];
%         
%     case 'nonlinear'
%         %Create sensor model
%         %Range/bearing measurement range
%         range_c = [-1000 1000;-pi pi];
%         sensor_model = modelgen.sensormodel(P_D,lambda_c,range_c);
%         
%         %Create nonlinear motion model (coordinate turn)
%         T = 1;
%         sigmaV = 1;
%         sigmaOmega = pi/180;
%         motion_model = motionmodel.ctmodel(T,sigmaV,sigmaOmega);
%         
%         %Create nonlinear measurement model (range/bearing)
%         sigma_r = 5;
%         sigma_b = pi/180;
%         s = [300;400];
%         meas_model = measmodel.rangebearingmeasmodel(sigma_r, sigma_b, s);
%         
%         %Creat ground truth model
%         nbirths = 4;
%         K = 100;
%         tbirth = zeros(nbirths,1);
%         tdeath = zeros(nbirths,1);
%         
%         initial_state(1).x = [0; 0; 5; 0; pi/180];       tbirth(1) = 1;   tdeath(1) = 50;
%         initial_state(2).x = [20; 20; -20; 0; pi/90];    tbirth(2) = 20;  tdeath(2) = 70;
%         initial_state(3).x = [-20; 10; -10; 0; pi/360];  tbirth(3) = 40;  tdeath(3) = 90;
%         initial_state(4).x = [-10; -10; 8; 0; pi/270];   tbirth(4) = 60;  tdeath(4) = K;
% 
%         birth_model = repmat(struct('w',log(0.03),'x',[],'P',diag([1 1 1 1*pi/90 1*pi/90].^2)),[1,4]);
%         birth_model(1).x = [0; 0; 5; 0; pi/180];
%         birth_model(2).x = [20; 20; -20; 0; pi/90];
%         birth_model(3).x = [-20; 10; -10; 0; pi/360];
%         birth_model(4).x = [-10; -10; 8; 0; pi/270];
% end
% 
% %% Generate true object data (noisy or noiseless) and measurement data
% ground_truth = modelgen.groundtruth(nbirths,[initial_state.x],tbirth,tdeath,K);
% ifnoisy = 1;
% objectdata = objectdatagen(ground_truth,motion_model,ifnoisy);
% measdata = measdatagen(objectdata,sensor_model,meas_model);
% 
% %% Object tracker parameter setting
% P_G = 0.999;            %gating size in percentage
% w_min = 1e-3;           %hypothesis pruning threshold
% merging_threshold = 2;  %hypothesis merging threshold
% M = 100;                %maximum number of hypotheses kept
% r_min = 1e-3;           %Bernoulli component pruning threshold
% r_recycle = 0.05;       %Bernoulli component recycling threshold
% r_estimate = 0.4;       %Threshold used to extract estimates from Bernoullis
% density_class_handle = feval(@GaussianDensity);    %density class handle
% 
% %Please check the provided script multiobjectracker_HA4.m for details of the function multiobjectracker
% tracker = multiobjectracker();
% tracker = tracker.initialize(density_class_handle,P_S,P_G,meas_model.d,w_min,merging_threshold,M,r_min,r_recycle,r_estimate);
% 
% %% GM-PHD filter
% GMPHDestimates = GMPHDfilter(tracker, birth_model, measdata, sensor_model, motion_model, meas_model);
% 
% %% PMBM filter
% PMBMestimates = PMBMfilter(tracker, birth_model, measdata, sensor_model, motion_model, meas_model);
% 
% %% Ploting
% true_state = cell2mat(objectdata.X');
% GMPHD_estimated_state = cell2mat(GMPHDestimates');
% PMBM_estimated_state = cell2mat(PMBMestimates');
% 
% true_cardinality = objectdata.N;
% GMPHD_estimated_cardinality = cellfun(@(x) size(x,2), GMPHDestimates);
% PMBM_estimated_cardinality = cellfun(@(x) size(x,2), PMBMestimates);
% 
% c = 100;
% p = 1;
% GMPHD_d_gospa = zeros(K,1);
% PMBM_d_gospa = zeros(K,1);
% for k = 1:K
%     if isempty(GMPHDestimates{k}); GMPHDestimates{k} = zeros(motion_model.d,0); end
%     if isempty(PMBMestimates{k}); PMBMestimates{k} = zeros(motion_model.d,0); end
%     %Evaluate kinematics estimation performance using GOSPA metric
%     [GMPHD_d_gospa(k), ~, GMPHD_decomposed_cost(k)] = GOSPA(objectdata.X{k}, GMPHDestimates{k}, p, c, 2);
%     [PMBM_d_gospa(k), ~, PMBM_decomposed_cost(k)] = GOSPA(objectdata.X{k}, PMBMestimates{k}, p, c, 2);
% end
% 
% %Trajectory plot
% figure
% hold on
% grid on
% 
% h1 = plot(true_state(1,:), true_state(2,:), 'bo');
% h2 = plot(PMBM_estimated_state(1,:), PMBM_estimated_state(2,:),'r+');
% 
% xlabel('x (m)'); ylabel('y (m)')
% legend([h1 h2],'Ground Truth','PMBM', 'Location', 'best')
% set(gca,'FontSize',12) 
% 
% %Cardinality plot
% figure
% grid on
% hold on
% 
% h1 = plot(1:length(true_cardinality),true_cardinality,'bo','linewidth',2);
% h2 = plot(1:length(GMPHD_estimated_cardinality),GMPHD_estimated_cardinality,'m+','linewidth',2);
% h3 = plot(1:length(PMBM_estimated_cardinality),PMBM_estimated_cardinality,'r+','linewidth',2);
% 
% xlabel('Time step')
% ylabel('Cardinality')
% legend([h1 h2 h3],'Ground Truth','PHD','PMBM', 'Location', 'best')
% set(gca,'FontSize',12) 
% 
% %GOSPA plot
% figure
% subplot(4,1,1)
% grid on
% hold on
% plot(1:K,GMPHD_d_gospa,'linewidth',2)
% plot(1:K,PMBM_d_gospa,'linewidth',2)
% legend('PHD','PMBM','location','best')
% ylabel('GOSPA')
% set(gca,'FontSize',12) 
% 
% subplot(4,1,2)
% grid on
% hold on
% plot(1:K,[GMPHD_decomposed_cost.localisation],'linewidth',2)
% plot(1:K,[PMBM_decomposed_cost.localisation],'linewidth',2)
% ylabel('Kinematics')
% set(gca,'FontSize',12) 
% 
% subplot(4,1,3)
% grid on
% hold on
% plot(1:K,[GMPHD_decomposed_cost.missed],'linewidth',2)
% plot(1:K,[PMBM_decomposed_cost.missed],'linewidth',2)
% ylabel('Missed')
% set(gca,'FontSize',12) 
% 
% subplot(4,1,4)
% grid on
% hold on
% plot(1:K,[GMPHD_decomposed_cost.false],'linewidth',2)
% plot(1:K,[PMBM_decomposed_cost.false],'linewidth',2)
% xlabel('Time Step'); ylabel('False')
% set(gca,'FontSize',12) 