classdef PMBMfilter
    %PMBMFILTER is a class containing necessary functions to implement the
    %PMBM tracker for the set of all trajectories, including both "alive"
    %and "dead" trajectories
    
    properties
        density %density class handle
        paras   %parameters specify a PMBM
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
            %       components x 1)
            %       obj.paras.PPP.states: parameters of mixture components
            %       in PPP intensity struct array of size (number of
            %       mixture components x 1)
            %       obj.paras.MBM.w: weights of MBs --- vector of size
            %       (number of MBs (global hypotheses) x 1)
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
            %       states: parameters specifying the object density
            
            obj.density = density_class_handle;
            obj.paras.PPP.w = log([birthmodel.w]');
            obj.paras.PPP.states = rmfield(birthmodel,'w')';
            obj.paras.MBM.w = [];
            obj.paras.MBM.ht = [];
            obj.paras.MBM.tt = {};
        end
        
        function Bern = Bern_predict(obj,Bern,motionmodel,P_S,r_min)
            %BERN_PREDICT performs prediction step for a Bernoulli component
            %INPUT: Bern: a struct that specifies a Bernoulli component,
            %             with fields: r: probability of existence --- scalar
            %             state: a struct contains parameters describing
            %             the object pdf 
            %             t_birth: birth time --- scalar
            %             t_death: predicted time of death --- vector
            %             w_death: weight of predicted time of death ---
            %             vector
            %       P_S: object survival probability --- scalar
            %       r_min: allowed minimum object's probability of
            %       existence --- scalar
            
            if Bern.w_death(end) >= r_min
                Bern.state = obj.density.predict(Bern.state,motionmodel);
                %Predicted time of death
                Bern.t_death = [Bern.t_death Bern.t_death(end)+1];
                Bern.w_death = [Bern.w_death(1:end-1) Bern.w_death(end)*(1-P_S) Bern.w_death(end)*P_S];
            end
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
            %                    t_birth: birth time --- scalar
            %                    t_death: predicted time of death ---
            %                    vector 
            %                    w_death: weight of predicted time of death
            %                    --- vector
            %       lik_undetected: missed detection likelihood --- scalar
            %       in logorithmic scale
            
            Bern = obj.paras.MBM.tt{tt_entry(1)}(tt_entry(2));
            l_nodetect = Bern.r*(1 - P_D*Bern.w_death(end));
            lik_undetected = 1 - Bern.r + l_nodetect;
            Bern.r = l_nodetect/lik_undetected;
            lik_undetected = log(lik_undetected);
            
            %Updated time of death
            Bern.w_death = [Bern.w_death(1:end-1) Bern.w_death(end)*(1-P_D)]/(1-Bern.w_death(end)*P_D);
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
            
            Bern = obj.paras.MBM.tt{tt_entry(1)}(tt_entry(2));
            lik_detected = obj.density.predictedLikelihood(Bern.state,z,measmodel) + log(P_D*Bern.r*Bern.w_death(end));
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
            %             state: a struct contains parameters describing
            %             the object pdf  
            %             t_birth: birth time --- scalar
            %             t_death: predicted time of death --- vector
            %             w_death: weight of predicted time of death ---
            %             vector
            
            Bern = obj.paras.MBM.tt{tt_entry(1)}(tt_entry(2));
            Bern.state = obj.density.update(Bern.state,z,measmodel);
            Bern.r = 1;
            
            %Updated time of death
            Bern.t_death = Bern.t_death(end);
            Bern.w_death = 1;
        end
        
        function obj = PPP_predict(obj,motionmodel,birthmodel,P_S)
            %PPP_PREDICT performs predicion step for PPP components
            %hypothesising undetected objects.
            %INPUT: P_S: object survival probability --- scalar
            
            %Predict existing PPP intensity
            obj.paras.PPP.w = obj.paras.PPP.w + log(P_S);
            obj.paras.PPP.states = arrayfun(@(x) obj.density.predict_state(x,motionmodel), obj.paras.PPP.states);
            %Incorporate PPP birth intensity into PPP intensity
            obj.paras.PPP.w = [obj.paras.PPP.w;log([birthmodel.w]')];
            obj.paras.PPP.states = [obj.paras.PPP.states;rmfield(birthmodel,'w')'];
        end
        
        function [Bern, lik_new] = PPP_detected_update(obj,k,indices,z,measmodel,P_D,clutter_intensity)
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
            %             t_birth: birth time --- scalar
            %             t_death: predicted time of death --- vector
            %             w_death: weight of predicted time of death ---
            %             vector
            %       lik_new: predicted likelihood of PPP --- scalar in
            %       logarithmic scale 
            
            state_upd = arrayfun(@(x) obj.density.update(x,z,measmodel), obj.paras.PPP.states(indices));
            w_upd = arrayfun(@(x) obj.density.predictedLikelihood(x,z,measmodel), obj.paras.PPP.states(indices)) + obj.paras.PPP.w(indices) + log(P_D);
            
            [w_upd,log_numerator] = normalizeLogWeights(w_upd);
            [~,lik_new] = normalizeLogWeights([log_numerator log(clutter_intensity)]);
            Bern.r = exp(log_numerator-lik_new);
            Bern.state = obj.density.momentMatching(w_upd,state_upd);
            
            Bern.t_birth = k;
            Bern.t_death = k;
            Bern.w_death = 1;
            
        end
        
        function obj = PPP_undetected_update(obj,P_D)
            %PPP_UNDETECTED_UPDATE performs PPP update for missed detection.
            %INPUT: P_D: object detection probability --- scalar
            
            obj.paras.PPP.w = obj.paras.PPP.w + log(1 - P_D);
        end
        
        function obj = PPP_prune(obj,threshold)
            %PPP_PRUNE prunes mixture components in the PPP intensity with
            %small weight. 
            %INPUT: threshold: pruning threshold --- scalar in logarithmic
            %scale 
            
            [obj.paras.PPP.w, obj.paras.PPP.states] = hypothesisReduction.prune(obj.paras.PPP.w, obj.paras.PPP.states, threshold);
        end
        
        function obj = Bern_prune(obj,prune_threshold)
            %BERN_PRUNE removes Bernoulli components with small probability
            %of existence and re-index the hypothesis table. If a track
            %contains no single object hypothesis after pruning, this track
            %is removed.
            %INPUT: prune_threshold: Bernoulli components with probability
            %of existence smaller than this threshold will be pruned --- scalar
            
            n_tt = length(obj.paras.MBM.tt);
            for i = 1:n_tt
                %Find all Bernoulli components needed to be pruned
                idx = arrayfun(@(x) x.r<prune_threshold, obj.paras.MBM.tt{i});
                %Prune theses Bernoulli components
                obj.paras.MBM.tt{i} = obj.paras.MBM.tt{i}(~idx);
                idx = find(idx);
                %Update hypothesis table, if a Bernoulli component is
                %pruned, set its corresponding entry to zero
                for j = 1:length(idx)
                    temp = obj.paras.MBM.ht(:,i);
                    temp(temp==idx(j)) = 0;
                    obj.paras.MBM.ht(:,i) = temp;
                end
            end
            
            %Remove tracks that contains no valid local hypotheses
            idx = sum(obj.paras.MBM.ht,1)~=0;
            obj.paras.MBM.ht = obj.paras.MBM.ht(:,idx);
            obj.paras.MBM.tt = obj.paras.MBM.tt(idx);
            if isempty(obj.paras.MBM.ht)
                %Ensure the algorithm still works when all Bernoullis are
                %pruned
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
        
        function estimates = PMBM_estimator(obj,threshold)
            %PMBM_ESTIMATOR performs object state estimation in the PMBM
            %filter
            %INPUT: threshold (if exist): object states are extracted from
            %       Bernoulli components with probability of existence no
            %       less than this threhold
            %OUTPUT:estimates: cell array of estimated object trajectories
            %%%
            %First, select the multi-Bernoulli with the highest weight.
            %Second, report the mean of the Bernoulli components whose
            %existence probability is above a threshold. The length of the
            %reported trajectory is determined by the predicted death time
            %with the highest weight
            
            estimates = {};
            n = 0;
            %Find the global hypothesis with the highest weight
            [~,I] = max(obj.paras.MBM.w);
            h_best = obj.paras.MBM.ht(I,:);
            for i = 1:length(h_best)
                %Check the validity of the local hypothesis
                if h_best(i)~=0
                    Bern = obj.paras.MBM.tt{i}(h_best(i));
                    %Report estimates from Bernoulli components
                    %with large enough probability of existence
                    if Bern.r >= threshold
                        n = n + 1;
                        [~,ideath] = max(Bern.w_death);
                        t_death = Bern.t_death(ideath);
                        tlen = t_death - Bern.t_birth + 1;
                        trajectory = obj.density.expectedValue(Bern.state);
                        estimates{n,1} = trajectory(:,1:tlen);
                    end
                end
            end
        end
        
        function obj = PMBM_predict(obj,P_S,motionmodel,birthmodel,r_min)
            %PMBM_PREDICT performs PMBM prediction step
            
            %PPP predict
            obj = PPP_predict(obj,motionmodel,birthmodel,P_S);
            %MBM predict
            for i = 1:length(obj.paras.MBM.tt)
                obj.paras.MBM.tt{i} = arrayfun(@(x) Bern_predict(obj,x,motionmodel,P_S,r_min), obj.paras.MBM.tt{i});
            end
        end
        
        function obj = PMBM_update(obj,k,z,measmodel,sensormodel,gating,w_min,M)
            %PMBM_UPDATE performs PMBM update step.
            %INPUT: z: measurements --- array of size (measurement
            %       dimension x number of measurements)
            %       gating: a struct with two fields that specifies gating
            %       parameters: P_G: gating size in decimal --- scalar;
            %                   size: gating size --- scalar.
            %       wmin: hypothesis weight pruning threshold --- scalar in
            %       logarithmic scale
            %       M: maximum global hypotheses kept
            
            m = size(z,2);                      %number of measurements received
            used_meas_u = false(m,1);           %measurement indices inside the gate of undetected objects
            nu = length(obj.paras.PPP.states);  %number of mixture components in PPP intensity
            gating_matrix_u = false(m,nu);
            for i = 1:nu
                %Perform gating for each mixture component in the PPP intensity
                [~,gating_matrix_u(:,i)] = ...
                    obj.density.ellipsoidalGating(obj.paras.PPP.states(i),z,measmodel,gating.size);
                used_meas_u = used_meas_u | gating_matrix_u(:,i);
            end
            
            n_tt = length(obj.paras.MBM.tt);    %number of pre-existing hypothesis trees
            likTable = cell(n_tt,1);            %initialise likelihood table, one for each hypothesis tree
            gating_matrix_d = cell(n_tt,1);
            used_meas_d = false(m,1);           %measurement indices inside the gate of detected objects
            for i = 1:n_tt
                %number of local hypotheses in hypothesis tree i
                num_hypo = length(obj.paras.MBM.tt{i});
                %construct gating matrix
                gating_matrix_d{i} = false(m,num_hypo);
                for j = 1:num_hypo
                    %Perform gating for each local hypothesis
                    [~,gating_matrix_d{i}(:,j)] = obj.density.ellipsoidalGating(obj.paras.MBM.tt{i}(j).state,z,measmodel,gating.size);
                    used_meas_d = used_meas_d | gating_matrix_d{i}(:,j);
                end
            end
            
            %measurement indices inside the gate
            used_meas = used_meas_d | used_meas_u;
            %find indices of measurements inside the gate of undetected
            %objects but not detected objects
            used_meas_u_not_d = used_meas > used_meas_d;
            
            %Update detected objects
            %obtain measurements that are inside the gate of detected objects
            z_d = z(:,used_meas_d);
            m = size(z_d,2);
            gating_matrix_d = cellfun(@(x) x(used_meas_d,:), gating_matrix_d, 'UniformOutput',false);
            n_tt_upd = n_tt + m;                %number of hypothesis trees
            hypoTable = cell(n_tt_upd,1);       %initialise hypothesis table, one for each hypothesis tree
            for i = 1:n_tt
                %number of local hypotheses in hypothesis tree i
                num_hypo = length(obj.paras.MBM.tt{i});
                %initialise likelihood table for hypothesis tree i
                likTable{i} = -inf(num_hypo,m+1);
                %initialise hypothesis table for hypothesis tree i
                hypoTable{i} = cell(num_hypo*(m+1),1);
                for j = 1:num_hypo
                    %Missed detection
                    [hypoTable{i}{(j-1)*(m+1)+1},likTable{i}(j,1)] = Bern_undetected_update(obj,[i,j],sensormodel.P_D);
                    %Update with measurement
                    likTable{i}(j,[false;gating_matrix_d{i}(:,j)]) = ...
                        Bern_detected_update_lik(obj,[i,j],z_d(:,gating_matrix_d{i}(:,j)),measmodel,sensormodel.P_D);
                    for jj = 1:m
                        if gating_matrix_d{i}(jj,j)
                            hypoTable{i}{(j-1)*(m+1)+jj+1} = Bern_detected_update_state(obj,[i,j],z_d(:,jj),measmodel);
                        end
                    end
                end
            end
            
            %Update undetected objects
            lik_new = -inf(m,1);
            gating_matrix_ud = gating_matrix_u(used_meas_d,:);
            %Create new hypothesis trees, one for each measurement inside
            %the gate 
            for i = 1:m
                if any(gating_matrix_ud(i,:))
                    [hypoTable{n_tt+i,1}{1}, lik_new(i)] = ...
                        PPP_detected_update(obj,k,gating_matrix_ud(i,:),z_d(:,i),measmodel,sensormodel.P_D,sensormodel.intensity_c);
                else
                    %For measurements not inside the gate of undetected
                    %objects, set likelihood to clutter intensity
                    lik_new(i) = log(sensormodel.intensity_c);
                end
            end
            used_meas_ud = sum(gating_matrix_ud, 2) >= 1; 
            
            %Cost matrix for first detection of undetected objects
            L2 = inf(m);
            L2(logical(eye(m))) = -lik_new;
            
            %Update global hypothesis
            w_upd = [];             
            ht_upd = zeros(0,n_tt_upd);
            H_upd = 0;
            
            %Number of global hypothesis
            H = length(obj.paras.MBM.w);
            if H == 0 %if there is no pre-existing hypothesis tree
                w_upd = 0;
                H_upd = 1;
                ht_upd = zeros(1,m);
                ht_upd(used_meas_ud) = 1;
            else
                for h = 1:H
                    %Cost matrix for detected objects
                    L1 = inf(m,n_tt);
                    lik_temp = 0;
                    for i = 1:n_tt
                        hypo_idx = obj.paras.MBM.ht(h,i);
                        if hypo_idx~=0
                            L1(:,i) = -(likTable{i}(hypo_idx,2:end) - likTable{i}(hypo_idx,1));
                            %we need add the removed weights back to
                            %calculate the updated global hypothesis weight
                            lik_temp = lik_temp + likTable{i}(hypo_idx,1);
                        end
                    end
                    %Cost matrix of size m-by-(n+m)
                    L = [L1 L2];
                    
                    if isempty(L)
                        %Consider the case that no measurements are inside
                        %the gate, thus missed detection
                        gainBest = 0;
                        col4rowBest = 0;
                    else
                        %Obtain M best assignments using Murty's algorithm
                        [col4rowBest,~,gainBest] = kBest2DAssign(L,ceil(exp(obj.paras.MBM.w(h)+log(M))));
                        %Obtain M best assignments using Gibbs sampling
%                       [col4rowBest,gainBest] = assign2DByGibbs(L,100,ceil(exp(obj.paras.MBM.w(h)+log(M))));
                    end
                    
                    %Restore weights
                    w_upd = [w_upd;-gainBest+lik_temp+obj.paras.MBM.w(h)];
                    
                    %Update global hypothesis look-up table
                    Mh = length(gainBest);
                    ht_upd_h = zeros(Mh,n_tt_upd);
                    for j = 1:Mh
                        ht_upd_h(j,1:n_tt_upd) = 0;
                        for i = 1:n_tt
                            if obj.paras.MBM.ht(h,i) ~= 0
                                idx = find(col4rowBest(:,j)==i, 1);
                                if isempty(idx)
                                    %missed detection
                                    ht_upd_h(j,i) = (obj.paras.MBM.ht(h,i)-1)*(m+1)+1;
                                else
                                    %measurement update
                                    ht_upd_h(j,i) = (obj.paras.MBM.ht(h,i)-1)*(m+1)+idx+1;
                                end
                            end
                        end
                        for i = n_tt+1:n_tt_upd
                            idx = find(col4rowBest(:,j)==i, 1);
                            if ~isempty(idx) && used_meas_ud(idx)
                                %measurement update for PPP
                                ht_upd_h(j,i) = 1;
                            end
                        end
                    end
                    H_upd = H_upd + Mh;
                    ht_upd = [ht_upd;ht_upd_h];
                end
                
                %Normalize global hypothesis weights
                 w_upd = normalizeLogWeights(w_upd);
                
            end
            
            %Append new hypothesis trees that created by measurements
            %inside the gate of undetected objects but not detected objects
            z_u_not_d = z(:,used_meas_u_not_d);
            num_u_not_d = size(z_u_not_d,2);
            gating_matrix_u_not_d = gating_matrix_u(used_meas_u_not_d,:);
            for i = 1:num_u_not_d
                [hypoTable{n_tt_upd+i,1}{1}, ~] = ...
                    PPP_detected_update(obj,k,gating_matrix_u_not_d(i,:),z_u_not_d(:,i),measmodel,sensormodel.P_D,sensormodel.intensity_c);
            end
            ht_upd = [ht_upd ones(H_upd,num_u_not_d)];
            
            %Update undetected objects with missed detection
            obj = PPP_undetected_update(obj,sensormodel.P_D);
            
            %Prune hypotheses with weight smaller than the specified
            %threshold 
            [w_upd, hypo_idx] = hypothesisReduction.prune(w_upd,1:H_upd,w_min);
            ht_upd = ht_upd(hypo_idx,:);
            w_upd = normalizeLogWeights(w_upd);
            
            %Keep at most M hypotheses with the highest weights
            [w_upd, hypo_idx] = hypothesisReduction.cap(w_upd,1:length(w_upd),M);
            ht_upd = ht_upd(hypo_idx,:);
            obj.paras.MBM.w = normalizeLogWeights(w_upd);
            
            %Remove empty hypothesis trees
            if ~isempty(ht_upd)
                idx = sum(ht_upd,1) >= 1;
                ht_upd = ht_upd(:,idx);
                hypoTable = hypoTable(idx);
                n_tt_upd = size(ht_upd,2);
            end
            
            %Prune local hypotheses that do not appear in maintained global
            %hypotheses 
            obj.paras.MBM.tt = cell(n_tt_upd,1);
            for i = 1:n_tt_upd
                temp = ht_upd(:,i);
                hypoTableTemp = hypoTable{i}(unique(temp(temp~=0), 'stable'));
                obj.paras.MBM.tt{i} = [hypoTableTemp{:}]';
            end
            
            %Re-index hypothesis table
            for i = 1:n_tt_upd
                idx = ht_upd(:,i) > 0;
                [~,~,ht_upd(idx,i)] = unique(ht_upd(idx,i),'rows','stable');
            end
            
            obj.paras.MBM.ht = ht_upd;
            
        end
        
    end
end
