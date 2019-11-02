classdef singleobjectracker
    %SINGLEOBJECTRACKER is a class containing functions to track a single
    %object in clutter. 
    %Model structures need to be called:
    %sensormodel: a structure specifies the sensor parameters
    %           P_D: object detection probability --- scalar
    %           lambda_c: average number of clutter measurements per time
    %           scan, Poisson distributed --- scalar 
    %           pdf_c: clutter (Poisson) density --- scalar
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
            %INITIATOR initializes singleobjectracker class
            %INPUT: density_class_handle: density class handle
            %       P_G: gating size in decimal --- scalar
            %            in this case, P_G = Pr[ d^2 <= G ] = significance
            %            level. The larger P_G, the G. P_G here is the
            %            significance level. Common values is 0.99
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
            %         --- scalar 
            
            obj.density = density_class_handle;
            obj.gating.P_G = P_G;
            obj.gating.size = chi2inv(obj.gating.P_G,m_d);
            obj.reduction.w_min = log(w_min);
            obj.reduction.merging_threshold = merging_threshold;
            obj.reduction.M = M;
        end
        
        function [estimates_x, estimates_P] = nearestNeighbourFilter(obj, state, Z, sensormodel, motionmodel, measmodel)
            %NEARESTNEIGHBOURFILTER tracks a single object using nearest
            %neighbor association 
            %INPUT: state: a structure with two fields:
            %                x: object initial state mean --- (object state
            %                dimension) x 1 vector 
            %                P: object initial state covariance --- (object
            %                state dimension) x (object state dimension)
            %                matrix  
            %       Z: cell array of size (total tracking time, 1), each
            %       cell stores measurements of  
            %            size (measurement dimension) x (number of
            %            measurements at corresponding time step) 
            %OUTPUT:estimates: cell array of size (total tracking time, 1),
            %       each cell stores estimated object state of size (object
            %       state dimension) x 1   
            
            N = numel(Z);
            
            % allocate memory
            estimates_x = cell(N,1);
            estimates_P = cell(N,1);
            
            % first predicted state is the prior distribution
            state_pred = state;
            
            for i=1:N
                z = Z{i};
                % Apply gating
                [z_ingate, ~] = obj.density.ellipsoidalGating(state_pred, z, measmodel, obj.gating.size);
                
                % number of hypothesis
                mk = size(z_ingate,2) + 1;
                
                % calculate predicted likelihood = N(z; zbar, S) => scalar value
                predicted_likelihood = exp(obj.density.predictedLikelihood(state_pred,z_ingate,measmodel));
                
                % EVALUATE HYPOTHESIS
                % detection
                w_theta_k = sensormodel.P_D * predicted_likelihood / sensormodel.intensity_c;
                % missdetection
                w_theta_0 = 1 - sensormodel.P_D;
                
                % UPDATE
                [max_w_theta, max_theta] = max(w_theta_k);
                if mk==1 || w_theta_0 > max_w_theta
                    state_upd = state_pred;
                else
                    z_NN = z_ingate(:,max_theta);    % nearest neighbour measurement
                    state_upd = obj.density.update(state_pred, z_NN, measmodel);
                end
                
                % ESTIMATES
                estimates_x{i} = state_upd.x;
                estimates_P{i} = state_upd.P;
                
                % PREDICTION
                state_pred = obj.density.predict(state_upd, motionmodel);
            end
        end
        
        
        function [estimates_x, estimates_P] = probDataAssocFilter(obj, state, Z, sensormodel, motionmodel, measmodel)
            %PROBDATAASSOCFILTER tracks a single object using probalistic
            %data association 
            %INPUT: state: a structure with two fields:
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
            %       state dimension) x 1  
            
            N = numel(Z);
            
            % allocate memory
            estimates_x = cell(N,1);
            estimates_P = cell(N,1);
            
            % first predicted state is the prior distribution
            p_old = state;
            w_old = 0;
            
            for i=1:N
                % read measurements
                z = Z{i};

                % number of old hypothesis
                mk_old = size(w_old,1);
                
                % Erase new weight and state vectors
                w_new = [];
                p_new = struct('x',{},'P',{});
                
                % Update
                for i_old = 1:mk_old
                    % for each hypothesis, perform ellipsoidal gating and only 
                    % create object detection hypotheses for detections inside the gate;
                    [z_ingate, ~] = obj.density.ellipsoidalGating(p_old(i_old), z, measmodel, obj.gating.size);
                    mk_new = size(z_ingate,2) + 1;
                    for i_new = 1:mk_new-1
                        pred_likelihood_log = obj.density.predictedLikelihood(p_old(i_old),z_ingate(:,i_new),measmodel);
                        w_new(end+1,1) = w_old(i_old) + pred_likelihood_log + log(sensormodel.P_D/sensormodel.intensity_c);
                        p_new(end+1,1) = obj.density.update(p_old(i_old), z_ingate(:,i_new), measmodel);
                    end
                    % for each hypothesis, create missed detection hypothesis;
                    w_new(end+1,1) = w_old(i_old) + log(1 - sensormodel.P_D);
                    p_new(end+1,1) = p_old(i_old);
                end
                
                % normalise hypothesis weights;
                w_new = normalizeLogWeights(w_new);

                % prune hypotheses with small weights, and then re-normalise the weights.
                [w_new,p_new] = hypothesisReduction.prune( w_new, p_new, obj.reduction.w_min );
                w_new = normalizeLogWeights(w_new);

                % merge different hypotheses using Gaussian moment matching;
                [w_new,p_new] = hypothesisReduction.merge( w_new, p_new, obj.reduction.merging_threshold, obj.density );

                % cap the number of the hypotheses, and then re-normalise the weights;
                [w_new,p_new] = hypothesisReduction.cap( w_new, p_new, 1 );

                % extract object state estimate using the most probably hypothesis estimation;
                [~,best_w_idx] = max(w_new);
                estimates_x{i} = p_new(best_w_idx).x;
                estimates_P{i} = p_new(best_w_idx).P;

                % update
                w_old = w_new;
                p_old = p_new;
                
                % for each hypothesis, perform prediction.
                for i_old = 1:length(w_old)
                    p_old(i_old) = obj.density.predict(p_old(i_old), motionmodel);
                end   
            end
        end
        
        function [estimates_x, estimates_P] = GaussianSumFilter(obj, state, Z, sensormodel, motionmodel, measmodel)
            %GAUSSIANSUMFILTER tracks a single object using Gaussian sum
            %filtering
            %INPUT: state: a structure with two fields:
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
            %       state dimension) x 1
            
            N = numel(Z);
            
            % allocate memory
            estimates_x = cell(N,1);
            estimates_P = cell(N,1);
            
            % first predicted state is the prior distribution
            p_old = state;
            w_old = 0;
            
            for i=1:N
                % read measurements
                z = Z{i};

                % number of old hypothesis
                mk_old = size(w_old,1);
                
                % Erase new weight and state vectors
                w_new = [];
                p_new = struct('x',{},'P',{});
                
                % Update
                for i_old = 1:mk_old
                    % for each hypothesis, perform ellipsoidal gating and only 
                    % create object detection hypotheses for detections inside the gate;
                    [z_ingate, ~] = obj.density.ellipsoidalGating(p_old(i_old), z, measmodel, obj.gating.size);
                    mk_new = size(z_ingate,2) + 1;
                    for i_new = 1:mk_new-1
                        pred_likelihood_log = obj.density.predictedLikelihood(p_old(i_old),z_ingate(:,i_new),measmodel);
                        % w_new = w_old*(P_D*N(z_theta;zbar,S)/lambda_c);
                        w_new(end+1,1) = w_old(i_old) + pred_likelihood_log + log(sensormodel.P_D/sensormodel.intensity_c);
                        % p_new = Kalman_Update(p_old)
                        p_new(end+1,1) = obj.density.update(p_old(i_old), z_ingate(:,i_new), measmodel);
                    end
                    % for each hypothesis, create missed detection hypothesis
                    % w_new = w_old*(1-P_D);
                    w_new(end+1,1) = w_old(i_old) + log(1 - sensormodel.P_D);
                    % p_new = p_old
                    p_new(end+1,1) = p_old(i_old);
                end
                
                % normalise hypothesis weights;
                w_new = normalizeLogWeights(w_new);

                % prune hypotheses with small weights, and then re-normalise the weights.
                [w_new,p_new] = hypothesisReduction.prune( w_new, p_new, obj.reduction.w_min );
                w_new = normalizeLogWeights(w_new);

                % merge different hypotheses using Gaussian moment matching;
                [w_new,p_new] = hypothesisReduction.merge( w_new, p_new, obj.reduction.merging_threshold, obj.density );

                % cap the number of the hypotheses, and then re-normalise the weights;
                [w_new,p_new] = hypothesisReduction.cap( w_new, p_new, obj.reduction.M );

                % extract object state estimate using the most probably hypothesis estimation;
                [~,best_w_idx] = max(w_new);
                estimates_x{i} = p_new(best_w_idx).x;
                estimates_P{i} = p_new(best_w_idx).P;

                % update
                w_old = w_new;
                p_old = p_new;
                
                % for each hypothesis, perform prediction.
                for i_old = 1:length(w_old)
                    p_old(i_old) = obj.density.predict(p_old(i_old), motionmodel);
                end   
            end
        end
        
        
    end
end

