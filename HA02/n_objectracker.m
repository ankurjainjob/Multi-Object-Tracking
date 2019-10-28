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
            %INITIATOR initializes n_objectrackersei se ela acabou operando hoje ou nao, ou sei la o que class
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
        
        function [estimates_x, estimates_P] = GNNfilter(obj, states, Z, sensormodel, motionmodel, measmodel)
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
            
            % STEPS FOR MOT
            % 0. perform prediction of each prior
            % 1. implement ellipsoidal gating for each predicted local hypothesis seperately, see Note below for details;
            % 2. construct 2D cost matrix of size (number of objects, number of measurements that at least fall inside the gates + number of objects);
            % 3. find the best assignment matrix using a 2D assignment solver;
            % 4. create new local hypotheses according to the best assignment matrix obtained;
            % 5. extract object state estimates;
            % 6. predict each local hypothesis.
            
            
            % number of time steps
            N = numel(Z);
            % number of objects
            n = numel(states);
            % allocate memory
            estimates_x = cell(N,1);
            estimates_P = cell(N,1);
            
            for k=1:N
                % measurements at time t=k
                z = Z{k};
                % number of measurements
                m = size(z,2);

                % 1. implement ellipsoidal gating for each predicted local hypothesis seperately, see Note below for details; 
                idx_z_ingate = zeros(n,m);
                for i=1:n
                    [~, idx_z_ingate(i,:)] = obj.density.ellipsoidalGating(states(i), z, measmodel, obj.gating.size);
                end
                % 1.1 disconsider measurements that do not fall inside any object gates
                idx_keep = sum(idx_z_ingate,1) > 0;
                z = z(:,idx_keep);
                idx_z_ingate = idx_z_ingate(:,idx_keep);
                m = sum(idx_keep);
                
                % 2. construct 2D cost matrix of size (number of objects, number of measurements that at least fall inside the gates + number of objects);
                L = inf(n,m+n);
                for i=1:n
                    for im = find(idx_z_ingate(i,:))
                        S_i_h    = measmodel.H(states(i).x) * states(i).P * measmodel.H(states(i).x).';
                        zbar_i_h = measmodel.h(states(i).x);
                        L(i,im) = -( log(sensormodel.P_D/sensormodel.intensity_c) ...
                                     -1/2*log(det(2*pi*S_i_h)) ...
                                     -1/2*(z(:,im) - zbar_i_h).' * inv(S_i_h) * (z(:,im) - zbar_i_h)  );
                    end
                    L(i,m+i) = - log(1-sensormodel.P_D);
                end
                
                % 3. find the best assignment matrix using a 2D assignment solver;
                % Murty's algorithm
                [col4row,~,gain] = assign2D(L);
                assert(gain~=-1, 'Assignment problem is unfeasible');
                
                % 4. create new local hypotheses according to the best assignment matrix obtained;
                for i=1:n
                    % object i was assigned to a meas. => KALMAN UPDATE.
                    % otherwise updated density = predicted density
                    if col4row(i) <= m   
                        states(i) = obj.density.update(states(i), z(:,col4row(i)), measmodel);
                    end
                end
                
                % 5. extract object state estimates;
                for i=1:n
                    estimates_x{k}(:,i) = states(i).x;
                    estimates_P{k}(:,:,i) = states(i).P;
                end
                
                % 6. predict each local hypothesis.
                states = arrayfun(@(s) obj.density.predict(s,motionmodel), states );
            end
            
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

