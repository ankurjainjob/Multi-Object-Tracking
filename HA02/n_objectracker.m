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

        % STEPS FOR GNN
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
                for j = find(idx_z_ingate(i,:))
                    S_i_h    = measmodel.H(states(i).x) * states(i).P * measmodel.H(states(i).x).';
                    zbar_i_h = measmodel.h(states(i).x);
                    L(i,j) = -( log(sensormodel.P_D/sensormodel.intensity_c) ...
                                 -1/2*log(det(2*pi*S_i_h)) ...
                                 -1/2*(z(:,j) - zbar_i_h).' / S_i_h * (z(:,j) - zbar_i_h)  );
                end
                L(i,m+i) = - log(1-sensormodel.P_D);
            end

            % 3. find the best assignment matrix using a 2D assignment solver;
            % Murty's algorithm
            [col4row,~,gain] = assign2D(L);
            assert(gain~=-1, 'Assignment problem is unfeasible');

            % 4. create new local hypotheses according to the best assignment matrix obtained;
            for i=1:n
                % if object i was assigned to a meas. => KALMAN UPDATE.
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


    function [estimates_x, estimates_P] = JPDAfilter(obj, states, Z, sensormodel, motionmodel, measmodel)
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

        % STEPS FOR JPDA
        % 1. implement ellipsoidal gating for each local hypothesis seperately;
        % 2. construct 2D cost matrix of size (number of objects, number of measurements that at least fall inside the gates + number of objects);
        % 3. find the M best assignment matrices using a M-best 2D assignment solver;
        % 4. normalise the weights of different data association hypotheses;
        % 5. prune assignment matrices that correspond to data association hypotheses with low weights and renormalise the weights;
        % 6. create new local hypotheses for each of the data association results;
        % 7. merge local hypotheses that correspond to the same object by moment matching;
        % 8. extract object state estimates;
        % 9. predict each local hypothesis.

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
            mk = size(z,2);

            % 1. implement ellipsoidal gating for each predicted local hypothesis seperately, see Note below for details; 
            idx_z_ingate = zeros(n,mk);
            for i=1:n
                [~, idx_z_ingate(i,:)] = obj.density.ellipsoidalGating(states(i), z, measmodel, obj.gating.size);
            end
            % 1.1 disconsider measurements that do not fall inside any object gates
            idx_keep = sum(idx_z_ingate,1) > 0;
            z = z(:,idx_keep);
            idx_z_ingate = idx_z_ingate(:,idx_keep);
            mk = sum(idx_keep);

            % 2. construct 2D cost matrix of size (number of objects, number of measurements that at least fall inside the gates + number of objects);
            L = inf(n,mk+n);
            for i=1:n
                for j = find(idx_z_ingate(i,:))
                    S_i_h    = measmodel.H(states(i).x) * states(i).P * measmodel.H(states(i).x).';
                    zbar_i_h = measmodel.h(states(i).x);
                    L(i,j) = -( log(sensormodel.P_D/sensormodel.intensity_c) ...
                                 -1/2*log(det(2*pi*S_i_h)) ...
                                 -1/2*(z(:,j) - zbar_i_h).' / S_i_h * (z(:,j) - zbar_i_h)  );
                end
                L(i,mk+i) = - log(1-sensormodel.P_D);
            end

            % 3. find the best assignment matrix using a 2D assignment solver (Murty's algorithm);
            M = obj.reduction.M;    % number of hypothesis at each step
            [Theta,~,gain] = kBest2DAssign(L,M);
            M = length(gain);       % there might be not enough hypothesis available
            assert( all(gain~=-1), 'Assignment problem is unfeasible');

            % 3.1 calculate weight for each hypothesis
            log_w = zeros(M,1);
            for iM =1:M
                tr_AL = sum(L(sub2ind(size(L),1:n,Theta(:,iM)')));     % same as trace(A'*L)
                % exp(-trace(A'*L)) gives in abs, but we want to keep w in log scale
                % this is equal to multiply each selected weights
                log_w(iM) = -tr_AL;    
            end
            Theta( Theta > mk ) = mk+1; % set misdetection hypothesis to index = mk+1

            % 4. normalise the weights of different data association hypotheses;
            log_w = normalizeLogWeights(log_w);

            % 5. prune assignment matrices that correspond to data association hypotheses with low weights and renormalise the weights;
            hyp = 1:M;  % indices for hypothesis (each column in col4row is one hypothesis)
            [log_w, hyp] = hypothesisReduction.prune( log_w, hyp, obj.reduction.w_min );
            % remove pruned hypothesis
            Theta = Theta(:,hyp);       % each column is a hypothesis
            log_w = normalizeLogWeights(log_w);


            % 6. create new local hypotheses for each of the data association results;
            beta = zeros(n,mk+1);   % marg. prob that a object i=1:n is associated to meas. j=0:m
            for i=1:n
                for i_Theta = 1:size(Theta,2)
                    j = Theta(i,i_Theta);      % j=1 means ass. to meas. 1, j=mk+1 means misdetection
                    beta(i,j) = beta(i,j) + exp( log_w(i_Theta)  );
                end
            end
            % sanity check: sum of beta over j = 1 (each row should sum 1)

            % 7. merge local hypotheses that correspond to the same object by moment matching;  
            for i=1:n
                P_pred = states(i).P;
                x_pred = states(i).x;
                H = measmodel.H(x_pred);
                S = H * P_pred * H' + measmodel.R;
                K = P_pred * H' / S;
                ksi_ij = cell(mk,1);            % innovation mean
                ksi_i  = zeros(size(z,1),1);    % expected innovation for each object
                aux    = zeros(size(z,1),size(z,1));
                for j=1:mk
                    ksi_ij{j} = z(:,j) - measmodel.h( x_pred );
                    ksi_i = ksi_i + beta(i,j) * ksi_ij{j};
                    aux = aux + beta(i,j) * ksi_ij{j} * ksi_ij{j}';
                end
                % update mean
                states(i).x = x_pred + K * ksi_i;
                states(i).P = beta(i,mk+1) * P_pred + ...
                              (1 - beta(i,mk+1)) * P_pred - K * S * K' + ...
                              K * ( aux - ksi_i*ksi_i' ) * K';
            end

            % 8. extract object state estimates;
            for i=1:n
                estimates_x{k}(:,i) = states(i).x;
                estimates_P{k}(:,:,i) = states(i).P;
            end

            % 9. predict each local hypothesis.
            states = arrayfun(@(s) obj.density.predict(s,motionmodel), states );
        end
    end 


    function [estimates_x, estimates_P] = TOMHT(obj, states, Z, sensormodel, motionmodel, measmodel)
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

        
        % STEPS FOR TOMHT
        % 1. for each local hypothesis in each hypothesis tree:
            % 1.1. implement ellipsoidal gating;
        % 2. disconsider measurements that do not fall inside any local hypothesis gate
        % 3. for each local hypothesis in each hypothesis tree:
            % 3.1. calculate missed detection and predicted likelihood for each measurement inside the gate and make sure to save these for future use; 
            % 3.2. create updated local hypotheses and make sure to save how these connects to the old hypotheses and to the new the measurements for future use;
        % 4. for each predicted global hypothesis: 
            % 4.1. create 2D cost matrix; 
            % 4.2. obtain M best assignments using a provided M-best 2D assignment solver; 
            % 4.3. update global hypothesis look-up table according to the M best assignment matrices obtained and use your new local hypotheses indexing;
        % 5. normalise global hypothesis weights and implement hypothesis reduction technique: pruning and capping;
        % 6. prune local hypotheses that are not included in any of the global hypotheses;
        % 7. Re-index global hypothesis look-up table;
        % 8. extract object state estimates from the global hypothesis with the highest weight;
        % 9. predict each local hypothesis in each hypothesis tree.


        % ------- ALLOCATE MEMORY -------
        
        % number of time steps
        N = numel(Z);
        % number of objects
        n = numel(states);   

        % allocate memory
        estimates_x = cell(N,1);
        estimates_P = cell(N,1);

        % global hypothesis lookup-table
        old_H = ones(1,n);
        % local hypothesis pdfs for each object - tree leaves in the local hypothesis trees
        old_H_i = cell(1,n);
        for i=1:n
            old_H_i{i} = states(i);
        end
        % hypothesis weights - at the start we only have one hypothesis = prior
        old_log_w = log(1);

        
        wbar = waitbar(0,sprintf('Calculating HO-MHT iterations - k=%d',0));
        
        
        % ------- START TRACKING -------
        
        for k=1:N

            % measurements at time t=k
            z = Z{k};
            % number of measurements
            m = size(z,2);

            DOIT = 1;
            
            % 1. for each local hypothesis in each hypothesis tree:
            idx_z_ingate = cell(1,n);  % structure: {object}[meas, local hyp. for object i]
            for i=1:n
                n_i = length(old_H_i{i});
                idx_z_ingate{i} = zeros(m,n_i);
                for lh=1:n_i
                    % 1.1. implement ellipsoidal gating;
                    [~, idx_z_ingate{i}(:,lh)] = obj.density.ellipsoidalGating(old_H_i{i}(lh), z, measmodel, obj.gating.size);
                end
            end
            % 2. disconsider measurements that do not fall inside any local hypothesis gate
            idx_clutter = sum(cell2mat(idx_z_ingate)') == 0;
            z = z(:,~idx_clutter);
            m = sum(~idx_clutter);
            idx_z_ingate = cellfun( @(idx_i) idx_i(~idx_clutter,:) ,idx_z_ingate, 'UniformOutput' ,false);
            
            
            % 3. for each local hypothesis in each hypothesis tree:
            H_i = cell(1,n);
            log_w_i = cell(1,n);
            for i=1:n
                n_i = length(old_H_i{i});  % number of local hypothesis for obj i                 
                log_w_i{i} = -inf(1,n_i*(m+1));     % init vector with (log(0) = -Inf)
                for lh=1:n_i       % lh = local hypothesis
                    for j=find(idx_z_ingate{i}(:,lh))'
                        newidx = (lh-1)*(m+1) + j;  % index of the new local hypothesis
                        
                        % 3.1. calculate predicted likelihood for each measurement inside the gate
                        S_i_h    = measmodel.H(old_H_i{i}(lh).x) * old_H_i{i}(lh).P * measmodel.H(old_H_i{i}(lh).x).';
                        zbar_i_h = measmodel.h(old_H_i{i}(lh).x);
                        log_w_i{i}(newidx) = log(sensormodel.P_D/sensormodel.intensity_c) ...
                                            -1/2*log(det(2*pi*S_i_h)) ...
                                            -1/2*(z(:,j) - zbar_i_h).' / S_i_h * (z(:,j) - zbar_i_h);
                        % 3.2. create updated local hypotheses
                        H_i{i}(newidx) = obj.density.update(old_H_i{i}(lh), z(:,j) , measmodel);
                    end
                    newidx = (lh)*(m+1);
                    % 3.1. calculate missed detection likelihood
                    log_w_i{i}(newidx) = log(1-sensormodel.P_D);
                    % 3.2. create updated misdetectionlocal hypotheses
                    H_i{i}(newidx) = old_H_i{i}(lh);
                end
            end
            
            % 4. for each old predicted global hypothesis:
            log_w = [];
            H     = [];
            for h=1:size(old_H,1)
                    
                % 4.1. create 2D cost matrix;
                L = inf(n,m+n);
                for i=1:n
                    % given a global hypothesis lets find the correspondent
                    % local hypothesis weights in the vector for each object
                    startingidx = (old_H(h,i)-1)*(m+1) +1;
                    finalidx    = old_H(h,i)*(m+1) -1;
                    L(i,1:m) = -log_w_i{i}( startingidx:finalidx );     % must have size == m
                    L(i,m+i) = -log_w_i{i}( finalidx+1 );               % size == 1 (misdetection hyp.)
                end

                % 4.2. obtain M best assignments using a provided M-best 2D assignment solver; 
                M = obj.reduction.M;    % number of hypothesis at each step
                [Theta,~,gain] = kBest2DAssign(L,M);
                M = length(gain);       % there might be not enough hypothesis available
                assert( all(gain~=-1), 'Assignment problem is unfeasible');

                % 4.3. update global hypothesis look-up table according to the M best assignment matrices obtained and use your new local hypotheses indexing;
                for iM =1:M
                    tr_AL = sum(L(sub2ind(size(L),1:n,Theta(:,iM)')));     % same as trace(A'*L)
                    % exp(-trace(A'*L)) gives in abs, but we want to keep w in log scale
                    % this is equal to multiply each selected weights
                    log_w(end+1,1) = old_log_w(h) + (-tr_AL);   % wnew = wold * w_i=1 * w_i=2 * ...
                    
                    Theta(Theta(:,iM)>m, iM) = m+1;     % set all misdetection hypothesis to the same index
                    
                    % update global hypothesis look-up table
                    H(end+1,1:n) = zeros(1,n);
                    for i=1:n
                        H(end,i) = (old_H(h,i)-1)*(m+1) + Theta(i,iM);
                    end
                end
            end

            % 5. normalise global hypothesis weights and implement hypothesis reduction technique: pruning and capping;
            log_w = normalizeLogWeights(log_w);
            % prune
            [log_w, hyp] = hypothesisReduction.prune( log_w, 1:length(log_w), obj.reduction.w_min );
            H = H(hyp,:);
            log_w = normalizeLogWeights(log_w);
            % cap
            [log_w, hyp] = hypothesisReduction.cap( log_w, 1:length(log_w), obj.reduction.M );
            H = H(hyp,:);
            log_w = normalizeLogWeights(log_w);
            
            
            for i=1:n
                % 6. prune local hypotheses that are not included in any of the global hypotheses;
                hyp_keep = unique(H(:,i));
                H_i{i} = H_i{i}(hyp_keep);
                log_w_i{i} = log_w_i{i}(hyp_keep);
                log_w_i{i} = normalizeLogWeights(log_w_i{i});
                % 7. Re-index global hypothesis look-up table;
                for ri=1:numel(hyp_keep)
                    H( H(:,i) == hyp_keep(ri), i) = ri;
                end
            end

            % 8. extract object state estimates from the global hypothesis with the highest weight;
            [~,idx_best] = max(log_w);
            for i=1:n
                estimates_x{k}(:,i)   = H_i{i}( H(idx_best,i) ).x;
                estimates_P{k}(:,:,i) = H_i{i}( H(idx_best,i) ).P;
            end
            
            % 9. predict each local hypothesis in each hypothesis tree.
            for i=1:n
                H_i{i} = arrayfun(@(s) obj.density.predict(s,motionmodel), H_i{i} );
            end
            
            old_H = H;
            old_H_i = H_i;
            old_log_w = log_w;
 
            
            waitbar(k/N, wbar, sprintf('Calculating HO-MHT iterations - k=%d/%d',k,N));
        end
        close(wbar);

    end
    
end
end


% backup
% 10. plot
%     x1 = -1000:20:900;
%     x2 = -1000:20:1000;
%     [X1,X2] = meshgrid(x1,x2);
%     for iM=1:M
%         for i=1:n
%             pX = 0*X1;
%             mean = H_i{i}(H(iM,i)).x(1:2);
%             cov  = H_i{i}(H(iM,i)).P(1:2,1:2);
%             cov = (cov+cov')/2;
%             w    = exp( log_w(i) );
%             for ix1=1:numel(x1)
%                 pX(:,ix1) = w * mvnpdf([X1(:,ix1),X2(:,ix1)],mean',cov);
%             end
%             subplot(2,2,i)
%             surf(X1,X2,pX);
%             shading interp 
%             view(2);
%         end
%     end

