classdef hypothesisReduction
    %HYPOTHESISREDUCTION is class containing different hypotheses reduction
    %method 
    %PRUNE: prune hypotheses with small weights.
    %CAP:   keep M hypotheses with the highest weights and discard the rest. 
    %MERGE: merge similar hypotheses in the sense of small Mahalanobis
    %       distance.
    %
    % Note: Usually, one need to re-normalise the weights after pruning/capping; 
    %       however, this is not done inside the functions you are going to implement.
    
    methods (Static)
        function [hypothesesWeight, multiHypotheses] = ...
                prune(hypothesesWeight, multiHypotheses, threshold)
            %PRUNE prunes hypotheses with small weights
            %INPUT: hypothesesWeight: the weights of different hypotheses
            %       in logarithmic scale --- (number of hypotheses) x 1
            %       vector
            %       multiHypotheses: (number of hypotheses) x 1 structure
            %       threshold: hypotheses with weights smaller than this
            %       threshold will be discarded --- scalar in logarithmic scale
            %OUTPUT:hypothesesWeight: hypotheses weights after pruning in
            %       logarithmic scale --- (number of hypotheses after
            %       pruning) x 1 vector   
            %       multiHypotheses: (number of hypotheses after pruning) x
            %       1 structure
            idx_keep = hypothesesWeight > threshold;
            % select
            hypothesesWeight = hypothesesWeight(idx_keep);
            % select indices
            multiHypotheses = multiHypotheses(idx_keep);
        end

        function [hypothesesWeight, multiHypotheses] = ...
                cap(hypothesesWeight, multiHypotheses, M)
            %CAP keeps M hypotheses with the highest weights and discard
            %the rest 
            %INPUT: hypothesesWeight: the weights of different hypotheses
            %       in logarithmic scale --- (number of hypotheses) x 1
            %       vector  
            %       multiHypotheses: (number of hypotheses) x 1 structure
            %       M: only keep M hypotheses --- scalar
            %OUTPUT:hypothesesWeight: hypotheses weights after capping in
            %       logarithmic scale --- (number of hypotheses after
            %       capping) x 1 vector 
            %       multiHypotheses: (number of hypotheses after capping) x
            %       1 structure 
            
            % sort by descending hypothesis weight
            [w_sorted,idx] = sort(hypothesesWeight,'descend');
            % keep best M hypothesis
            M = min(M,length(hypothesesWeight));
            idx_keep = idx(1:M);
            % select indices
            hypothesesWeight = hypothesesWeight(idx_keep);
            multiHypotheses = multiHypotheses(idx_keep);
        end

        function [hypothesesWeight,multiHypotheses] = ...
                merge(hypothesesWeight,multiHypotheses,threshold,density)
            %MERGE merges hypotheses within small Mahalanobis distance
            %INPUT: hypothesesWeight: the weights of different hypotheses
            %       in logarithmic scale --- (number of hypotheses) x 1
            %       vector  
            %       multiHypotheses: (number of hypotheses) x 1 structure
            %       threshold: merging threshold --- scalar
            %       density: a class handle
            %OUTPUT:hypothesesWeight: hypotheses weights after merging in
            %       logarithmic scale --- (number of hypotheses after
            %       merging) x 1 vector  
            %       multiHypotheses: (number of hypotheses after merging) x
            %       1 structure 
            
            [hypothesesWeight,multiHypotheses] = ...
                density.mixtureReduction(hypothesesWeight,multiHypotheses,threshold);

        end  

    end
end



%% test 

% %Note that multiHypotheses does not really need to be a struct array to make the function work properly
% multiHypotheses = 1:100;
% %Generate some random weights
% hypothesesWeight = rand(100,1);
% hypothesesWeight = log(hypothesesWeight/sum(hypothesesWeight));
% %Pruning
% [hypothesesWeight_hat, multiHypotheses_hat] = hypothesisReduction.prune(hypothesesWeight, multiHypotheses, log(1e-2)); 
% %Capping
% [hypothesesWeight_hat, multiHypotheses_hat] = hypothesisReduction.cap(hypothesesWeight, multiHypotheses, 50); 
