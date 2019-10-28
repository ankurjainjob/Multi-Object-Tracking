classdef multiobjectracker
    %MULTIOBJECTRACKER is a class containing functions to track multiple
    %object in clutter. 
    %Model structures need to be called:
    %sensormodel: a structure specifies the sensor parameters
    %           P_D: object detection probability --- scalar
    %           lambda_c: average number of clutter measurements per time
    %           scan, Poisson distributed --- scalar 
    %           pdf_c: value of clutter pdf --- scalar
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
        P_S         %object survival probability
    end
    
    methods
        
        function obj = initialize(obj,density_class_handle,P_S,P_G,m_d,w_min,merging_threshold,M)
            %INITIATOR initializes singleobjectracker class
            %INPUT: density_class_handle: density class handle
            %       P_S: object survival probability
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
            %         --- scalar 
            %         obj.P_S: object survival probability
            obj.density = density_class_handle;
            obj.gating.P_G = P_G;
            obj.gating.size = chi2inv(obj.gating.P_G,m_d);
            obj.reduction.w_min = log(w_min);
            obj.reduction.merging_threshold = merging_threshold;
            obj.reduction.M = M;
            obj.P_S = P_S;
        end
        
        function estimates = GMPHDfilter(obj, birthmodel, Z, sensormodel, motionmodel, measmodel)
            %GMPHDFILTER tracks multiple objects using Gaussian mixture
            %probability hypothesis density filter 
            %INPUT: birthmodel: object birth model: structure array of size
            %       (1, number of hypothesised new born objects per time
            %       step) with three fields:  
            %                   w: object birth weight --- scalar in
            %                   logarithm domain
            %                   x: object initial state mean --- (object
            %                   state dimension) x 1 vector 
            %                   P: object initial state covariance ---
            %                   (object state dimension) x (object state
            %                   dimension) matrix  
            %       Z: cell array of size (total tracking time, 1), each
            %       cell stores measurements of size (measurement
            %       dimension) x (number of measurements at corresponding
            %       time step)  
            %OUTPUT:estimates: cell array of size (total tracking time, 1),
            %       each cell stores estimated object state of size (object
            %       state dimension) x (number of objects)  
            
            %Preallocate memory
            K = length(Z);
            estimates = cell(K,1);
            
            %Create a class instance
            PPP = PHDfilter();
            %Initialize the PPP
            PPP = initialize(PPP,obj.density,birthmodel);
            
            for k = 1:K
                %PPP update
                PPP = update(PPP,Z{k},measmodel,sensormodel.intensity_c,sensormodel.P_D,obj.gating);
                %PPP approximation
                PPP = componentReduction(PPP,obj.reduction);
                %Extract state estimates from the PPP
                estimates{k} = PHD_estimator(PPP);
                %PPP prediction
                PPP = predict(PPP,motionmodel,obj.P_S,birthmodel);
            end

        end
        
    end
end

