function [object_tracks] = trackgen(K,measmodel,motionmodel,sensormodel,birthmodel,P_S)
%TRACKGEN draw samples from RFS of objects (birth model, dynamic equations) 
%to simulate object trajectories. Assume that a 2D measurement model is used. 
%INPUT: K: total tracking time --- scalar
%       sensormodel: a structure specifies the sensor parameters
%           P_D: object detection probability --- scalar
%           lambda_c: average number of clutter measurements per time scan, 
%                     Poisson distributed --- scalar
%           pdf_c: clutter (Poisson) intensity --- scalar
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
%           w: weights of mixture components
%           x: mean of mixture components
%           P: covariance of mixture components
%       object survival probability --- scalar
%OUTPUT:object_tracks: a structure array specifies object tracks ---
%       (number of tracks x 1) 
%           tbirth: track initiate (object appear) time --- scalar
%           tdeath: track end time (object disappear) time --- scalar
%           x: object trajectory --- (object state dimension x time steps 
%              object exists in the scen)
%Note that if the number of tracks is zero, set the output to empty

n = 0;

%Convert birthmodel.w
[log_w,log_sum_w] = normalizeLogWeights([birthmodel.w]);

%surveillance area
range_c = sensormodel.range_c;

for k = 1:K
    %randomly sample number of births
    nb = poissrnd(exp(log_sum_w));
    for i = 1:nb
        %randomly sample birth component
        idx = find(rand<cumsum(exp(log_w)),1,'first');
        n = n+1;
        %randomly sample kinematic state
        object_tracks(n,1).x = mvnrnd(birthmodel(idx).x, birthmodel(idx).P)';
        %birth time
        object_tracks(n,1).tbirth = k;
        %death time
        object_tracks(n,1).tdeath = k;
        
        time = k;
        xk = object_tracks(n,1).x;
        
        xypos = mvnrnd(measmodel.h(xk)', measmodel.R)';
        xpos = xypos(1);
        ypos = xypos(2);
        
        ifalive = 1;
        while ifalive && (xpos>=range_c(1,1))&&(xpos<=range_c(1,2))...
                &&(ypos>=range_c(2,1))&&(ypos<=range_c(2,2))&&(time+1<=K)

            %simulate the kinematic state
            xk = mvnrnd(motionmodel.f(xk), motionmodel.Q)';
            object_tracks(n,1).x = [object_tracks(n,1).x xk];
            object_tracks(n,1).tdeath = object_tracks(n,1).tdeath + 1;
            
            %simulate measurements
            time = time + 1;
            xypos = mvnrnd(measmodel.h(xk)', measmodel.R)';
            xpos = xypos(1);
            ypos = xypos(2);
            
            if rand > P_S
                ifalive = 0;
            end
        end
    end
end

if n==0
    object_tracks = [];
end 

end