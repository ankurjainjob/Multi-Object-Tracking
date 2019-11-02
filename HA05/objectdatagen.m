function objectdata = objectdatagen(groundtruth,motionmodel,ifnoisy)
%TARGETDATA generates groundtruth object data
%INPUT:  groundtruth specifies the parameters used to generate groundtruth
%        motionmodel: a structure specifies the motion model parameters
%        ifnoisy: boolean value indicating whether to generate noisy object
%        state sequence or not 
%OUTPUT: objectdata.X:  (K x 1) cell array, each cell stores object states
%        of size (object state dimension) x (number of objects at
%        corresponding time step)  
%        objectdata.N:  (K x 1) cell array, each cell stores the number of
%        objects at corresponding time step 

%Preallocate memory
K = groundtruth.K;
objectdata.X = cell(K,1);
objectdata.N = zeros(K,1);

for i = 1:groundtruth.nbirths
    objectstate = groundtruth.xstart(:,i);
    for k = groundtruth.tbirth(i):min(groundtruth.tdeath(i),K)
        if ifnoisy
            objectstate = mvnrnd(motionmodel.f(objectstate), motionmodel.Q)';
        else
            objectstate = motionmodel.f(objectstate);
        end
        objectdata.X{k} = [objectdata.X{k} objectstate];
        objectdata.N(k) = objectdata.N(k) + 1;
    end
end

end