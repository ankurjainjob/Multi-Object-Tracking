function measdata = measdatagen(objectdata, sensormodel, measmodel)
%MEASDATAGEN generates object-generated measurements and clutter
%INPUT:     objectdata: a structure contains object data
%           sensormodel: a structure specifies sensor model parameters
%           measmodel: a structure specifies the measurement model
%           parameters 
%OUTPUT:    measdata: cell array of size (total tracking time, 1), each
%           cell stores measurements of size (measurement dimension) x
%           (number of measurements at corresponding time step) 

%Initialize memory
measdata = cell(length(objectdata.X),1);

%Generate measurements
for k = 1:length(objectdata.X)
    if objectdata.N(k) > 0
        idx = rand(objectdata.N(k),1) <= sensormodel.P_D;
        %Only generate object-originated observations for detected objects
        if ~isempty(objectdata.X{k}(:,idx))
            objectstates = objectdata.X{k}(:,idx);
            for i = 1:size(objectstates,2)
                meas = mvnrnd(measmodel.h(objectstates(:,i))', measmodel.R)';
                measdata{k} = [measdata{k} meas];
            end
        end
    end
    %Number of clutter measurements
    N_c = poissrnd(sensormodel.lambda_c);
    %Generate clutter
    if measmodel.d == 2
        C = repmat(sensormodel.range_c(:,1),[1 N_c])+ diag(sensormodel.range_c*[-1; 1])*rand(measmodel.d,N_c);
    elseif measmodel.d == 1
        C = (sensormodel.range_c(1,2)-sensormodel.range_c(1,1))*rand(measmodel.d,N_c)-sensormodel.range_c(1,2);
    end
    %Total measurements are the union of object detections and clutter
    measdata{k}= [measdata{k} C];                                                                  
end

end