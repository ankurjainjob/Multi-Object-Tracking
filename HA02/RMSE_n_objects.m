function rmse = RMSE_n_objects(X,Y)

% Compute RMSE distance between X and Y estimates
% Inputs: X,Y - cell array, each cell contains matrices of column vectors
% Output: scalar distance between X and Y
% Note: the Euclidean 2-norm is used as the "base" distance on the region

K = length(X);
dist = zeros(K,1);
for k = 1:K
    %Calculate sizes of the input point patterns
    n = size(X{k},2);
    
    %Calculate cost/weight matrix for pairings - fast method with vectorization
    XX= repmat(X{k},[1 n]);
    YY= reshape(repmat(Y{k},[n 1]),[size(Y{k},1) n*n]);
    D = reshape(sqrt(sum((XX-YY).^2)),[n n]);
    
    %Compute optimal assignment and cost
    [~,~,cost]= assign2D(D);
    
    %Calculate final distance
    dist(k)= cost/n;
end
rmse = mean(dist);

end