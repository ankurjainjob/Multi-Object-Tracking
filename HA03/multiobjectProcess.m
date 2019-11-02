classdef multiobjectProcess
    %MULTIOBJECTPROCESS is a class used to draw samples from RFSs
    %PROPERTIES: 
    %spatial_distribution: a struct with two fields:
    %   name: a string specifing a spatial distribution by its name
    %   paras: a cell array storing spatial distribution parameters:
    %          uniform distribution, size (1 x i): interval 
    %          Gaussian distribution, size (2 x i): mean and covariance
    %          Gaussian mixture, size (3 x i): means, covariances and
    %          mixing proportions 
    %          i: number of RFS components, e.g., Poisson and Bernoulli, i
    %          = 1,  multi-Bernoulli and multi-Bernoulli mixture, i =
    %          number of Bernoulli components 
    %RFS_type: a string specifying the type of RFS, e.g., 'Poisson',
    %          'Bernoulli' 
    %card_pmf: a function handle that returns the cardinality probability 
    %          mass function of the cardinality of some RFS with the given
    %          support. 
    
    %NOTE: one can use prod(1-r)*poly(-r./(1-r)) to calculate the
    %      cardinality pmf of a multi-Bernoulli process with length(r)
    %      number of Bernoulli component; the ith Bernoulli has parameter
    %      r(i). 
    
    properties
        spatial_distribution
        RFS_type
        card_pmf
    end
    
    methods
        
        function obj = spatialDistribution(obj,inputArg1,inputArg2,inputArg3)
            %SPATIALDISTRIBUTION assigns the 'spatial_distribution' property
            %INPUT: inputArg1, inputArg2, inputArg3: input arguments
            switch nargin
                case 2
                    obj.spatial_distribution.name = 'uniform';
                    for i = 1:length(inputArg1)
                        obj.spatial_distribution.paras{1}{i} = inputArg1{i};
                    end
                case 3
                    obj.spatial_distribution.name = 'Gaussian';
                    for i = 1:length(inputArg1)
                        obj.spatial_distribution.paras{1}{i} = inputArg1{i};
                        obj.spatial_distribution.paras{2}{i} = inputArg2{i};
                    end
                case 4
                    obj.spatial_distribution.name = 'Gaussian mixture';
                    for i = 1:length(inputArg1)
                        obj.spatial_distribution.paras{1}{i} = inputArg1{i};
                        obj.spatial_distribution.paras{2}{i} = inputArg2{i};
                        obj.spatial_distribution.paras{3}{i} = inputArg3{i};
                    end
            end
        end
        
        function cardStemPlot(obj, support)
            %CARDSTEMPLOT creates a stem plot of the cardinality pmf
            %INPUT: support: the support of the given cardinality pmf ---
            %vector 
            figure
            grid on
            stem(support,obj.card_pmf(support))
            xlabel('Cardinality')
            ylabel('Probability')
            title_str = ['Cardinality PMF of a ', obj.RFS_type, ' RFS'];
            title(title_str)
            set(gca,'FontSize',12) 
        end
        
        function instance = drawSamples(obj,idxParas,v)
            %DRAWSAMPLES draw samples from the given spatial distribution
            %INPUT: idxParas: index of the parameters used to specify the
            %                 spatial distribution --- scalar
            %       v: number of samples to be drawn --- scalar
            %OUTPUT:instance: drawn samples of size (d, v), d: dimension
            %                 of the spatial distribution
            %%%
            %For uniform distribution, inputArg1{i} specifies the interval
            %of a uniform distribution, a d-by-2 matrix, where d is the
            %dimension. The interval of the ith dimension is from
            %interval(i,1) to interval(i,2). 
            %%%
            %For Gaussian distribution, inputArg1{i} specifies the mean of
            %multivariate Gaussian distribution components, a m length
            %numeric vector, where m is the number of variables in each
            %component. inputArg2{i} specifies the covariances of
            %multivariate Gaussian distribution components, a m-by-m matrix.
            %%%
            %For Gaussian mixture distribution, inputArg1{i} specifies the
            %means of multivariate Gaussian distribution components, a
            %k-by-m numeric matrix, where k is the number of components and
            %m is the number of variables in each component.
            %inputArg1{i}(i,:) is the mean of component i. inputArg2{i}
            %specifies the covariances of multivariate Gaussian
            %distribution components, either a m-by-m-by-k array, where
            %inputArg2{i}(:,:,i) is the covariance matrix of component i,
            %or a m-by-m array, where covaraince can be the same across
            %components. inputArg3{i} specifies the mixing proportions of
            %mixture components, specified as a numeric vector with length
            %k. 
            switch obj.spatial_distribution.name
                case 'uniform'
                    dim = size(obj.spatial_distribution.paras{1}{idxParas},1);
                    instance = zeros(dim,v);
                    for i = 1:dim
                        instance(i,:) = unifrnd(obj.spatial_distribution.paras{1}{idxParas}(i,1),...
                            obj.spatial_distribution.paras{1}{idxParas}(i,2),1,v);
                    end
                case 'Gaussian'
                    instance = mvnrnd(obj.spatial_distribution.paras{1}{idxParas},...
                        obj.spatial_distribution.paras{2}{idxParas},v)';
                case 'Gaussian mixture'
                    gm = gmdistribution(obj.spatial_distribution.paras{1}{idxParas},...
                        obj.spatial_distribution.paras{2}{idxParas},obj.spatial_distribution.paras{3}{idxParas});
                    if v == 0
                        instance = [];
                    else
                        instance = random(gm,v)';
                    end
            end
        end
        
        function instance2Dplot(obj,instance)
            %INSTANCE2DPLOT creates a scatter plot of samples instances
            %NOTE: this function only supports 2D distributions
            figure
            if ~isempty(instance)
                plot(instance(1,:),instance(2,:),'o','LineWidth',2,'MarkerSize',5)
                legend('Samples')
                grid on
            end
            xlabel('x')
            ylabel('y')
            title_str = ['Samples of a ', obj.RFS_type, ' RFS with ', ...
                obj.spatial_distribution.name, ' PDF'];
            title(title_str)
            set(gca,'FontSize',12) 
        end
        
        function [obj, instance] = PoissonRFSs(obj,lambda)
            %POISSONRFSS draw samples from and calculate the cardinality 
            %pmf of some Poisson RFS
            %INPUT: lambda: mean of cardinality pmf --- scalar
            %OUTPUT:instance: sample instances of size (d, v), d: dimension
            %                 of the spatial distribution, v: number of
            %                 samples 
            obj.RFS_type = 'Poisson';
            obj.card_pmf = @(x) poisspdf(x, lambda);
            %draw an integer v from Poisson distirbution with parameter
            %lambda 
            v = poissrnd(lambda);
            %draw v elements from the given spatial distribution
            instance = drawSamples(obj,1,v);
        end
        
        function [obj, instance] = BernoulliRFSs(obj,r)
            %BERNOULLIRFSS draw samples from and calculate the cardinality 
            %pmf of some Bernoulli RFS
            %INPUT: r: with probability r, the RFS cardinality = 1 ---
            %scalar 
            %OUTPUT:instance: sample instances of size (d, v), d: dimension
            %                 of the spatial distribution, v: number of
            %                 samples 
            obj.RFS_type = 'Bernoulli';
            obj.card_pmf = @(x) binopdf(x,1,r);
            %draw an integer v from Bernoulli distirbution with parameter r
            v = binornd(1,r);
            %draw v elements from the given spatial distribution
            instance = drawSamples(obj,1,v);
        end
        
        function [obj, instance] = multiBernoulliRFSs(obj,r)
            %MULTIBERNOULLIRFSS draw samples from and calculate the
            %cardinality pmf of some multi-Bernoulli RFS
            %INPUT: r: vector of size (number of Bernoulli components x 1)
            %OUTPUT:instance: sample instances of size (d, v), d: dimension
            %                 of the spatial distribution, v: number of
            %                 samples 
            obj.RFS_type = 'Multi-Bernoulli';
            M = length(r);
            %draw at most M samples
            instance = cell(M,1);
            for i = 1:M
                %draw an integer v from each Bernoulli distirbution with
                %parameter r 
                v = binornd(1,r(i));
                %draw v elements from the given spatial distribution
                instance{i} = drawSamples(obj,i,v);
            end
            instance = cell2mat(instance');
            
            %calculate the cardinality pmf of a multi-Bernoulli RFS
            lr1 = length(find(r==1));
            r = r(r~=1);
            pcard = [zeros(lr1,1);prod(1-r)*poly(-r./(1-r))];
            obj.card_pmf = @(x) multiobjectProcess.createPMF(pcard,x);

        end
        
        function [obj, instance] = multiBernoulliMixtureRFSs(obj,r,p)
            %MULTIBERNOULLIRFSS draw samples from and calculate the
            %cardinality pmf of some multi-Bernoulli RFS
            %INPUT: r: cell array of size (number of multi-Bernoullis x 1),
            %          each cell contains a vector of size (number of
            %          Bernoulli components x 1)
            %       p: weights of multi-Bernoullis --- (number of
            %          multi-Bernoullis x 1) vector 
            %OUTPUT:instance: sample instances of size (d, v), d: dimension
            %                 of the spatial distribution, v: number of
            %                 samples 
            %For input convenience, here we assume that Bernoulli component
            %with the same index in different multi-Bernoulli has the same
            %pdf but possibly different probability of existence
            
            obj.RFS_type = 'Multi-Bernoulli mixture';
            M = cellfun('length',r);

            %normalise weights of different multi-Bernoullis
            p_norm = p/sum(p);
            mbidx = mnrnd(1,p_norm)==1;
            %draw an integer v from each Bernoulli distirbution with
            %parameter r of the selected multi-Bernoulli
            [~,instance] = multiBernoulliRFSs(obj,r{mbidx});
            
            num_mb = length(M);
            pcard = zeros(num_mb,max(M)+1);
            %calculate the cardinality pmf for each multi-Bernoulli RFS
            for i = 1:num_mb
                lr1 = length(find(r{i}==1));
                r{i} = r{i}(r{i}~=1);
                %append zeros to make each multi-Bernoulli's cardinality
                %pmf have equal support 
                pcard(i,:) = [zeros(1,lr1) prod(1-r{i})*poly(-r{i}./(1-r{i})) zeros(1,max(M)-M(i))];
            end
            if size(p_norm,1) > 1
                p_norm = p_norm';
            end
            %calculate the cardinality pmf of the multi-Bernoulli mixture
            pcard = sum(pcard.*p_norm',1);
            obj.card_pmf = @(x) multiobjectProcess.createPMF(pcard,x);
        end
        
    end
    
    methods (Static)
        
        function pmf = createPMF(pcard, pf)
            %CREATEPMF helps to create the card_pmf function handle
            len = length(pf);
            pmf = zeros(1,len);
            for i = 1:len
                if pf(i)+1 <= length(pcard)
                    pmf(i) = pcard(pf(i)+1);
                end
            end
        end
        
    end
    
end

