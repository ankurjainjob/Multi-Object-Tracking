classdef Animate_2D_tracking
    %ANIMATE_2D_TRACKING Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
    end
    
    methods(Static)
        function [ xy ] = sigmaEllipse2D( mu, Sigma, level, npoints )
            %SIGMAELLIPSE2D generates x,y-points which lie on the ellipse describing
            % a sigma level in the Gaussian density defined by mean and covariance.
            %
            %Input:
            %   MU          [2 x 1] Mean of the Gaussian density
            %   SIGMA       [2 x 2] Covariance matrix of the Gaussian density
            %   LEVEL       Which sigma level curve to plot. Can take any positive value, 
            %               but common choices are 1, 2 or 3. Default = 3.
            %   NPOINTS     Number of points on the ellipse to generate. Default = 32.
            %
            %Output:
            %   XY          [2 x npoints] matrix. First row holds x-coordinates, second
            %               row holds the y-coordinates. First and last columns should 
            %               be the same point, to create a closed curve.

            %Setting default values, in case only mu and Sigma are specified.
            if nargin < 3
                level = 3;
            end
            if nargin < 4
                npoints = 32;
            end

            % Procedure:
            % - A 3 sigma level curve is given by {x} such that (x-mux)'*Q^-1*(x-mux) = 3^2
            %      or in scalar form: (x-mux) = sqrt(Q)*3
            % - replacing z= sqrtm(Q^-1)*(x-mux), such that we have now z'*z = 3^2
            %      which is now a circle with radius equal 3.
            % - Sampling in z, we have z = 3*[cos(theta); sin(theta)]', for theta=1:2*pi
            % - Back to x we get:  x = mux  + 3* sqrtm(Q)*[cos(theta); sin(theta)]'

            xy = [];
            for ang = linspace(0,2*pi,npoints)
                xy(:,end+1) = mu + level * sqrtm(Sigma) * [cos(ang) sin(ang)]';
            end
        end
        
        function [zk,S] = estimated_meas(x,P,measmodel)
            % measurement covariance
            Hx = measmodel.H(x);
            % Innovation covariance
            S = Hx * P * Hx' + measmodel.R;
            % ensure it is positive definite
            S = (S+S')/2;
            % object measurement
            zk = measmodel.h(x);
        end
    end
    
    methods
        function obj = animate_2D_tracking(obj)
        end
        
        function animate(obj, est, initial_state, measdata, measmodel, range_c)
            fig = figure;
            hold on;
            xlim(range_c(1,:));
            ylim(range_c(2,:));
            
            n = numel(initial_state);
            
            % plot meas data
            meas_x = measdata{1}(1,:);
            meas_y = measdata{1}(2,:);
            pl_meas = scatter(meas_x, meas_y, 50, 'o', 'filled');
            pl_meas.XDataSource = 'meas_x';
            pl_meas.YDataSource = 'meas_y';
            
            % initialize objects
            for i=1:n
                [zk,S] = obj.estimated_meas(initial_state(i).x,initial_state(i).P,measmodel);

                % plot estimated covariance
                ellipse{i} = obj.sigmaEllipse2D(zk, S);
                S_x{i} = ellipse{i}(1,:); 
                S_y{i} = ellipse{i}(2,:); 
                pl_cov  = plot(S_x{i},S_y{i});
                pl_cov.XDataSource = ['S_x{',num2str(i),'}']; 
                pl_cov.YDataSource = ['S_y{',num2str(i),'}'];

                % plot estimated mean
                zk_x{i} = zk(1);
                zk_y{i} = zk(2);
                pl_mean = scatter(zk_x{i}, zk_y{i},100, 'o', 'filled');
                pl_mean.XDataSource = ['zk_x{',num2str(i),'}'];
                pl_mean.YDataSource = ['zk_y{',num2str(i),'}'];
            end
            
            
            % --- ITERATE ANIMATION ---
            for k=1:numel(measdata)
                
                % meas
                meas_x = measdata{k}(1,:);
                meas_y = measdata{k}(2,:);
                
                for i=1:n
                    [zk,S] = obj.estimated_meas(est(k).x(:,i),est(k).P(:,:,i),measmodel);
                    ellipse = obj.sigmaEllipse2D(zk, S);
                    % variance
                    S_x{i} = ellipse(1,:); 
                    S_y{i} = ellipse(2,:);
                    % mean
                    zk_x{i} = zk(1);
                    zk_y{i} = zk(2);
                end
                
                refreshdata(fig,'caller');
                drawnow;
%                 pause(0.05);
            end
        end
        
    end
end





