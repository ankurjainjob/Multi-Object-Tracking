classdef trajectoryAccumulatedStateDensity
    
    methods (Static)
        
        function expected_value = expectedValue(density)
            expected_value = density.x;
        end
        
        function covariance = covariance(density)
            covariance = density.P;
        end
        
        function density_pred = predict_state(density, motionmodel)
            %PREDICT_STATE performs linear/nonlinear (Extended) Kalman
            %prediction step (without accumulated state density)
            %INPUT: state: a structure with two fields:
            %                   x: object state mean --- (state dimension) x trajectory length matrix
            %                   P: object state covariance --- (state dimension x trajectory length) x (state dimension x trajectory length) matrix
            %       motionmodel: a structure specifies the motion model parameters
            %OUTPUT:state_pred: a structure with two fields:
            %                   x: predicted object state mean --- (state dimension) x trajectory length matrix
            %                   P: predicted object state covariance --- (state dimension x trajectory length) x (state dimension x trajectory length) matrix
            
            density_pred.x = motionmodel.f(density.x(:,end));
            density_pred.P = motionmodel.F(density.x(:,end))*density.P(end-motionmodel.d+1:end,end-motionmodel.d+1:end)*motionmodel.F(density.x(:,end))'+motionmodel.Q;
            
        end
        
        function density_pred = predict(density, motionmodel)
            %PREDICT performs linear/nonlinear (Extended) Kalman prediction step
            %INPUT: state: a structure with two fields:
            %                   x: object state mean --- (state dimension) x trajectory length matrix
            %                   P: object state covariance --- (state dimension x trajectory length) x (state dimension x trajectory length) matrix
            %       motionmodel: a structure specifies the motion model parameters
            %OUTPUT:state_pred: a structure with two fields:
            %                   x: predicted object state mean --- (state dimension) x trajectory length matrix
            %                   P: predicted object state covariance --- (state dimension x trajectory length) x (state dimension x trajectory length) matrix
            
            density_pred.x = [density.x motionmodel.f(density.x(:,end))];
            density_pred.P = [density.P density.P(:,end-motionmodel.d+1:end)*motionmodel.F(density.x(:,end))';...
                motionmodel.F(density.x(:,end))*density.P(end-motionmodel.d+1:end,:) ...
                motionmodel.F(density.x(:,end))*density.P(end-motionmodel.d+1:end,end-motionmodel.d+1:end)*motionmodel.F(density.x(:,end))'+motionmodel.Q];
            
        end
        
        function density_upd = update(density_pred, z, measmodel)
            %UPDATE performs linear/nonlinear (Extended) Kalman update step
            %INPUT: z: measurements --- (measurement dimension) x 1 vector
            %       state_pred: a structure with two fields:
            %                   x: predicted object state mean --- (state dimension) x trajectory length matrix
            %                   P: predicted object state covariance --- (state dimension x trajectory length) x (state dimension x trajectory length) matrix
            %       measmodel: a structure specifies the measurement model parameters
            %OUTPUT:state_upd: a structure with two fields:
            %                   x: updated object state mean --- (state dimension) x trajectory length matrix
            %                   P: updated object state covariance --- (state dimension x trajectory length) x (state dimension x trajectory length) matrix
            
            nx = size(density_pred.x,1);
            
            %Measurement model Jacobian
            Hx = measmodel.H(density_pred.x(:,end));
            %Innovation covariance
            S = Hx*density_pred.P(end-nx+1:end,end-nx+1:end)*Hx' + measmodel.R;
            %Make sure matrix S is positive definite
            S = (S+S')/2;
            
            K = (density_pred.P(:,end-nx+1:end)*Hx')/S;
            
            %State update
            density_upd.x = density_pred.x + reshape(K*(z - measmodel.h(density_pred.x(:,end))),[nx,size(density_pred.x,2)]);
            %Covariance update
            density_upd.P = density_pred.P - K*Hx*density_pred.P(end-nx+1:end,:);
            
        end
        
        function predict_likelihood = predictedLikelihood(density_pred,z,measmodel)
            %PREDICTLIKELIHOOD calculates the predicted likelihood in logarithm domain
            %INPUT:  z: measurements --- (measurement dimension) x (number of measurements) matrix
            %        state_pred: a structure with two fields:
            %                   x: predicted object state mean --- (state dimension) x trajectory length matrix
            %                   P: predicted object state covariance --- (state dimension x trajectory length) x (state dimension x trajectory length) matrix
            %        measmodel: a structure specifies the measurement model parameters
            %OUTPUT: meas_likelihood: measurement update likelihood for each measurement in logarithm domain --- (number of measurements) x 1 vector
            
            nx = size(density_pred.x,1);
            
            %Measurement model Jacobian
            Hx = measmodel.H(density_pred.x(:,end));
            %Innovation covariance
            S = Hx*density_pred.P(end-nx+1:end,end-nx+1:end)*Hx' + measmodel.R;
            %Make sure matrix S is positive definite
            S = (S+S')/2;
            %Calculate predicted likelihood
            predict_likelihood = log_mvnpdf(z',(measmodel.h(density_pred.x(:,end)))',S);
        end
        
        function [z_ingate, meas_in_gate] = ellipsoidalGating(density_pred, z, measmodel, gating_size)
            %ELLIPSOIDALGATING performs ellipsoidal gating for a single target
            %INPUT:  z: measurements --- (measurement dimension) x (number of measurements) matrix
            %        state_pred: a structure with two fields:
            %                   x: predicted object state mean --- (state dimension) x trajectory length matrix
            %                   P: predicted object state covariance --- (state dimension x trajectory length) x (state dimension x trajectory length) matrix
            %        measmodel: a structure specifies the measurement model parameters
            %        gating_size: gating size --- scalar
            %OUTPUT: z_ingate: measurements in the gate --- (measurement dimension) x (number of measurements in the gate) matrix
            %        meas_in_gate: boolean vector indicating whether the
            %        corresponding measurement is in the gate or not ---
            %        (number of measurements) x 1
            
            if isempty(z)
                z_ingate = z;
                meas_in_gate = false(0,1);
                return
            end
            zlength = size(z,2);
            meas_in_gate = false(zlength,1);
            
            nx = size(density_pred.x,1);
            
            S = measmodel.H(density_pred.x(:,end))*density_pred.P(end-nx+1:end,end-nx+1:end)*measmodel.H(density_pred.x(:,end))' + measmodel.R;
            %Make sure matrix S is positive definite
            S = (S+S')/2;
            
            nu = z - repmat(measmodel.h(density_pred.x(:,end)),[1 zlength]);
            dist = diag(nu.'*(S\nu));
            
            meas_in_gate(dist<gating_size) = true;
            z_ingate = z(:,meas_in_gate);
        end
        
        function density = momentMatching(w, mixture_density)
            %MOMENTMATCHING: approximate a Gaussian mixture density as a single Gaussian using moment matching
            %INPUT: w: normalised weight of Gaussian components in logarithm domain --- (number of Gaussians) x 1 vector
            %       states: structure array of size (number of Gaussian components x 1), each structure has two fields
            %               x: means of Gaussian components --- (variable dimension) x (number of Gaussians) matrix
            %               P: variances of Gaussian components --- (variable dimension) x (variable dimension) x (number of Gaussians) matrix
            %OUTPUT:state: a structure with two fields:
            %               x_hat: approximated mean --- (variable dimension) x 1 vector
            %               P_hat: approximated covariance --- (variable dimension) x (variable dimension) matrix
            
            if length(w) == 1
                density = mixture_density;
                return;
            end
            
            w = exp(w);
            %Moment matching
            density.x = [mixture_density(:).x]*w;
            numGaussian = length(w);
            density.P = zeros(size(mixture_density(1).P));
            for i = 1:numGaussian
                %Add spread of means
                x_diff = mixture_density(i).x - density.x;
                density.P = density.P + w(i).*(mixture_density(i).P + x_diff*x_diff');
            end
        end
        
    end
end

