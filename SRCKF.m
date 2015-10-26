classdef SRCKF
    
    properties
        
        Q;
        R;
        Qsqrt;
        Rsqrt;
        
        nx;
        
        xkk1;
        Skk1;
        
        xkk;
        Skk;
        
    end
    
    methods
        
        % Constructor
        function obj = SRCKF(Q, R, nx)
            obj.Q = Q;
            obj.R = R;
            obj.Qsqrt = sqrt(Q);
            obj.Rsqrt = sqrt(R);
            obj.nx = nx;
        end
        
        
        % Propagation
        function obj = Predict(obj, xkk, Skk)
            
            obj.xkk1 = xkk;
            
            [~, Skk1] = qr([Skk obj.Qsqrt]',0);
            
            obj.Skk1 = Skk1';
            
        end
        
        % Correction
        function obj = Update(obj, z)
            
            Rsqrt = obj.Rsqrt;
            xkk1 = obj.xkk1;
            Skk1 = obj.Skk1;
            
            %%%========================================================================
            %%% Genrate a set of Cubature Points
            %%%========================================================================
            
            nx = obj.nx; % state vector dimension
            
            nPts = 2*nx;   % No. of Cubature Points
            
            CPtArray = sqrt(nPts/2)*[eye(nx) -eye(nx)];
            
            %%%========================================================================
            
            Xi =  repmat(xkk1,1,nPts) + Skk1*CPtArray;
            
            Zi = MstEq(Xi);
            
            zkk1 = sum(Zi,2)/nPts;      % predicted Measurement
            
            X = (Xi-repmat(xkk1,1,nPts))/sqrt(nPts);
            
            Z = (Zi-repmat(zkk1,1,nPts))/sqrt(nPts);
            
            [~, Szz] = qr([Z Rsqrt]',0);
            
            Szz = Szz';                 % Square-root Innovations Covariance
            
            Pxz = X*Z';
            
            G = (Pxz/Szz')/Szz;         % Cubature Kalman Gain
            
            xkk = xkk1 + G*(z - zkk1);
            
            [~, Skk] = qr([(X - G*Z)  G*Rsqrt]',0);
            
            Skk = Skk';
            
            % Save objects
            obj.xkk = xkk;
            obj.Skk = Skk;
            
        end
    end
    
    
    
    methods(Static)
        
        function z = MstEq(x)
            
            z = [ 0.8*cos(x(1,:)) - 0.2*cos(x(1,:) + x(2,:)) ;
                0.8*sin(x(1,:)) - 0.2*sin(x(1,:) + x(2,:)) ];
        end
    end
    
end