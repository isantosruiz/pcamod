%PCAMOD  Data-driven process modeling for monitoring purposes using PCA.
%   Author:
%      2019, Ildeberto de los Santos Ruiz
%      idelossantos@ittg.edu.mx

classdef pcamod
    properties
        Z             % normalized training data matrix
        T             % transformed data matrix
        P             % eigenvectors
        lambda        % eigenvalues
        mu            % mean vector of data matrix
        sigma         % standard-deviation vector of data matrix
    end
    methods
        function obj = pcamod(X)
            [obj.Z,obj.mu,obj.sigma] = zscore(X);
            [obj.P,obj.T,obj.lambda] = pca(obj.Z);
        end
        function n = pdimension(obj,percentage)
            tev = cumsum(obj.lambda)/sum(obj.lambda);
            n = find(tev >= percentage/100,1);
        end
        function z = normalize(obj,x)
            z = (x - obj.mu')./obj.sigma';
        end
        function [zp,zr,T2,SPE] = map(obj,x,q,normalize)
            if nargin < 4
                normalize = true;
            end
            if normalize
                z = obj.normalize(x);
            else
                z = x;
            end
            Pp = obj.P(:,1:q);
            Pr = obj.P(:,q+1:end);
            Lq = diag(1./obj.lambda(1:q));
            zp = Pp*Pp'*z;
            zr = Pr*Pr'*z;
            T2 = z'*Pp*Lq*Pp'*z;
            SPE = zr'*zr;
        end
        function [T2,SPE] = indices(obj,q)
            m = size(obj.Z,1);
            T2 = zeros(m,1);
            SPE = zeros(m,1);
            for k = 1:m
                [~,~,T2(k),SPE(k)] = obj.map(obj.Z(k,:)',q,false);
            end
        end
        function [uT2,uSPE] = ucl(obj,q,alpha)
            m = size(obj.Z,1);
            uT2 = q*(m^2-1)/m/(m-q)*finv(alpha,q,m-q);
            [~,SPE] = obj.indices(q);
            m = mean(SPE);
            v = var(SPE);
            uSPE = v/2/m*chi2inv(alpha,2*m^2/v);
        end
    end
end