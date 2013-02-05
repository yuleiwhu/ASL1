function [X] = SolveGroupMMV_CoSaMP(Phi,Y,option)
if nargin == 3
    if ~isfield(option,'K')
        option.K = floor(size(Phi,1)/2-1);
    end
elseif nargin<3
    option.K = floor(size(Phi,1)/2-1);
end
K = option.K;

[ M,N] = size(Phi);
[M,m] = size(Y);
Omega = [];
v = Y;
for i = 1:30
    Prox_s = Phi'*v;
    prox_rows = sum(Prox_s.^2,2);
    
    [prox_order,ind] = sort(prox_rows,'descend');
    omega = ind(1:2*K);
    
    Omega = union(omega,Omega);
    sP = Phi(:,Omega);
    X = zeros(N,m);
    Prox_s = inv(sP'*sP+eye(length(Omega))*0.0001)*sP'*Y;
    prox_rows = sum(Prox_s.^2,2);
    X(Omega,:) = Prox_s;
    [prox_order,ind] = sort(prox_rows,'descend');
    
    X(Omega(ind(K+1:end)),:)=0;
    X(X<0) = 0;
    Omega(ind(K+1:end)) = [];
    
    v = Y - Phi(:,Omega)*X(Omega,:);
    
    pause(0.2);
    figure(1);subplot(1,2,1);imagesc(X);
    %     subplot(1,2,2),plot(gamma)
    disp(['Iterations:',num2str(i),' Residual: ',num2str(norm(Y-Phi*X)^2/2)]);
end

end