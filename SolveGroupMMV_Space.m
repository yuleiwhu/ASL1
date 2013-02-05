function [S,curve] = SolveGroupMMV_Space(G,Y,option)
% Function SolveGroupMMV is to solve the sparse vector from multiple
% measurements vector (MMV) with assumption that sparse vectors share same
% sparse pattern.
%   
%   G is the dictionary, with identity matrix at the end;
%   Y is the multiple measurements matrix
%   s is the length of signal
%   e is the length of measurements, i.e. the rows of matrix G
%
%   Parameters should be fixed before using:
%       lambda : the importance of sparsity (between 0 and 1)
%       gamma  : the weight of sparsity and structure
%       a2     : the importance of the error term
%   Please refer the following papers for details:
%   [1] A. Gramfort, D. Strohmeier, J. Haueisen, M. Hamalainen, and M. 
%       Kowalski, Functional brain imaging with M/EEG using structured 
%       sparsity in time-frequency dictionaries., in IPMI, 2011, vol. 22,
%       pp. 600-11.
%   [2] R. Jenatton, J. Mairal, G. Obozinski, and F. Bach, Proximal methods
%       for hierarchical sparse coding, arXiv preprint arXiv:1009.2139,
%       vol. 12, pp. 2297-2334, 2010.
%--------------------------------------------------------------------------
if nargin == 3
    if ~isfield(option,'noiselevel')    option.noiselevel = 'l';    end
    if ~isfield(option,'maxIter')       option.maxIter = 1000;      end
    if ~isfield(option,'display')       option.display = 0;         end
elseif nargin < 3
    option.noiselevel = 'l';
    option.maxIter = 1000;
    option.display = 0;
end

noiselevel = option.noiselevel;
maxIter = option.maxIter;
display = option.display;

[m,n] = size(G);
[h,w] = size(Y);
ls = n;
X = zeros(n,w);
S = G'*Y;%zeros(ls,w);
% avgS = zeros(size(S));
% a1S = ones(ls,1); a1S = repmat(a1S,1,w);

% Parameters to be fixed
lambda = 0.1; % Sparse factor, bigger the sparser
switch lower(noiselevel)
    case 's'
        gamma = 0.005;
    case 'm'
        gamma = 0.05;
    case 'l'
        gamma = max(S(:))*0.2;
end

% gamma = max(0.5,mean(std(Y,1,2))*1.02);
a1 = gamma*lambda; % regularization of code 's'
a3 = gamma*(1-lambda); % regularization of group
% disp([gamma]);
sig = max(abs(eig(G'*G)));
beta = 1/2/sig;
tau = 1;

scaleAug = ones(maxIter,1);%*log(1:maxIter)/log(maxIter);
% scaleAug = (maxIter:-1:1)/maxIter;

curve = zeros(maxIter,1);
for i = 1:maxIter
    
    X0 = [S];
    prox = X+beta*(G)'*(Y-G*X);% + a4*(avg-X);
    pS = prox(:,:);
    % Soft Threshold S Term
    Sg = pS - beta*repmat(a1,n,w).*scaleAug(i); % sparse structure
    Sg(Sg<0) = 0;
    pSg = 1-beta*repmat(a3,n,w).*scaleAug(i)./repmat(sqrt(sum(Sg.^2,2))/sqrt(w),1,w); % group structure
    pSg(pSg<0)=0;
    S = pSg.*Sg;

    % Form new estimates
    Xstar = [S];
    
    % Compute steps
    tau0 = tau; tau = (1+sqrt(1+4*tau^2))/2;
    X = Xstar + (tau0-1)/tau*(Xstar-X0);
    
    % stop 
    curve(i) = norm(X(:)-X0(:));
    if i > 3
        if abs(curve(i)-curve(i-1))/abs(curve(i-1)-curve(i-2)) < 1e-4
            break;
        end
    end
  
    if display
    pause(0.2);
    figure(1);subplot(1,2,1);imagesc(S);
    subplot(1,2,2),plot(gamma)
    disp(['Iterations:',num2str(i),' Residual: ',num2str(norm(Y-G*X)^2/2)]);
    end
end

% S = X(:,:);
end

function x = soft(xin,h)
    temp = abs(xin)-h;
    temp(temp<0)=0;
    x = sign(xin).*temp;
end



