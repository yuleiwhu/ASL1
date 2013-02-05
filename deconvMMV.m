function [cbf,spectrum] = deconvMMV(data,M0,Phi,Amp)
% DECONVMMV [FunctionAbstract]
% ------------------------------------------------------------------
%    Details:
%       deconvMMV is a function to compute the CBF and AAT distribution via
%       the deconvolutional method.
%    ---------------
%    Inputs:
%    ---------------
%       asl     -   cPASL object, contains the diverse-ti data (M0 and
%                   perfusion data at different inversion time);
%       Phi     -   The dictionary (normalized)
%       Amp     -   The amplitude of each column of dictionary
%    ---------------
%    Outputs:
%    ---------------
%       cbf     -   CBF
%       spectrum-   The distribution of AAT
%    ---------------
%    Example:
%    ---------------
%       %% Loading data
%       load /Users/lyu/Documents/MATLAB/Data/DASL_YL/rtasl_slice_7.mat
% 
%       %% Dictionary
%       DicAAT = (1.2-0.7):0.05:2.2-0.05;
%       N = length(DicAAT); % dictionary width
%       T = 1.2:0.1:2.2;
%       M = length(T); % dictionary height
%       Phi = zeros(M,N);
%       Amp = zeros(N,1);
%       for i = 1:N
%           Phi(:,i) = asl.paslModel(DicAAT(i),'buxton',T); % f is set to 1 by default
%           Amp(i) = norm(Phi(:,i));
%           Phi(:,i) = Phi(:,i)/Amp(i);
%       end
% 
%       [cbf,spectrum] = deconvMMV(asl,Phi,Amp);
%       spect = cMRI(spectrum);
%       spect.showMRI;
% ------------------------------------------------------------------
% Created on Jan 29, 2013 by Lei Yu
%


[M,N] = size(Phi);
[h,w,M] = size(data);


spectrum = zeros(h,w,N);
count = zeros(h,w);
for i = 1:h
    for j = 1:w
        % ------------------------------------------------------------------
        [locy,locx,convImg] = voisin(i,j,w,h);
        % ------------------------------------------------------------------
        % Collect data
        Y = zeros(M,length(locy));
        for t=1:length(locy)
            Y(:,t) = data(locy(t),locx(t),:);
        end
        
        % ------------------------------------------------------------------
        % Solve mmv problem
        option.noiselevel = 'l';
        option.maxIter = 1000;
        [S,D] = SolveGroupMMV(Phi,Y,option);
        S = S*D;
%         S = TMSBL(Phi, Y,'PRUNE_GAMMA',1e-5,'PRINT',0);
        % ------------------------------------------------------------------
        
        for t = 1:length(locy)
            spectrum(locy(t),locx(t),:) = squeeze(spectrum(locy(t),locx(t),:)) + S(:,t);
            count(locy(t),locx(t)) = count(locy(t),locx(t)) + 1;% counting times that been computed
        end
        
        fprintf('The (%d,%d) position.\n',i,j);
    end
end

spectrum = spectrum./repmat(count,[1,1,N]);

% ------------------------------------------------------------------
% Scale the M0 and Amplitude of dictionary
%
for i = 1:N
    M0 = M0(:);
    m0 = mean(M0(M0>100));
    spectrum(:,:,i) = spectrum(:,:,i)/Amp(i)/m0;
end
cbf = sum(spectrum,3); 
end

function [locy,locx,convImg] = voisin(i,j,w,h)
% Find neighbors
orderN = 2;
tempNShape = ones(2*orderN+1);
tempPulse = zeros(w,h);tempPulse(i,j) = 1;
convImg = conv2(tempPulse,tempNShape,'smae');
[locy,locx] = find(convImg~=0);
end