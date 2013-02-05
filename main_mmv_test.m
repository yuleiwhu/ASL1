% MAIN_MMV_TEST [FunctionAbstract]
% ------------------------------------------------------------------
%    Details:
%    ---------------
%    Inputs:
%    ---------------
% 
%    ---------------
%    Outputs:
%    ---------------
%
%    ---------------
%    Example:
%    ---------------
% 
% ------------------------------------------------------------------
% Created on Jan 31, 2013 by Lei Yu
% 


load Phi_0_02_buxton.mat
% Phi= randn(size(Phi));
% Phi = normcol(Phi);

step = 0.02;
mu1 = 1; mu2 = 1.9;
daat = (1.2-0.7):step:(2.2-step);
rho = 0.999;
eta1 = pdf('norm',daat,mu1,0.1)*rho;
eta2 = pdf('norm',daat,mu2,0.01)*(1-rho);
eta = (eta1+eta2)/sum(eta1+eta2);

y = Phi*eta(:);

m = 10;
Yo = repmat(y,1,m);
Y = Yo + randn(size(Yo))*.2;
% .2 then snr = 2-4dB
% .1 then snr = 8-10dB
% .05 then snr = 13-14dB
disp(snr(Yo,Y))
%%
option.noiselevel = 'l';
option.maxIter = 500;
option.display = 0;
option.K = 2;
X = Phi'*Y;
% [X1] = SolveGroupMMV_CoSaMP(Phi,Y,option);
tic;[X1,D,tau0] = SolveGroupMMV(Phi,Y,option);toc
X1 = D*X1;
disp(norm(Y-Phi*X1));

xeta1 = X1(:,5);
disp(sum(xeta1))

xeta1 = xeta1/sum(xeta1);

% for i = 1:m
% X(:,i) = SolveCoSaMP(Phi,Y(:,i),1);
% end
tic
X = TMSBL(Phi, Y,'PRUNE_GAMMA',1e-3,'PRINT',0);
toc
disp(norm(Y-Phi*X));


figure(1000),plot(tau0);
xeta = X(:,5);%mean(X,2);
disp(sum(xeta));
xeta = xeta/sum(xeta);


figure(11),clf;
hold on;
plot(daat,eta)
plot(daat,xeta,'r')
plot(daat,xeta1,'.-g');
hold off;
title([num2str(sum(xeta1(:).*daat(:))),', ',num2str(sum(eta(:).*daat(:)))])
return;
%%


clear all
close all

% ---------------
% The simulated data will be the two gaussian centered pulsive
% set two centers
step = 0.02;
mu1 = 1.5; mu2 = 1;
daat = (1.2-0.7):step:(2.2-step);
rho = 0.8;
eta1 = pdf('norm',daat,mu1,0.05)*rho;
eta2 = pdf('norm',daat,mu2,0.05)*(1-rho);
eta = (eta1+eta2)/sum(eta1+eta2);

T = 1.2:0.1:2.2;
asl = cPASL();
ccbf = 100;
M0 = 800;
perfusion = asl.paslModel([eta(:),daat(:)],'buxton',T)*M0*ccbf;

figure(110),plot(daat,eta);
figure(111),plot(T,perfusion);
%%
% --------------
% Data prepare
m = 10;
Yorg = repmat(perfusion,1,m);
sigma = 0;
Y = Yorg + randn(size(Yorg))*sigma;


% Dictionary
DicAAT = (1.2-0.7):step:(2.2-step);
N = length(DicAAT); % dictionary width
M = length(T); % dictionary height
Phi = zeros(M,N);
Amp = zeros(N,1);
for i = 1:N
    Phi(:,i) = asl.paslModel(DicAAT(i),'norm',T); % f is set to 1 by default
    Amp(i) = norm(Phi(:,i));
    Phi(:,i) = Phi(:,i)/Amp(i);
    disp(i)
end
%%
% -------------
%
S = SolveGroupMMV(Phi,Y);
spectrum = sum(S,2)/m;
for i = 1:N
    M0 = M0(:);
    m0 = mean(M0(M0>100));
    spectrum(i) = spectrum(i)/Amp(i)/m0;
end

%%
etaest = spectrum;
dat1 = perfusion;
dat2 = mean(Y,2);
cbff = sum(etaest);
etaest = etaest/cbff;
perfusion_est= asl.paslModel([etaest,DicAAT'],'buxton',T)*cbff*mean(M0(M0>100));
figure(111);
hold on
plot(T,perfusion_est,'r');
hold on 
plot(T,dat1,'g');
plot(T,dat2,'b');
hold off
figure(110)
hold on
plot(DicAAT,etaest,'r'); axis tight;
xlabel('AAT');
ylabel('Distribution of AAT');
