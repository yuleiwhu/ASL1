% MAIN_DIVERSE_MULTI_VOXEL_SIMU [FunctionAbstract]
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
% Created on Jan 30, 2013 by Lei Yu
% 

clear all
close all

% ---------------
% The simulated data will be the two gaussian centered pulsive
% set two centers
step = 0.05;
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

% Data preparing
w = 5;h=5;
dataorg = zeros(h,w,length(perfusion));
for i = 1:length(perfusion)
    dataorg(:,:,i) = perfusion(i);
end
data = dataorg + randn(h,w,length(perfusion))*1.5; 
fprintf('Measurement SNR: %2f\n',snr(dataorg(:),data(:)))
data(data<0) = 0;
fprintf('Measurement SNR after zero threshold: %2f\n',snr(dataorg(:),data(:)))

% Dictionary
DicAAT = (1.2-0.7):step:(2.2-step);
N = length(DicAAT); % dictionary width
M = length(T); % dictionary height
Phi = zeros(M,N);
Amp = zeros(N,1);
for i = 1:N
    Phi(:,i) = asl.paslModel(DicAAT(i),'buxton',T); % f is set to 1 by default
    Amp(i) = norm(Phi(:,i));
    Phi(:,i) = Phi(:,i)/Amp(i);
end

%%
[cbf,spectrum] = deconvMMV(data,M0,Phi,Amp);
%% Show results
spect = cMRI(spectrum);
spect.showMRI;
cbf = cMRI(sum(spectrum(:,:,10:end-10),3));
cbf.showMRI;
%%
% [p] = [22,45];
[p] = [2,2];
etaest = squeeze(spectrum(p(1),p(2),:));
dat = squeeze(data(p(1),p(2),:));
cbff = sum(etaest);
etaest = etaest/cbff;
perfusion_est= asl.paslModel([etaest,DicAAT'],'buxton',T)*cbff*mean(M0(M0>100));
figure(111);
hold on
plot(T,perfusion_est,'r');
hold on 
plot(T,dat,'g');
hold off
figure(110)
hold on
plot(DicAAT,etaest,'r'); axis tight;
xlabel('AAT');
ylabel('Distribution of AAT');

%%
est_snr = 0;
for i = 1:w
    for j= 1:h
        etaest = squeeze(spectrum(i,j,:));
        etaest = etaest/sum(etaest);
        est_snr = est_snr + snr(eta,etaest);
    end
end
disp(est_snr/w/h);

