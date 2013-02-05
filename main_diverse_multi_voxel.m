% MAIN_DIVERSE_MULTI_VOXEL [FunctionAbstract]
% ------------------------------------------------------------------
%    Details: 
%    Compute the distribution of AAT from real data
%    ---------------
%    Parameters:
%    ---------------
%       Some parameters should be given before running:
%           numSlice    -   index of slice
%    ---------------
%    Outputs:
%    ---------------
%           spect       -   DAAT with format of cMRI
%           cbf         -   CBF with format of cMRI
%
%       After running, you can use, 
%           spect.showMRI to check the daat
%           cbf.showMRI to check the cbf
%
%    ---------------
%    Example:
%    ---------------
%   
% ------------------------------------------------------------------
% Created on Jan 30, 2013 by Lei Yu
% 

close all
clear all
% Loading data
numSlice = 7;
step = 0.05;
[path2file] = dataFetch(numSlice);
load(path2file);

perfusion = asl.dMRI;
m0 = asl.dM0;
w = perfusion(1).wImg;
h = perfusion(1).hImg;
data = zeros(h,w,length(perfusion));
dM0 = zeros(h,w,length(perfusion));
for i = 1:length(perfusion)
    data(:,:,i) = perfusion(i).dMRI;
    dM0(:,:,i) = m0(i).dMRI;
end
M0 = mean(dM0,3);

% Dictionary


T = 1.2:0.1:2.2;
T = T + (numSlice-1)*0.045;
M = length(T); % dictionary height
DicAAT = (min(T)-0.7):step:(max(T)-0.1);
N = length(DicAAT); % dictionary width
Phi = zeros(M,N);
Amp = zeros(N,1);
for i = 1:N
    Phi(:,i) = asl.paslModel(DicAAT(i),'buxton',T); % f is set to 1 by default
    Amp(i) = norm(Phi(:,i));
    Phi(:,i) = Phi(:,i)/Amp(i);
end

[cbf,spectrum] = deconvMMV(data,M0,Phi,Amp);
spect = cMRI(spectrum);
spect.showMRI;
cbf = cMRI(sum(spectrum(:,:,10:end-10),3));
cbf.showMRI;
save([path2file(1:end-4),'_results_TMSBL.mat'],'cbf','spectrum','spect','M0','T','DicAAT','data');
return;
%%
[p] = [22,40];
% [p] = [30,30];
daat = squeeze(spectrum(p(1),p(2),:));
dat = squeeze(data(p(1),p(2),:));
cbff = sum(daat);
daat = daat/cbff;
perfusion= asl.paslModel([daat,DicAAT'],'buxton',T)*cbff*mean(M0(M0>100));
figure(100);
plot(T,perfusion);
hold on 
plot(T,dat,'r');
hold off
figure(101)
plot(DicAAT,daat); axis tight;
xlabel('AAT');
ylabel('Distribution of AAT');

return;
%%
AAAT = zeros(h,w);
aat_spect = sum(spectrum,3);
for i = 1:N
    c_spect = spectrum(:,:,i);
    AAAT(aat_spect>5) = AAAT(aat_spect>5) +  c_spect(aat_spect>5) ./ aat_spect(aat_spect>5) .* DicAAT(i);
end
AAT = cMRI(AAAT);
AAT.showMRI;
return;