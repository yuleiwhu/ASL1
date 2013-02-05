function [path2file] = dataFetch(numSlice)
% MAIN_DATAFETCH [FunctionAbstract]
% ------------------------------------------------------------------
%    Details:
%       Read data from nii files, and save the data slice by slice with
%       different inversion time.
%
%       The data is saved with the format of cPASL, so if you don't know
%       the format, please go to see the documents of cPASL.
%
%       In cPASL, the field dMRI is the place where perfusion data
%       (difference between labeled and control images) is located, and the
%       field dM0 is the place where the M0 (the first repetition of the
%       data) is located.
%
%       If you only want to fecth data of a fixed slice, please set the
%       parameter numSlice to the index of the slice.
%
% see also
%    cPASL
% ------------------------------------------------------------------
% Created on Jan 29, 2013 by Lei Yu
%
dataPath = '/Users/lyu/Documents/MATLAB/Data/DASL_YL/';
prefix = 'rtasl_';
ext = '.nii';
range = [1200:100:2200];
numTI = length(range);

% Choose another slice
if nargin < 1
    numSlice = 7;
end
path2file = [dataPath,prefix,'slice_',num2str(numSlice),'.mat'];
if exist(path2file,'file')
    reply = lower(input('File already exists. Do you want to continue? Y/N [Y]:','s'));
    if isempty(reply)
        return;
    elseif isequal(reply,'n')
        return;
    end
end


asl = cPASL();
asl.dMRI(numTI) = cMRI;
asl.dM0(numTI) = cMRI;

for i = 1:numTI
    fullPath = [dataPath,prefix,num2str(range(i)),ext];
    obj = cMRI(fullPath);
    control = cMRI(obj.dMRI(:,:,:,2:2:end));
    label = cMRI(obj.dMRI(:,:,:,3:2:end));
    perfusion = label-control;
    objAvg = perfusion.mean(numSlice);
    
    asl.dMRI(i) = objAvg;
    asl.dM0(i) = cMRI(obj.dMRI(:,:,numSlice,1));
end



save(path2file,'asl');

end