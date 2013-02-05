function [D,Amp] = normcol(Phi)
% NORMCOL [FunctionAbstract]
% ------------------------------------------------------------------
%    Details:
%       Normalize the columns of matrix.
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
% Created on Feb 04, 2013 by Lei Yu
% 
D = zeros(size(Phi));
Amp = zeros(size(Phi,2),1);
for i = 1:size(Phi,2)
    Amp(i) = norm(Phi(:,i));
    D(:,i) = Phi(:,i) / (Amp(i)+eps);
end
end