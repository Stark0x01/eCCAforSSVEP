function [output1, output2] = CCA(signal1, signal2)
% Canonical Correlation Analysis, CCA
% Input:
%   signal: #channels, #points
% Output:
% output1: spatial filter of signal1
% output2: spatial filter of signal2
%
X=signal1;
Y=signal2;
T=size(signal2, 2);
%compute covariance matrix 
meanx=mean(X,2);
meany=mean(Y,2);
s11=0;s22=0;s12=0;s21=0;
for i1=1:T
  s11=s11+(X(:,i1)-meanx)*(X(:,i1)-meanx)';
  s22=s22+(Y(:,i1)-meany)*(Y(:,i1)-meany)';
  s12=s12+(X(:,i1)-meanx)*(Y(:,i1)-meany)';
  s21=s21+(Y(:,i1)-meany)*(X(:,i1)-meanx)';
end
s11=s11/(T-1);         
s22=s22/(T-1); 
s12=s12/(T-1); 
s21=s21/(T-1);
%compute eigvalue and eigvector
[eigvectora,eigvaluea]=eig(inv(s11)*s12*inv(s22)*s21);
[eigvectorb,eigvalueb]=eig(inv(s22)*s21*inv(s11)*s12);

evaluea= diag(eigvaluea);
evalueb = diag(eigvalueb);
% correlation coefficient & canonical variates of signal1
[corrcoef1, index]= max(sqrt(evaluea));
output1 = eigvectora(:, index);
% correlation coefficient & canonical variates of signal2
[corrcoef2, index]= max(sqrt(evalueb));
output2 = eigvectorb(:, index);    
end

