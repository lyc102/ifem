function latexerrtable(N,err,formatstr)

%% 
% display error table in latex formate.

if ~exist('formatstr','var')
    formatstr = '%0.5e';
end
k = length(N);
ts = zeros(k,5); ts(:,3) = 38; ts = char(ts); td = char(92*ones(k,2));
n = size(err,2);
dispstr = num2str(N);
for i=1:n
    dispstr = [dispstr ts num2str(err(:,i),formatstr)];
end
display([dispstr ts(:,1:2) td]);