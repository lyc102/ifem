function MN = mytimes(M,N,varargin)
% M is a K-by-I-by-J array and N is a K-by-J-by-L array
% compute a K-by-I-by-L array MN such that MN(k,:,:) = M(k,:,:)*N(k,:,:)
% example for reshape and permute
% N = zeros(9,3,3);
% for i=1:9
%     N(i,:,:) = [1,2,3;4,5,6;7,8,9];
% end
% N =  reshape(permute(N,[2,1,3]),27,3) will give
% [1,2,3;4,5,6;7,8,9;1,2,3;4,5,6;7,8,9...]
% N = permute(reshape(N,3,9,3),[2,1,3]) change it back

CNP = size(M,1);
NV = size(N,3);

if ~isempty(varargin)
    Inum = varargin{1}(1); Jnum = varargin{1}(2);
else
    Inum = 3; Jnum = 3;
end

CI = repmat(repmat([1:Inum]',Jnum,1),CNP,1) + reshape(Inum*repmat(0:1:CNP-1,Inum*Jnum,1),[],1);
CJ = repmat(reshape(repmat([1:Jnum],Inum,1),[],1),CNP,1) + reshape(Jnum*repmat(0:1:CNP-1,Inum*Jnum,1),[],1);
CX = zeros(size(CI));
for j = 1:Jnum
    vj = (j-1)*Inum + (1:Inum)';
    ii = repmat(vj,CNP,1) + Inum*Jnum*reshape(repmat(0:1:CNP-1,Inum,1),[],1);
    CX(ii) = reshape(M(:,:,j)',[],1)';
end
Mspa = sparse(CI,CJ,CX,Inum*CNP,Jnum*CNP);

if NV == 1
    
    MN = reshape(Mspa*reshape(N',[],1),Inum,[])';
    
elseif NV >= 2
    
    MN = Mspa*reshape(permute(N,[2,1,3]),Inum*CNP,NV); % Inum maynot be correct
    
    MN = permute(reshape(MN,Inum,CNP,NV),[2,1,3]);
    
end

% CI = repmat(repmat([1;2;3],3,1),CNP,1) + reshape(3*repmat(0:1:CNP-1,9,1),[],1);
% CJ = repmat(reshape(repmat([1,2,3],3,1),[],1),CNP,1) + reshape(3*repmat(0:1:CNP-1,9,1),[],1);
% CX = zeros(size(CI));
% i1 = repmat([1;2;3],CNP,1) + 9*reshape(repmat(0:1:CNP-1,3,1),[],1);
% i2 = repmat([4;5;6],CNP,1) + 9*reshape(repmat(0:1:CNP-1,3,1),[],1);
% i3 = repmat([7;8;9],CNP,1) + 9*reshape(repmat(0:1:CNP-1,3,1),[],1);
% CX(i1) = reshape(M(:,:,1)',[],1)';
% CX(i2) = reshape(M(:,:,2)',[],1)';
% CX(i3) = reshape(M(:,:,3)',[],1)';
% Mspa = sparse(CI,CJ,CX,3*CNP,3*CNP);
% if NV == 1
%     
%     MN = reshape(Mspa*reshape(N',[],1),3,[])';
%     
% elseif NV >= 2
%     
%     MN = Mspa*reshape(permute(N,[2,1,3]),3*CNP,NV);
%     
%     MN = permute(reshape(MN,3,CNP,NV),[2,1,3]);
%     
% end


end