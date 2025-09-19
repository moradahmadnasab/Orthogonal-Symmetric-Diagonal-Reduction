function [XX,ei111] = RSSEP1(A,B)
% RSSEP1 is an implementation of Algorithm 2 
% from "M. Ahmadnasab, Symmetric-diagonal reductions as preprocessing for symmetric positive definite
% generalized eigenvalue solvers, Journal of Mathematical Modeling, 11 (2023) 301-322.", which is used with 
% the same name in the numerical experiment section of the paper.   
% Inputs: Square matrices A and B with the same size.
% Outputs: Eigenpair set of (A,B), namely, XX for eigenvectors and 
% a diagonal matrix ei111 from the set of eigenvalues.
% Written by: Morad Ahmadnasab, 5/10/2022, University of Kurdistan, Sanandaj, Iran.
%%%%%%%%%%%

%%%% Step 1:

n=length(B);
[V,D]=eig(B);
d=diag(D);
maxd=max(d);
mind=min(d);
maxdmind=maxd/mind;
if abs(maxdmind) <1e+10
   [di,Id]=sort(d,'descend');
else
   [di,Id]=sort(d,'ascend');
end

u1=V(:,Id);

%%%% Step 2 together with symmetrizing discussed in Section 5.3:

s1=diag(sqrt(di));  
s2=zeros([n n]);
for i = 1:n
    s2(i,i)=1/s1(i,i);
end

aa=(s2*u1')*A*(u1*s2);

%
for k1=2:n
    for k2=1:k1-1
        aa(k2,k1)=aa(k1,k2);
    end
end


%%%% Step 3:

[X,ei111]=eig(aa);

%%%% Step 4: 

XX=(u1*s2)*X;

end





