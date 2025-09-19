function [XX,ei111] = OSDR1(A,B)
% OSDR1 is an implementation of Algorithm 1 
% from "M. Ahmadnasab, Symmetric-diagonal reductions as preprocessing for symmetric positive definite
% generalized eigenvalue solvers, Journal of Mathematical Modeling, 11 (2023) 301-322.",  which is used with 
% the same name in the numerical experiment section of the paper. 
% Inputs: Square matrices A and B with the same size.
% Outputs: Eigenpair set of (A,B), namely, XX for the set of eigenvectors 
% and a diagonal matrix ei111 from the set of eigenvalues.
% Written by: Morad Ahmadnasab, 5/10/2022, University of Kurdistan, Sanandaj, Iran.
%%%%%%%%%%%

%%%% Step 1:

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
s1=diag(di);
u1=V(:,Id);

%%%%  Step 2 together with symmetrizing discussed in Section 5.3:

aa=u1'*A*u1;

%
n=length(A);
for k1=2:n
    for k2=1:k1-1
         aa(k2,k1)=aa(k1,k2);
    end
end

%%%%  Step 3: 

[X,ei111]=eig(aa,s1,'chol');

%%%%  Step 4:

XX=u1*X;
end






