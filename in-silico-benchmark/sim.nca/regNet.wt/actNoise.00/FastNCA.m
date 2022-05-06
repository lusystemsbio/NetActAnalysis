function [Ae,Se] = FastNCA(X,A)

%% FastNCA algorithm
%% Model: X = AS + G 
%% (X: the microarray data, A: connectivity matrix, S: TFA matrix, G: noise)
%% (A is sparse and satisfies the NCA criteria)
%% Problem: estimate A and S from X, with a known structure of A
%%
%% Ref: 
%% 1. CQ Chang, Z Ding, YS Hung, PCW Fung, "Fast Network Component Analysis
%% (FastNCA) for gene regulatory network reconstruction from microarray data,"
%% submitted to Bioinformatics, Oct 2007.
%% 2. CQ Chang, YS Hung, PCW Fung, Z Ding, "Network component analysis for 
%% blind source separation," In Proc. 2006 International Conference on
%% Communications, Circuits and Systems (ICCCAS2006), volume 1, pages 323-326,
%% Guilin, China, June 2006.
%% 3. CQ Chang, Z Ding, YS Hung, PCW Fung, "Fast Network Component Analysis
%% for gene regulation networks," In Proc. 2007 IEEE International Workshop
%% on Machine Learning for Signal Processing, Thessaloniki, Greece, Aug 2007.
%%
%% Input:
%% X: the data matrix, gene expression microarray data
%% A: the network topology matrix (Z0 in the Bioinformatics paper), the structure of A 
%% (zero or not) is used to partition the data matrix
%%
%% Output:
%% Ae: estimate of the connectivity matrix
%% Se: estimate of the TFA matrix
%%
%% Note: NCA criteria are checked before applying this algorithm

%% Coded: Apr 2006; revised: Apr 2007; documented: Oct 2007 
%% CQ Chang, 
%% Dept. of Electrical & Electronic Engineering, 
%% The University of Hong Kong
%% Email: cqchang@eee.hku.hk

M = size(A,2);
[U,S,V] = svd(X,'econ');
XM = U(:,1:M)*S(1:M,1:M)*V(:,1:M)';
U=U(:,1:M);	% M-dimensional signal subspace of the data matrix

Ae = A;
l=1;
for l=1:M
    U0 = U(find(A(:,l)==0),:);	% partition, Xr in Bioinformatics paper
    if size(U0,1) < M
        [UU,SS,VV] = svd(U0);
        t = VV(:,end);
    else  % modified for speed
        [UU,SS,VV] = svd(U0','econ');
        t = UU(:,end);
    end  %% get the M-th right singular vector of Xr (U0)

    a = U*t; 
    a = a .* (A(:,l)~=0);  %% get the estimate of l-th column of A
    break
    Ae(:,l) = a/sum(abs(a))*length(find(A(:,l)~=0)); %% scale the column
end

Se = Ae\XM;  %% get the estimate of the TFA matrix
