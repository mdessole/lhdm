%% compile
mex -v -R2017b lhdm.c -lmwblas -lmwlapack


%% generate matrix with given sv distribution
%% generate matrix with given sv distribution
m = 512;
n = 1024;
M = max(m,n);
N = floor(min(m,n)*9/10);
% generate random matrix with
% condition number = O(gap)
[U,~] = qr(randn(M));
[V,~] = qr(randn(M)); 
gap = 1e-9;
eigv = 2.^(1-[1:M]) +gap;
%eigv(1:floor(M/3)) = eigv(1:floor(M/3)) + 1e+0;
eigv(1:N) = eigv(1:N) + 0.01*(N-[1:N]);
A = U * diag(eigv) * V';
A = A(1:m,1:n);

%% generate sparse rhs >= 0
perc = 10; % percentage of nonzeros
cardS = max(1,floor(n*perc/100));
x = zeros(n,1);
step = floor(n/cardS);
i = 1;
while (i<=n)
    x(i) = rand(1,1);
    i=i+step;
end
b = A*x;

[m,n] = size(A);

%% set lhdm parameters
tau_c = 3;
tau_u=1;
tau_w=5;
k_max = 32;

thres = [0.1*tau_c 0.1*tau_w 0.1*tau_u]; 
blocksize = [k_max k_max];

%% solve nnls
[x_lhdm, iter] = lhdm(A,b,thres,blocksize);


fprintf("residual ||A*x_lhdm-b||_2 = %g \n", norm(A*x_lhdm-b));
fprintf("distance from optimum ||x - x_lhdm||_2 = %g \n", norm(x-x_lhdm));
fprintf("optimal support cardinality ||x||_0 = %d \n", length(find(abs(x)>eps)));
fprintf("computed support cardinality ||x_lhdm||_0 =%d \n", length(find(abs(x_lhdm)>eps)));