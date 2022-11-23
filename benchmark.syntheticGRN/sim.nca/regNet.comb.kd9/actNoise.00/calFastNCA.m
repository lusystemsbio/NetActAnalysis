%X = csvread('subnet1_data.xls');
%A = csvread('subnet1_top.xls');

X = csvread('net30tf.states.exp.csv');
A = csvread('net30tf.matrix.csv');

[m,n] = size(X); 
X = X(2:m, 2:n);
#[m,n] = size(X) 

[m,n] = size(A); 
A = A(2:m, 2:n);
#[m,n] = size(A) 

csvwrite('data.csv', X) 
csvwrite('top.csv', A)

[Ae, Se] = FastNCA(X, A);

% save data to the file:

fname_ae = 'Ae.csv'; 
fname_se = 'Se.csv';

#[m,n] = size(Ae)
#[m,n] = size(Se)

csvwrite(fname_ae, Ae) 
csvwrite(fname_se, Se)
