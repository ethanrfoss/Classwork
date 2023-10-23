
L = [1 2;7 8;21 -1];
P = L*L';
Z = [1 -2 3 -2;-1 -1 -1 -1;4 4 5 4;7 8 9 2];
X = Z*Z';

disp('Eigenvalues of P:');
eig(P)

disp('EigenValues of X:');
eig(X)

disp('Orthogonal of L:');
null(L')

disp('Orthogonal of Z:');
null(Z')

disp('Rank of P:');
rank(P)

disp('Rank of X:');
rank(X)

disp('Singular Values of P:');
svd(P)

disp('SIngular Values of X:');
svd(X)
