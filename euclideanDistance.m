X = [14; 2; 50; 11; 8; 7]
Y = [14]

D = sqrt(sum(X.^2,2) - 2 * X*Y.' + sum(Y.^2,2).')