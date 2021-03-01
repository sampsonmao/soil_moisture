# Prior
P = I - D %*% (F.inv %*% t(F.inv)) %*% (t(FF) %*% L.inv)
K = rho - I
r = as.matrix(dist(D[,1:2], upper=T, diag=T))
J = theta/phi^2 * r * rho
V.inv = t(L.inv) %*% L.inv
W1 = K %*% V.inv %*% P
W2 = J %*% V.inv %*% P
trW1 = sum(diag(W1))
trW2 = sum(diag(W2))
trW1sq = sum(diag(W1 %*% W1))
trW2sq = sum(diag(W2 %*% W2))
trW1W2 = sum(diag(W1 %*% W2))
ref.mat = matrix(c(n-k, trW1, trW2,
                   trW1, trW1sq, trW1W2,
                   trW2, trW1W2, trW2sq), ncol=3, nrow=3, byrow=T)
prior = determinant(ref.mat)
