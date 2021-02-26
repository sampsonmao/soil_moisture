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




$$
  \begin{aligned}
V &= LL^T \\
F &= L^{-1}D \rightarrow F^TF = DV^{-1}D \\
F &= QR \rightarrow F^{-1} = Q^TR^{-1} \rightarrow RF^{-1} = Q^T \\
F^{-1}F^{-T} &= (DV^{-1}D)^{-1} \\
Z &= L^{-1}I \rightarrow LZ = I \rightarrow Z=L^{-1} \\
F^TL^{-1} &= D^TL^{-T}L^{-1} = D^TV^{-1} = F^TZ 
\end{aligned}
$$
  $$
  \begin{aligned}
V(\theta, \phi) &=\theta\rho(\phi) + (1-\theta)I \\
P(\theta,\phi) &= I-DF^{-1}F^{-T}F^TZ\\
W_1 &= \frac{\partial V(\theta, \phi)}{\partial \theta} V(\theta, \phi)^{-1}P(\theta,\phi)
=KV^{-1}P \\
W_2 &= \frac{\partial V(\theta, \phi)}{\partial \phi} V(\theta, \phi)^{-1}P(\theta,\phi)
=JV^{-1}P \\
\pi^R(\theta,\phi) &\propto 
\begin{vmatrix}
n-k & trW_1 & trW_2 \\
trW_1 & trW_1^2 & trW_1W_2 \\
trW_2 & trW_1W2 & trW_2^2
\end{vmatrix}^{\frac{1}{2}}
\end{aligned}
$$
  
  $$
  \begin{aligned}
K(\phi) &= \rho(\phi) - I \\
J(\theta, \phi) &= \theta \frac{r}{\phi^2}\rho(\phi)
\end{aligned}
$$
  
  $$
  \begin{aligned}
V(\theta, \phi) &=\theta\rho(\phi) + (1-\theta)I \\
P(\theta,\phi) &= I-D(D^TV(\theta,\phi)^{-1}D)^{-1}D^TV(\theta,\phi)^{-1}\\
W_1 &= \frac{\partial V(\theta, \phi)}{\partial \theta} V(\theta, \phi)^{-1}P(\theta,\phi)
=\rho(\phi) -I\\
W_2 &= \frac{\partial V(\theta, \phi)}{\partial \phi} V(\theta, \phi)^{-1}P(\theta,\phi)
=\frac{r}{\phi^2}\rho(\phi) \\
\pi^R(\theta,\phi) &\propto 
\begin{vmatrix}
n-k & trW_1 & trW_2 \\
trW_1 & trW_1^2 & trW_1W_2 \\
trW_2 & trW_1W2 & trW_2^2
\end{vmatrix}^{\frac{1}{2}}
\end{aligned}
$$