## **************************************************************************************************** 
## ************************************************** Laplace Approximation - Newton-Raphson
## **************************************************************************************************** 

## Putting the break before updating "converged" ensures that another iteration is performed before leaving the loop, so that W, B, L are updated and computed according to h.vect
LA.sparse = function(STATE)
  {
    h.vect = rep(0, DATA$n+1)

    converged = F
    for(iiii in 1:1000)
      {
        grad.wrt.f = - 0.5 * DATA$y^2 * exp(- (h.vect[1:DATA$n] + h.vect[DATA$n+1]))
        
        diag.neg.hess.p.y.giv.h = -c(grad.wrt.f, sum(grad.wrt.f))
        
        neg.hess = STATE$Q.mat.sparse
        
        diag(neg.hess) = diag(neg.hess) + diag.neg.hess.p.y.giv.h
        neg.hess[1:DATA$n, DATA$n+1] = neg.hess[DATA$n+1, 1:DATA$n] = -grad.wrt.f
        if(converged == T) break
        
        grad.wrt.h = -0.5 - grad.wrt.f
        grad.wrt.h = c(grad.wrt.h, sum(grad.wrt.h)) - STATE$Q.mat.sparse %*% h.vect
        
        h.new = h.vect + solve(neg.hess, grad.wrt.h)

        if(mean((h.new - h.vect)^2) < 1e-4) converged = T

        h.vect = h.new
      }

    res = list()
    res$mu = h.vect
    res$P = neg.hess
    
    res
  }
