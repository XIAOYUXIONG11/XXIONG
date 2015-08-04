## Elliptical Slice sampling transition operator - for f only! 
ELLIPTICALSS = function(current.state, block.index)
{
  block.index = SAMPLES$BLOCK.INDEX

  current.position = c()
  for(i in 1:length(BLOCKING[[block.index]]$group))
    current.position = c(current.position, CURRENT.STATE[[SAMPLER$blocks[[block.index]]$group[i]]])

  for(iii in ls(CURRENT.STATE))
    PROPOSED.STATE[[iii]] = CURRENT.STATE[[iii]]

  auxiliary.nu = exp(CURRENT.STATE$psi.sigma/2) * CURRENT.STATE$L.chol.Q.mat %*% rnorm(DATA$n)

  logy = CURRENT.STATE$log.p.y.giv.f + log(runif(1))
  auxiliary.theta = runif(1, 0, 2*pi)
  auxiliary.thetamin = auxiliary.theta - 2*pi
  auxiliary.thetamax = auxiliary.theta

  while(1) {
    PROPOSED.STATE$f = CURRENT.STATE$f * cos(auxiliary.theta) + auxiliary.nu * sin(auxiliary.theta)
    
    SAMPLER$blocks[[block.index]]$compute.utilities(PROPOSED.STATE)
    compute.logjointlik(PROPOSED.STATE)

    if(PROPOSED.STATE$log.p.y.giv.f > logy)
      break
    
    if(auxiliary.theta < 0) auxiliary.thetamin = auxiliary.theta
    if(auxiliary.theta >= 0) auxiliary.thetamax = auxiliary.theta      
    auxiliary.theta = runif(1, auxiliary.thetamin, auxiliary.thetamax)
  }
  
  A = 0
  
  A
}
