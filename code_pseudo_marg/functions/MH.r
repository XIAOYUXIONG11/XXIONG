## Metropolis-Hastings transition operator 
MH = function(current.state, block.index)
{
  block.index = SAMPLES$BLOCK.INDEX

  current.position = c()
  for(i in 1:length(BLOCKING[[block.index]]$group))
    current.position = c(current.position, CURRENT.STATE[[SAMPLER$blocks[[block.index]]$group[i]]])

  for(iii in ls(CURRENT.STATE))
    PROPOSED.STATE[[iii]] = CURRENT.STATE[[iii]]

  proposed.position = SAMPLER$blocks[[block.index]]$L.chol.cov.proposal %*% rnorm(length(current.position)) + current.position

  for(i in 1:length(BLOCKING[[block.index]]$group))
    {
      PROPOSED.STATE[[SAMPLER$blocks[[block.index]]$group[i]]] = proposed.position[1:MODEL$SIZE.GROUPS[SAMPLER$blocks[[block.index]]$group[i]]]
      proposed.position = proposed.position[-c(1:MODEL$SIZE.GROUPS[SAMPLER$blocks[[block.index]]$group[i]])]
    }

  SAMPLER$blocks[[block.index]]$compute.utilities(PROPOSED.STATE)
  compute.logjointlik(PROPOSED.STATE)
  
  A = min(0, PROPOSED.STATE[[paste("LOGJOINTLIK.AR", block.index, sep="")]] - CURRENT.STATE[[paste("LOGJOINTLIK.AR", block.index, sep="")]])

  if(is.nan(A)) A = -Inf
  
  A
}
