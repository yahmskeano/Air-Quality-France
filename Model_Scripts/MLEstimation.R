exponential <- function(a, b, dist){
  exp = a * exp(-(dist/b)^2)
  return(dist)
}

sep <- function(spat.a, spat.b, spat.dist, temp.a, temp.b, temp.dist){
  temp = exponential(temp.a, temp.b, temp.dist)
  spat = exponential(spat.a, spat.b, spat.dist)
  cov = temp * spat
  return(cov)
}

get_mu <- function(X, cov, Z){
  cz_inv = solve(cov)
  beta = solve(t(X) %*% cz_inv %*% X) %*% t(X) %*% cz_inv %*% Z
  return(X %*% beta)
}



gp_ll <- function(vec){
  spat.a = vec[1]
  spat.b = vec[2]
  temp.a = vec[3]
  temp.b = vec[4]
  
  cov = sep(spat.a, spat.b, dist, temp.a, temp.b, dist)
  
  mu = get_mu(X, cov, Z)
  
  ll = dmvnorm(Z, mu, cov, log=TRUE)
  
  return(ll)
}