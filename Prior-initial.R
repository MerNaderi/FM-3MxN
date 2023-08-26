
Prior.initial = function(data, method = c("kmeans", "fast-kmed",
                                          "inc-kmed", "rank-kmed",
                                          "Simple-kmed"), G = G ){

  if(is.na(dim(data)[3])){
    stop("data should be an array with dimension dimension p * n * N")
  }

  p = dim(data)[1]
  n = dim(data)[2]
  N = dim(data)[3]

  vect = matrix(0, N, p * n)
  for(i in 1:N) vect[i, ] = c(data[,, i])
  CT = data[,,1]
  for(j in 2:N) CT = cbind(CT, data[,, j])

  RT = data[,,1]
  for(j in 2:N) RT = rbind(RT, data[,, j])

  if(method == "kmeans"){
    KMEAN = kmeans(vect, G)
    Group = KMEAN$cluster
  }

  if(method == "fast-kmed"){
    mrwdist = kmed :: distNumeric(vect, vect, method = "mrw")
    result = kmed :: fastkmed(mrwdist, ncluster = G, iterate = 50)
    Group = result$cluster
  }
  if(method == "inc-kmed"){
    mrwdist = kmed :: distNumeric(vect, vect, method = "mrw")
    result = kmed :: inckmed(mrwdist, ncluster = G, iterate = 50, alpha = 1.5)
    Group = result$cluster
  }
  if(method == "rank-kmed"){
    mrwdist = kmed :: distNumeric(vect, vect, method = "mrw")
    result = kmed :: rankkmed(mrwdist, ncluster = G, iterate = 50)
    Group = result$cluster
  }
  if(method == "Simple-kmed"){
    mrwdist = kmed :: distNumeric(vect, vect, method = "mrw")
    result = kmed :: skm(mrwdist, ncluster = G, seeding = 50)
    Group = result$cluster
  }

  nu0 = rep(4, G)
  pp = table(Group)/N

  M0 = L0 = array(0, c(p, n, G))
  S0 = array(0, c(p, p, G))
  Ps0 = array(0, c(n, n, G))
  for(i in 1:G){
    M0[,, i] = apply(data[,, Group == i], c(1,2), mean)
    L0[,, i] = apply(data[,, Group == i], c(1,2), moments::skewness)
    Ps0[,, i] = cov(RT[Group == i, ])
    S0[,, i] = cov(t(CT[, Group == i]))
  }

  out = list(M = M0, L= L0, S = S0, Ps = Ps0, nu = nu0, pp = pp)

  return(out)
}







