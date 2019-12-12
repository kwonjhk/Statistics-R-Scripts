#occupancy matrix of DTMC

omatrix <- function(x, n) {
  r = dim(x)[1]
  c = dim(x)[2]
  final = matrix(0L, nrow=r, ncol=c)
  for(i in 0:n) {
    y = diag(1, nrow=r, ncol=r)
    if(i>0) {
      for(j in 1:i) {
        y = y %*% x
      }
    }
    print(y)
    final = final + y
  }
  return(final)
}

#calculate p.hat matrix of CTMC given rate matrix matrix.r

p.hat <- function(matrix.r) {
  r = max(rowSums(matrix.r))
  row = dim(matrix.r)[1]
  col = dim(matrix.r)[2]
  matrix.p = matrix(0L, nrow=row, ncol=col)
  for(i in 1:row) {
    for(j in 1:col) {
      if(i==j) {
        matrix.p[i,j] = (r-rowSums(matrix.r)[i])/r
      } else {
        matrix.p[i,j] = matrix.r[i,j]/r
      }
    }
  }
  return(matrix.p)
}

x.test = rbind(c(0, 2, 3, 0),
            c(4, 0, 2, 0),
            c(0, 2, 0, 2),
            c(1, 0, 3, 0))
print(p.hat(x.test))

#transient probability matrix given x with rt = rmax.t

tmatrix.c <- function(x, n, rmax.t) {
  r = dim(x)[1]
  c = dim(x)[2]
  final = matrix(0L, nrow=r, ncol=c)
  for(i in 0:n) {
    y = diag(1, nrow=r, ncol=r)
    if(i>0) {
      for(j in 1:i) {
        y = y %*% x
      }
    }
    final = final + y*(exp(-rmax.t)*(rmax.t^i)/factorial(i))
  }
  print(final)
}

#occupancy matrix of CTMC
omatrix.c <- function(x, m, r.1, t) {
  r = dim(x)[1]
  c = dim(x)[2]
  rmax.t = r.1*t
  final = matrix(0L, nrow=r, ncol=c)
  for(i in 0:m) {
    y = diag(1, nrow=r, ncol=r)
    if(i>0) {
      for(j in 1:i) {
        y = y %*% x
      }
    }
    final = final + (y*ppois(i, lambda=rmax.t, lower=FALSE))
  }
  print(final/r.1)
}
