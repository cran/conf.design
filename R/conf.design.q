conf.design <- function(G, p, block.name = "Blocks",
                        treatment.names =
                        if(length(nam <- dimnames(G)[[2]]) == 0)
                        paste("T", 1:ncol(G), sep = "") else nam)
{
  list.mat <- function(M, f)
    {
      l <- list()
      for(i in 1:ncol(M))
        l[[i]] <- f(M[, i])
      l
    }
  if(!is.matrix(G))
    G <- matrix(G, nrow = 1)
  nf <- ncol(G)
  D <- as.matrix(0:(p - 1))
  if(nf > 1)
    for(j in 2:nf) {
      E <- D
      D <- NULL
      for(i in 0:(p - 1))
        D <- rbind(D, cbind(E, i))
    }
  m <- (D %*% t(G)) %% p
  B <- do.call("paste", c(list.mat(m, format), list(sep = "")))
  D <- cbind(B, format(D))[sort.list(B),  ]
  D <- as.data.frame(list.mat(D, as.factor))
  names(D) <- c(block.name, treatment.names)
  class(D) <- c("design", class(D))
  D
}

conf.set <- function(G, p) {
  space <- function(G, p)
    {
#
# Generate all distinct linear combinations of the rows of G over GF(p)
#
      x <- 0:(p - 1)
      M <- as.matrix(x)
      k <- nrow(G)
      if(k > 1) {
        for(i in 2:k) {
          N <- NULL
          for(j in x)
            N <- rbind(N, cbind(M, j))
          M <- N
        }
      }
      M <- (M %*% G) %% p
#
# if the rows of G can be assumed linearly independent the rest can be omitted.
#
      m <- 0
      for(j in 1:ncol(M))
        m <- p * m + M[, j]
      M[!duplicated(m),  , drop = F]
    }
  S <- space(G, p)
  S[apply(S, 1, function(x)
          any(t <- x > 0) && x[t][1] == 1),  ]
}

direct.sum <- function(D1, D2, ..., tiebreak = letters) {
  l <- list(...)
  if(length(l))
    return(Recall(D1, Recall(D2, ..., tiebreak = tiebreak[-1]), 
                  tiebreak = tiebreak))
  E1 <- lapply(D1, function(x, n2)
               rep(x, rep(n2, length(x))), nrow(D2))
  E2 <- lapply(D2, function(x, n1)
               rep(x, n1), nrow(D1))
  D <- c(E1, E2)
  if(any(i <- duplicated(names(D))))
    names(D)[i] <- paste(names(D)[i], tiebreak[1], sep = "")
  D <- as.data.frame(D)
  class(D) <- c("design", class(D))
  D
}
factorize <- function(x, ...)
  UseMethod("factorize")

factorize.default <- function (n) {
  if (!is.numeric(n)) 
    stop("cannot factorize non-numeric arguments")
  if (length(n) > 1) {
    l <- list()
    for (i in seq(along = n))
      l[[i]] <- Recall(n[i])
    return(l)
  }
  if (n != round(n) || n < 2) 
    return(n)
  tab <- primes(n)
  fac <- numeric(0)
  while(length(tab <- tab[n %% tab == 0]) > 0) {
    n <- n/prod(tab)
    fac <- c(fac, tab)
  }
  sort(fac)
}

factorize.factor <- function(x, name = deparse(substitute(x)),
                             extension = letters, drop = T, sep = "") {
  llev <- factorize.default(length(levels(x)))
  if(length(llev) == 1)
    return(if(drop) x else {
      x <- design(x)
      names(x) <- name
      x
    }
           )
  D <- NULL
  for(i in llev) {
    E <- D
    D <- NULL
    for(j in 1:i)
      D <- rbind(D, cbind(E, j))
  }
  l <- list()
  for(i in seq(along = llev))
    l[[i]] <- factor(D[, i][x] - 1)
  l <- as.data.frame(l)
  names(l) <- paste(name, extension[1:length(llev)], sep = sep)
  class(l) <- c("design", class(l))
  l
}

join <- function(...) {
  m <- list(...)
  l <- list()
  for(i in seq(along = m))
    l <- c(l, if(is.list(k <- m[[i]])) k else list(k))
  l <- lapply(l, function(f)
              format(as.character(f)))
  as.factor(do.call("paste", c(l, list(sep = ""))))
}

rjoin <- function(..., part.name = "Part") {
  l <- lapply(list(...), as.data.frame)	# for some safety...
  bf <- factor(paste(part.name,
                     rep(1:length(l), sapply(l, nrow)), sep = ""))
  D <- as.data.frame(c(list(bf), as.list(do.call("rbind", l))))
  names(D) <- c(part.name, names(l[[1]]))
  class(D) <- c("design", class(D))
  D
}

primes <- function(n) {
# Find all primes less than n (or max(n) if length(n) > 1).
# Uses an obvious sieve method.  Nothing flash.
#
  if ((M2 <- max(n)) <= 1)
    return(numeric(0))
  x <- 1:M2
  x[1] <- 0
  p <- 1
  M <- floor(sqrt(M2))
  while((p <- p + 1) <= M)
    if(x[p] != 0)
      x[seq(p^2, n, p)] <- 0
  x[x > 0]
}



