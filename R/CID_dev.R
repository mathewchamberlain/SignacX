CID.GetDistMat <- function(edges, n)
{
  "%^%" <- function(A, n) {if(n == 1) A else A %*% (A %^% (n-1)) }
  # create adjacency matrix
  m <- new("ngTMatrix", 
           i = c(as.integer(edges$V1)-1L, as.integer(edges$V2)-1L), 
           j = c(as.integer(edges$V2)-1L, as.integer(edges$V1)-1L),
           Dim = as.integer(c(max(edges), max(edges))))
  
  dm = list("") # initialize distance matrix
  for (j in 1:n)
      dm[[j]] = j * m %^% j
return(dm)
}