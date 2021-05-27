
# Discretizes each value of vector Y into bins of length 1/N
# If y in [0,1], then this corresponds to discretizing Y into (N+1) values
# Returns a vector Y with the discretized values
# If the original value were discretized, they are returned unchanged
discretize_Y <- function(Y, N) {
  if (length(unique(Y)) > (N + 1)) {
    return(round(N * Y) / N)
  } else {
    return(Y)
  }
}
