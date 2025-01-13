# source: https://towardsdatascience.com/a-new-coefficient-of-correlation-64ae4f260310

xicor <- function(X, Y, ties = TRUE){
    n <- length(X)
    r <- rank(Y[order(X)], ties.method = "random")
    set.seed(42)
    if(ties){
        l <- rank(Y[order(X)], ties.method = "max")
        return( 1 - n*sum( abs(r[-1] - r[-n]) ) / (2*sum(l*(n - l))) )
    } else {
        return( 1 - 3 * sum( abs(r[-1] - r[-n]) ) / (n^2 - 1) )    
    }
}