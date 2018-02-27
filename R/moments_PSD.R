#' Consistent estimators of high-order moments of the population spectral distribution for high-dimensional data
#' @description  The function calculates consistent estimators of moments of the spectral distribution
#' of the population covariance matrix given the spectral of the sample covariance matrix.
#' @keywords population spectral moments estimators
#' @name moments_PSD
#' @param eigenvalues all eigenvalues of the sample covariance matrix including 0's.
#' @param n degree of freedom of the sample covariance matrix.
#' @param mom_degree the maximum order of moments.
#' @return Estimators of moments from the first to the \code{mom_degree} -th order.
#' @references  Bai, Z., Chen, J., & Yao, J. (2010).
#'\emph{On estimation of the population spectral distribution from a high-dimensional sample covariance matrix.}
#' Australian & New Zealand Journal of Statistics, 52(4), 423-437.
#' @examples
#' set.seed(10086)
#' n = 400; p= 500
#' pop_eig = seq(10,1,length = p)
#' # Data with covariance matrix diag(pop_eig)
#' Z = matrix(rnorm(n*p),n,p)
#' X = Z %*% diag(sqrt(pop_eig))
#' raw_eig = svd(cov(X))$d
#' emp_eig = raw_eig[raw_eig>=0]
#' # Moments of population spectral distribution
#' colMeans(outer(pop_eig, 1:4, "^"))
#' # Estimators
#' moments_PSD(emp_eig, n-1, 4)
#' @export
moments_PSD = function( eigenvalues,
                        n,
                        mom_degree){
        if(!is.vector(eigenvalues, mode = "numeric")){
                stop("eigenvalues must be a numeric vector")
        }else{
                if( any(eigenvalues<0)){
                        stop("eigenvalues must be nonnegative.")
                }else{
                        eigenvalues = sort(eigenvalues, decreasing = TRUE)
                        p = length(eigenvalues)
                }
        }
        if(!is.numeric(n)){
                stop("n must be a numeric atomic")
        }
        if(length(n)!= 1){
                stop("length of n should be 1.")
        }
        if(n<=0){
                stop("n must be positive.")
        }
        if(ceiling(n) != n){
                stop("n must be integer.")
        }
        if(!is.numeric(mom_degree)){
                stop("mom_degree must be numeric")
        }
        if(length(mom_degree)!= 1){
                stop("length of mom_degree should be 1.")
        }
        if(mom_degree <1){
                stop("mom_degree must be no less than 1.")
        }
        if(ceiling(mom_degree) != mom_degree){
                stop("mom_degree must be integer.")
        }
        gamma = p/n
        if(p >= n){
                emp_dual = eigenvalues[1:n]
        }else{
                emp_dual = c(eigenvalues, rep(0,n-p))
        }
        emp_moments = colMeans(outer( emp_dual, 1:mom_degree, FUN = '^'))
        pop_moments = numeric(mom_degree)
        pop_moments[1] = mean(eigenvalues)

        if(mom_degree > 1){
                # recursive formulas; cannot be parallel
                for( kk in 2:mom_degree){
                        max_values = kk %/% (1:(kk-1)) # The maximum possible value of i1,..., i_(kk-1) because j * i_{j}<= kk
                        possible_values = mapply(seq, 0, max_values, SIMPLIFY = FALSE) # Possible values of i
                        partitions = t(expand.grid(possible_values)) # all possible partitions
                        valid_partitions = partitions[, colSums(partitions * (1:(kk-1))) == kk, drop = FALSE] # valid partitions
                        # Coefficients
                        fac1 = (gamma^(colSums(valid_partitions)))  # gamma^(i1+i2+i3+...+ij)
                        fac2 = factorial(kk) / apply(factorial(valid_partitions), 2, prod) # fac/fac3 is the phi in equation (8) of the reference
                        fac3 = factorial(kk+1 - colSums(valid_partitions))
                        offset = sum( fac1 * fac2 /fac3 * apply( (pop_moments[1:(kk-1)])^valid_partitions, 2, prod ))
                        pop_moments[kk]= (emp_moments[kk] - offset) / gamma
                }
        }
        return(pop_moments)
}
