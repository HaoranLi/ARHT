# This is a function to calculate a consistent estimator of $p^{-1}tr(\Sigma^m)$ for any positive integer m,
# given the spectral of the sample covariance matrix $S_n$.
# Reference Lemma 1 of Bai, Chen, Yao (2010),
# ON ESTIMATION OF THE POPULATION SPECTRAL DISTRIBUTIONFROM A HIGH-DIMENSIONAL SAMPLE COVARIANCE MATRIX
# Australian & New Zealand Journal of Statistics.

# @emp_eig: eigenvalues of S_n.
#
# @gamma: p/n.
#
# @degree: the order of the moment.
#
# @return: the 1st, 2nd, ..., degree-th, order of moments of the distribution spectral of Sigma.

moments_PSD = function( eigenvalues,
                        n,
                        degree){
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
                stop("n must be numeric")
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

        if(!is.numeric(degree)){
                stop("degree must be numeric")
        }
        if(length(degree)!= 1){
                stop("length of degree should be 1.")
        }
        if(degree <1){
                stop("degree must be no less than 1.")
        }
        if(ceiling(degree) != degree){
                stop("degree must be integer.")
        }

        gamma = p/n
        if(p >=n){
                emp_dual = eigenvalues[1:n]
        }else{
                emp_dual = c(eigenvalues, rep(0,n-p))
        }

        emp_moments = colMeans(outer( emp_dual, 1:degree, FUN = '^'))

        pop_moments = numeric(degree)

        pop_moments[1] = mean(eigenvalues)

        if(degree > 1){
                # recursive formulas; cannot be parallel
                for( kk in 2:degree){
                        max_values = kk %/% (1:(kk-1)) # The maximum possible value of i1,..., i_(kk-1)

                        possible_values = mapply(seq, 0, max_values, SIMPLIFY = FALSE) # Possible values of i

                        partitions = t(expand.grid(possible_values)) # all possible partitions

                        valid_partitions = partitions[, colSums(partitions * (1:(kk-1))) == kk, drop = FALSE] # valid partitions

                        # Coefficients
                        fac1 = (gamma^(colSums(valid_partitions)))

                        fac2 = factorial(kk) / apply( factorial(valid_partitions), 2, prod)

                        fac3 = factorial(kk+1 - colSums(valid_partitions))

                        offset = sum( fac1 * fac2 /fac3 * apply( (pop_moments[1:(kk-1)])^valid_partitions, 2, prod ))

                        pop_moments[kk]= (emp_moments[kk] - offset)/gamma
                }
        }
        return(pop_moments)
}
