#' 3-variate positively correlated chi-squared sample generation when degrees of freedom are large

#' @export
r3chisq = function(size, df, corr_mat){

        if(!is.numeric(size)){
                stop("size must be numeric")
        }
        if(length(size) != 1L){
                stop("length(size) must be 1")
        }
        if(!is.vector(df, mode = "numeric")){
                stop("df must be a numeric vector")
        }
        if(length(df) != 3L){
                stop("length(df) must be 3")
        }
        if(any(df < 0)){
                stop("df must be positive")
        }
        df = ceiling(df)

        if( (length(dim(corr_mat)) > 2L) || !(is.numeric(corr_mat))){
                stop("corr_mat must be a numeric matrix")
        }
        if(!is.matrix(corr_mat)){
                corr_mat = as.matrix(corr_mat)
        }
        if(!isSymmetric(corr_mat)){
                stop("corr_mat must be symmetric")
        }
        corr_mat = corr_mat * (corr_mat >= 0)

        half_covariances = floor(c((df[1])^(1/2) * (df[2])^(1/2) * corr_mat[1,2],
                                   (df[1])^(1/2) * (df[3])^(1/2) * corr_mat[1,3],
                                   (df[2])^(1/2) * (df[3])^(1/2) * corr_mat[2,3]))

        # Build three diag matrices Q1, Q2, Q3.  z^TQ_1z, z^TQ_2z, z^TQ_3z would be the chi-square vector
        # See supplementary material for algorithm explanation
        mincov = min(half_covariances)
        part1 = c(rep(1, mincov), rep(1, half_covariances[1]-mincov), rep(1, half_covariances[2]-mincov), rep(0, half_covariances[3]-mincov))
        part2 = c(rep(1, mincov), rep(1, half_covariances[1]-mincov), rep(0, half_covariances[2]-mincov), rep(1, half_covariances[3]-mincov))
        part3 = c(rep(1, mincov), rep(0, half_covariances[1]-mincov), rep(1, half_covariances[2]-mincov), rep(1, half_covariances[3]-mincov))
        ramainder1 = max(df[1] - sum(part1), 0)
        ramainder2 = max(df[2] - sum(part2), 0)
        ramainder3 = max(df[3] - sum(part3), 0)
        Q1 = c(part1, rep(1, ramainder1), rep(0, ramainder2), rep(0, ramainder3))
        Q2 = c(part2, rep(0, ramainder1), rep(1, ramainder2), rep(0, ramainder3))
        Q3 = c(part3, rep(0, ramainder1), rep(0, ramainder2), rep(1, ramainder3))
        Qs = cbind(Q1, Q2, Q3)

        #chi-squared(1) sample
        bootstrap_sample = matrix((rnorm(length(Q1) * size))^2, nrow = size, ncol = length(Q1))
        chisq = bootstrap_sample %*% Qs
        approx_cov = 2* matrix(c(sum(Q1), half_covariances[1], half_covariances[2],
                                 half_covariances[1], sum(Q2), half_covariances[3],
                                 half_covariances[2], half_covariances[3], sum(Q3)),
                                nrow = 3, ncol = 3)
        return(list(sample = chisq, approx_cov = approx_cov))
}


