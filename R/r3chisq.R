#' 3-variate positively correlated chi-squared sample generation when degrees of freedom are large
#'@description Generate samples approximately from three positively correlated chi-squared random variables
#' \eqn{(\chi^2(d_1), \chi^2(d_2), \chi^2(d_3))}
#' when the degrees of freedom \eqn{(d_1, d_2, d_3)} are large.
#' @name r3chisq
#' @details It is generally hard to sample from \eqn{(\chi^2(d_1), \chi^2(d_2), \chi^2(d_3))} with a designed correlation matrix.
#'  In the algorithm, we approximate the random vector by \eqn{(z^T Q_1 z, z^T Q_2 z, z^T Q_3 z)}
#'  where \eqn{z} is a standard norm random vector and \eqn{Q_1,Q_2,Q_3} are diagonal matrices
#'  with diagonal elements 1's and 0's. The designed positive correlations is approximated by carefully
#'  selecting common locations of 1's on the diagonals. The generated sample may have slightly larger marginal degrees
#'  of freedom than the inputted \code{df}, also slightly different covariances.
#' @param size sample size.
#' @param df the degree of freedoms of the marginal distributions. Must be non-negative, but can be non-integer.
#'  The function uses \code{ceiling(df)} if non-integer.
#' @param corr_mat the target correlation matrix; negative elements will be set to 0.
#' @return
#' \itemize{
#' \item{\code{sample}}: a \code{size}-by-3 matrix contains the generated sample.
#' \item{\code{approx_cov}}: the true covariance matrix of \code{sample}.}
#' @references Li, H., Aue, A., Paul, D., Peng, J., & Wang, P. (2016). \emph{
#' An adaptable generalization of Hotelling's  \eqn{T^2}
#'  test in high dimension.} arXiv preprint <arXiv:1609.08725>.
#' @examples
#' set.seed(10086)
#' cor_examp = matrix(c(1,1/6,2/3,1/6,1,2/3,2/3,2/3,1),3,3)
#' a_sam = r3chisq(size = 10000,
#'                 df = c(80,90,100),
#'                 corr_mat = cor_examp)
#' cov(a_sam$sample) - a_sam$approx_cov
#' cov2cor(a_sam$approx_cov) - cor_examp
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


