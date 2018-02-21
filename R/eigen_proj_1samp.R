#' An internal function
#' @name eigen_proj_1samp
#' @keywords internal
#' @description This is an internal function for \code{ARHT}, computing summary statistics for the one- sample mean test problem.
#'  The function returns positive eigenvalues of the sample covariance matrix,
#'  and the scaled projection onto the eigenspace of the sample covariance matrix of the distance between the sample
#'   mean \eqn{\bar{X}} and \code{mu_0}, that is,
#'   \deqn{n^{1/2} P^T (\bar{X} - mu_0),}
#'  where \eqn{P} is the eigenvector matrix of the sample covariance matrix.
#'  The function is designed to handle the situation when
#'  singular value algorithm in \code{svd()} does not converge, by adding ridge regularization term to
#'  the sample covariance matrix.
#' @param X \code{X} as in \code{\link{ARHT}}.
#' @param mu_0 \code{mu_0} as in \code{\link{ARHT}}.
#' @param lower_lambda the lower bound of the lambda sequence, that is \code{lambda_range[1]} as in \code{\link{ARHT}} if specified;
#'  otherwise \code{NULL}.
#' @seealso \code{\link{ARHT}}, \code{\link{eigen_proj_2samp}}
#' @return
#' \itemize{\item{\code{n}}: degree of freedom of the sample covariance matrix.
#'          \item{\code{pos_eig_val}}: positive eigenvalues of the sample covariance matrix.
#'          \item{\code{proj_shift}}: the scaled projection of the distance between the sample mean and \code{mu_0}; see Description.
#'          \item{\code{lower_bound}}: If \code{lower_lambda} specified, returns its value; if not,
#'          returns the generated lower bound of the lambda sequence.
#'}
#'@examples
#' set.seed(10086)
#' # One-sample test
#' n1 = 300; p = 500
#' X = matrix(rnorm(n1 *p), nrow = n1, ncol = p )
#' eigen_proj_1samp(X, rep(0.01, times = p))
#' @export
eigen_proj_1samp = function(X,
                            mu_0,
                            lower_lambda = NULL){
        N = nrow(X)
        p = ncol(X)
        nn = N-1L
        #gamma = p/nn
        X_bar = colMeans(X)
        half_S = 1 / sqrt(nn) * t(X) %*%  (diag(1,nrow = N, ncol = N) - 1/N * matrix(1,N,N)) # S = half_S %*% t(half_S)

        svd_half_S = try(svd(half_S, nu = p ,  nv = 0), silent = TRUE)

        # Handle the situation where svd fails to converge
        if(inherits(svd_half_S,"try-error")){
                S = half_S %*% t(half_S)
                # If lower_lambda specified, add the lower bound to S, then svd.
                # If not specified, generate recommended lower bound.
                # Initially (mean(diag(S))/100), if not enough, ridge = 1.5 * ridge.
                if(!is.null(lower_lambda)){
                        ridge = lower_lambda
                        svdofS_ridge = try(svd(S + diag(ridge, nrow = p), nu = p , nv = 0), silent = TRUE)
                        if(inherits(svdofS_ridge, "try-error")){
                                stop("The lower bound of lambda sequence is too small.")
                        }
                }else{  # the generated lower bound of lambda sequence
                        ridge = (mean(diag(S))/100)
                        svdofS_ridge = try(svd(S + diag(ridge, nrow = p), nu = p , nv = 0), silent = TRUE)
                        loop_counter = 0L
                        while(inherits(svdofS_ridge, "try-error") && loop_counter<=20L){
                                ridge = ridge * 1.5
                                loop_counter = loop_counter +1L
                                svdofS_ridge = try(svd(S + diag(ridge, nrow = p), nu = p , nv = 0), silent = TRUE)
                        }
                        if(loop_counter > 20L)
                                stop("singular value algorithm in svd() did not converge")
                }
                # projection of X_bar - mu_0 to the eigenspace of S.
                emp_evec = svdofS_ridge$u
                # To speed up the computation of Stieltjes transform, separate positive eigenvalues and negative ones.
                emp_eig  = (svdofS_ridge$d - ridge) * (svdofS_ridge$d >= ridge)
                positive_emp_eig = emp_eig[emp_eig > (1e-8) * mean(emp_eig)]
                #num_zero_emp_eig = p - length(positive_emp_eig) # number of 0 eigenvalues
        }else{
                emp_evec = svd_half_S$u
                eig_raw = svd_half_S$d^2
                positive_emp_eig = eig_raw[eig_raw > (1e-8 ) * mean(eig_raw)]
                num_zero_emp_eig = p - length(positive_emp_eig)
                emp_eig = c(positive_emp_eig, rep(0, times = num_zero_emp_eig))
                if(!is.null(lower_lambda)){
                        ridge = lower_lambda
                }else{
                        ridge =  (mean(emp_eig)/100)
                }
        }
        proj = as.vector(sqrt(N) * t(emp_evec) %*% (X_bar - mu_0))
        return(list(n = nn, pos_eig_val = positive_emp_eig, proj_shift = proj, lower_lambda = ridge))
}


