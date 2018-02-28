#' An adaptable generalized Hotelling's \eqn{T^2} test for high dimensional data
#' @export
#' @import stats
#' @description This function performs the adaptable regularized Hotelling's \eqn{T^2} test (ARHT) (Li et al., (2016) <arXiv:1609.08725>) for the one-sample
#'              and two-sample test problem, where we're interested in detecting the mean vector in the one-sample problem or the difference
#'              between mean vectors in the two-sample problem in a high dimensional regime.
#'
#' @details The method incorporates ridge-regularization in the classic Hotelling's \eqn{T^2} test with the regularization parameter
#'          chosen such that the asymptotic power under a class of probabilistic alternative prior models is maximized. ARHT combines
#'          different prior models by taking the maximum of statistics under all models. ARHT is distributed as the maximum
#'          of a correlated multivariate normal random vector. We estimate its covariance matrix and bootstrap its distribution. The
#'          returned p-value is a Monte Carlo approximation to its true value using the bootstrap sample, therefore not deterministic.
#'          Various methods are available to calibrate the slightly inflated Type 1 error rate of ARHT, including Cube-root transformation,
#'          square-root transformation and chi-square approximation.
#' @param X the n1-by-p observation matrix with numeric column variables.
#' @param Y an optional n2-by-p observation matrix; if \code{NULL}, a one-sample test is conducted on \code{X}; otherwise, a two-sample test
#'          is conducted on \code{X} and \code{Y}.
#' @param mu_0 the null hypothesis vector to be tested; if \code{NULL}, the default value is the 0 vector of length p.
#' @param prob_alt_prior a non-empty list; Each field is a numeric vector with sum 1. The default value is the "canonical weights"
#'        \code{list(c(1,0,0), c(0,1,0), c(0,0,1))}; Each field represents a probabilistic prior model specified by weights of \eqn{I_p},
#'        \eqn{\Sigma}, \eqn{\Sigma^2}, etc, where \eqn{\Sigma} is the population covariance matrix of the observations.
#' @param Type1error_calib the method to calibrate Type 1 error rate of ARHT. Choose its first element when more than one are specified.
#'        Four values are allowed:
#'        \itemize{\item{\code{cube_root}} The default value; cube-root transformation;
#'                 \item{\code{sqrt}} Square-root transformation;
#'                 \item{\code{chi_sq}} Chi-square approximation, not available when more than three models are specified in \code{prob_alt_prior};
#'                 \item{\code{none}} No calibration.
#' }
#' @param lambda_range optional user-supplied lambda range; If \code{NULL}, ARHT chooses its own range.
#' @param nlambda optional user-supplied number of lambda's in grid search; default to be \code{2000}; the grid is progressively coarser.
#' @param bs_size positive numeric with default value \code{1e5}; only effective when more than one prior models are specified in \code{prob_alt_prior};
#'        control the size of the bootstrap sample used to approximate the ARHT p-value.
#' @references Li, H. Aue, A., Paul, D. Peng, J., & Wang, P. (2016). \emph{An adaptable generalization of Hotelling's \eqn{T^2} test in high dimension.}
#'             <arXiv:1609:08725>.
#' @references Chen, L., Paul, D., Prentice, R., & Wang, P. (2011). \emph{A regularized Hotelling's \eqn{T^2} test for pathway analysis in proteomic studies.}
#'             Journal of the American Statistical Association, 106(496), 1345-1360.
#' @return \itemize{
#'  \item{\code{ARHT_pvalue}}: The p-value of ARHT test.
#'                             \itemize{
#'                              \item If \code{length(prob_alt_prior)==1}, it is identical to \code{RHT_pvalue}.
#'                              \item If \code{length(prob_alt_prior)>1}, it is the p-value after combining results from all prior models. The value is
#'                                    bootstrapped, therefore not deterministic.
#'                              }
#'  \item{\code{RHT_opt_lambda}}: The optimal lambda's chosen under each of the prior models in \code{prob_alt_prior}. It has the same length and order as
#'                              \code{prob_alt_prior}.
#'  \item{\code{RHT_pvalue}}: The p-value of RHT tests with the lambda's in \code{RHT_opt_lambda}.
#'  \item{\code{RHT_std}}: The standardized RHT statistics with the lambda's in \code{RHT_opt_lambda}.
#'  Take its maximum to get the statistic of ARHT test.
#'  \item{\code{Theta1}}: As defined in Li et al. (2016) <arXiv:1609.08725>, the estimated asymptotic means of RHT statistics with the lambda's in \code{RHT_opt_lambda}.
#'  \item{\code{Theta2}}: As defined in Li et al. (2016) <arXiv:1609.08725>, \code{2*Theta2} are the estimated asymptotic variances of RHT statistics the lambda's in \code{RHT_opt_lambda}.
#'  \item{\code{Corr_RHT}}: The estimated correlation matrix of the statistics in \code{RHT_std}.
#'}
#' @examples
#' set.seed(10086)
#' # One-sample test
#' n1 = 300; p =500
#' dataX = matrix(rnorm(n1 * p), nrow = n1, ncol = p)
#' res1 = ARHT(dataX)
#'
#' # Two-sample test
#' n2= 400
#' dataY = matrix(rnorm(n2 * p), nrow = n2, ncol = p )
#' res2 = ARHT(dataX, dataY, mu_0 = rep(0.01,p))
#'
#' # Specify probabilistic alternative priors model
#' res3 = ARHT(dataX, dataY, mu_0 = rep(0.01,p),
#'      prob_alt_prior = list(c(1/3, 1/3, 1/3), c(0,1,0)))
#'
#' # Change Type 1 error calibration method
#' res4 = ARHT(dataX, dataY, mu_0 = rep(0.01,p),
#'      Type1error_calib = "sqrt")
#'
#' RejectOrNot = res4$ARHT_pvalue < 0.05
#'

ARHT = function(X,
                Y = NULL,
                mu_0 = NULL,
                prob_alt_prior = list(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1)),
                Type1error_calib = c("cube_root", "sqrt", "chi_sq", "none"),
                lambda_range = NULL,
                nlambda = 2000,
                bs_size = 1e5){

        if(length(dim(X)) > 2L || !(is.numeric(X)))
                stop("X must be a numeric matrix with column variables")
        if(!is.matrix(X))
                X = as.matrix(X)
        if(nrow(X) <= 1L){
                stop("The number of rows in X must be larger than 1")
        }

        if(is.null(Y)){
                mode = "one_sample"
        }else{
                if(length(dim(Y)) > 2L || !(is.numeric(Y)))
                        stop("Y must be a numeric matrix with column variables")
                if(!is.matrix(Y))
                        Y <- as.matrix(Y)
                if(nrow(Y) <= 1L){
                        stop("The number of rows in Y must be larger than 1")
                }
                if(ncol(Y) != ncol(X)){
                        stop("The dimensions of X and Y differ")
                }
                mode = "two_sample"
        }

        if(!is.null(mu_0)){
                if(!is.vector(mu_0, mode = "numeric"))
                        stop("mu_0 must be a numeric vector")
                if(length(mu_0) != ncol(X))
                        stop("The dimension of X doesn't match with that of mu_0")
        }else{
                mu_0 = numeric(ncol(X))
        }

        if(!is.list(prob_alt_prior))
                stop("prob_alt_prior must be a list of numeric vectors")

        if(!all(sapply(prob_alt_prior, is.vector, mode = "numeric")))
                stop("prob_alt_prior must be a list of numeric vectors")

        valid_prob_alt_prior = sapply(prob_alt_prior, function(a){ round(sum(a), 5) != 1})
        if(any(valid_prob_alt_prior)){
                stop(paste("In Model", paste( which(valid_prob_alt_prior), collapse = ", "),
                           "specified in prob_alt_prior, the sum of prior weights is not 1"))
        }

        # throw away meaningless 0's and shorten prob_alt_prior
        max_nonzero_index = max(sapply(prob_alt_prior, function(xxx) max(which(xxx != 0))))
        prob_alt_prior = lapply(prob_alt_prior, function(xxx){
                xxx[1: min(length(xxx), max_nonzero_index)]})

        if(!is.null(lambda_range)){
                if(!is.vector(lambda_range, mode = "numeric"))
                        stop("lambda_range must be a numeric vector of two elements")
                if(length(lambda_range)!=2L)
                        stop("The length of lambda_range must be 2.")
                if(lambda_range[1]<=0)
                        stop("The lower bound of lambda sequence must be positive")
                if(lambda_range[2]<= lambda_range[1])
                        stop("The upper bound of lambda sequence must be larger than the lower bound")
        }

        if( (!is.numeric(nlambda)) || (length(nlambda)!= 1) )
                stop("nlambda must be numeric of length 1")
        if(nlambda<=0)
                stop("nlambda must be postive")
        nlambda = ceiling(nlambda)

        if(!(Type1error_calib[1] %in% c("cube_root", "sqrt", "chi_sq", "none"))){
                Type1error_calib = "cube_root"
                warning('Unknown value for Type1error_calib; default value "cube_root" is chosen instead')
        }
        if( (length(prob_alt_prior) >3L) && (Type1error_calib[1] == "chi_sq")){
                stop("Chi-square calibration of Type 1 error is not available when the number of prior models in prob_alt_prior
                     is larger than 3")
        }
        if(length(prob_alt_prior)>1L ){
                if( (!is.numeric(bs_size)) || (length(bs_size)!= 1) )
                        stop("bs_size must be numeric of length 1")
                if(bs_size <=0 )
                        stop("bs_size must be postive")
                if(bs_size < 1e3)
                        warning("Bootstrap sample size is too small; ARHT_pvalue is not reliable.")
        }
        bs_size = ceiling(bs_size)

        if((Type1error_calib[1] != "chi_sq") && (length(prob_alt_prior)>1L)){
                bootstrap_sample = matrix(rnorm(length(prob_alt_prior)*bs_size),
                                          ncol = bs_size)
        }

        if(mode == "one_sample"){
                eig_proj = eigen_proj_1samp(X, mu_0, lower_lambda = lambda_range[1])
        }else{
                eig_proj = eigen_proj_2samp(X, Y, mu_0, lower_lambda = lambda_range[1])
        }
        p = ncol(X)
        n = eig_proj$n
        gamma = p/n
        proj_diff = eig_proj$proj_shift
        ridge = eig_proj$lower_lambda

        # To speed up the computation of Stieltjes transform, separate positive eigenvalues and negative ones.
        positive_emp_eig = eig_proj$pos_eig_val
        num_zero_emp_eig = p - length(positive_emp_eig)
        emp_eig = c(positive_emp_eig, rep(0, times = num_zero_emp_eig))

        ## specify the lambda's net. Use log-scale. Progressively coarser.
        if(is.null(lambda_range)){
                lambda = exp(seq(from = log(ridge),
                                 to = log(20 * emp_eig[1] + (ridge - mean(emp_eig)/100) * (ridge - mean(emp_eig)/100 >0)),
                                 length = nlambda))
        }else{
                lambda = exp(seq(from = log(lambda_range[1]),
                                 to = log(lambda_range[2]),
                                 length = nlambda))
        }

        ## Stieltjes transform, its derivative, Theta_1, Theta_2
        mF = 1/p * ( rowSums(1/outer(lambda, positive_emp_eig, FUN = "+"))
                     + num_zero_emp_eig/lambda )

        mFprime = 1/p * (rowSums(1/(outer(lambda, positive_emp_eig, FUN = "+"))^2)
                         + num_zero_emp_eig/lambda^2)
        Theta1 = (1 - lambda*mF)/(1 - gamma*(1 - lambda * mF))

        Theta2 = (1 + gamma*Theta1)^2 * (Theta1 - lambda *(mF - lambda * mFprime)/(1 - gamma*(1 - lambda * mF))^2)

        # Calculate the power under each prior model
        prior_max_order = max(sapply(prob_alt_prior,length))
        unified_prob_alt_prior = lapply(prob_alt_prior, function(i) c(i, rep(0, times = max(prior_max_order,2) - length(i))))
        matrix_prob_alt_prior = do.call(rbind, unified_prob_alt_prior)
        if(prior_max_order <= 2L){
                rhos = rbind(mF, Theta1)
        }else{
                pop_moments = moments_PSD(emp_eig, n, prior_max_order-2)
                rhos = matrix(NA, nrow = prior_max_order, ncol = length(mF))
                rhos[1, ] = mF
                rhos[2, ] = Theta1
                # recursive formulas; cannot be parallel
                for(ii in 3:prior_max_order){
                        rhos[ii,] = (1 + gamma * Theta1) * (pop_moments[ii-2] - lambda * rhos[ii-1, ])
                }
        }
        powers = t(matrix_prob_alt_prior %*% rhos) / sqrt(2*gamma*Theta2) # Column: prior model; Row: lambda

        opt_lambda_index = apply(powers, 2, which.max) # optimal lambda index under each prior model

        ## Estimated covariance matrix of standardized RHT statistics with optimal lambda's
        G = matrix( apply( expand.grid(opt_lambda_index, opt_lambda_index), 1,
                           function(ddd){
                                aaa = ddd[1]
                                bbb = ddd[2]
                                if( abs(aaa - bbb) < 1e-8){
                                        return(1)
                                }else{
                                        return( (1 + gamma * Theta1[aaa]) * (1 + gamma * Theta1[bbb]) * (
                                                lambda[aaa] * Theta1[aaa] - lambda[bbb] * Theta1[bbb]) / (
                                                        (lambda[aaa] - lambda[bbb]) * sqrt(Theta2[aaa] * Theta2[bbb]))
                                                )
                                }
                                }),
                    nrow = length(opt_lambda_index), ncol = length(opt_lambda_index) )

        ## square root of G ##
        G_eigen = eigen(G,symmetric=T) ### project G to the ''closest'' nonnegative definite matrix
        G_evec = G_eigen$vectors
        G_eval = G_eigen$values
        G_eval_plus = G_eval * (G_eval >= 0)
        G_sqrt = G_evec %*% diag(sqrt(G_eval_plus))

        # standardized statistics
        RHT = sapply(lambda[opt_lambda_index], function(xx){
                (1/p) * sum( proj_diff^2 / (emp_eig + xx))}
        )
        if(Type1error_calib[1] != "chi_sq"){
                if(Type1error_calib[1] == "cube_root"){
                        RHT_std = {sqrt(p) * ( RHT^(1/3) - (Theta1[opt_lambda_index])^(1/3)) /
                                        sqrt(2*Theta2[opt_lambda_index]) / (1 / 3 * Theta1[opt_lambda_index]^(-2/3))}
                }
                if(Type1error_calib[1] == "sqrt"){
                        RHT_std = {sqrt(p) * (sqrt(RHT) - sqrt(Theta1[opt_lambda_index]))/
                                        sqrt(Theta2[opt_lambda_index] / 2 / Theta1[opt_lambda_index])}
                }
                if(Type1error_calib[1] == "none"){
                        RHT_std = (RHT - Theta1[opt_lambda_index]) / sqrt(2 * Theta2[opt_lambda_index] / p)
                }
                # p-values
                if(length(prob_alt_prior) == 1){
                        p_value = 1 - pnorm(RHT_std)
                        composite_p_value = p_value
                }else{
                        p_value = 1 - pnorm(RHT_std)
                        Tmax = apply(G_sqrt %*% bootstrap_sample,2,max)
                        composite_p_value = 1 - mean(max(RHT_std)>Tmax)
                }
        }

        if(Type1error_calib[1] == "chi_sq"){
                if(length(prob_alt_prior) == 1L){
                        # when one prior model is specified, no need for bootstrap
                        constant_coef = Theta2[opt_lambda_index] / Theta1[opt_lambda_index]
                        degree_freedom =  p * (Theta1[opt_lambda_index])^2 / Theta2[opt_lambda_index]
                        p_value = 1 - pchisq( p * RHT / constant_coef, df = degree_freedom)
                        composite_p_value = p_value
                }else{
                        if( length(prob_alt_prior) == 2L){
                                # Trick: add dummy variables to make the length of opt_lambda_index when less than 3 priors are specified
                                # max(RHT(lambda_1), RHT(lambda_2), RHT(lambda_1)) = max(RHT(lambda_1), RHT(lambda_2))
                                length3_opt_lambda_index = c(opt_lambda_index, opt_lambda_index[1])
                                # expand G to 3 variables
                                G_tmp = G_sqrt %*% t(G_sqrt)
                                G_expand = rbind( cbind(G_tmp, c(1, G_tmp[1,2])), c(1, G_tmp[1,2], 1))
                        }else{  # when 3 prior models are specified
                                length3_opt_lambda_index = opt_lambda_index
                                G_expand = G_sqrt %*% t(G_sqrt)
                        }
                        # Call r3chisq to get the generated 3-vairate chi-square bootstrap sample
                        constant_coef = Theta2[length3_opt_lambda_index] / Theta1[length3_opt_lambda_index]
                        degree_freedom = ceiling( p * (Theta1[length3_opt_lambda_index])^2 / Theta2[length3_opt_lambda_index])
                        chisq = r3chisq(size = bs_size, df = degree_freedom, corr_mat = G_expand)$sample
                        # standardize the bootstrap sample
                        T1sample = (1/p*constant_coef[1] * chisq[,1] - Theta1[opt_lambda_index[1]])/sqrt(2*Theta2[opt_lambda_index[1]]/p)
                        T2sample = (1/p*constant_coef[2] * chisq[,2] - Theta1[opt_lambda_index[2]])/sqrt(2*Theta2[opt_lambda_index[2]]/p)
                        T3sample = (1/p*constant_coef[3] * chisq[,3] - Theta1[opt_lambda_index[3]])/sqrt(2*Theta2[opt_lambda_index[3]]/p)
                        Tmax = pmax(T1sample, T2sample, T3sample)
                        p_value = 1 - pchisq( p * RHT / constant_coef[1:length(RHT)], df = degree_freedom[1:length(RHT)])
                        RHT_std = (RHT - Theta1[opt_lambda_index])/sqrt(2*Theta2[opt_lambda_index]/p)
                        composite_p_value = 1 - mean(max(RHT_std) > Tmax)
                }
        }
        return(list(ARHT_pvalue = composite_p_value,
                    RHT_pvalue = p_value,
                    RHT_std = RHT_std,
                    RHT_opt_lambda = lambda[opt_lambda_index],
                    Theta1 = Theta1[opt_lambda_index],
                    Theta2 = Theta2[opt_lambda_index],
                    Corr_RHT = G
                    ))
}


