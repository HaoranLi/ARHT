#' Adaptable Regularized Hotelling's T^2 test for one-sample testing problem
#' @param X n-by-p observation matrix with  numeric column variables.
#' @param mu_0 the null hypothesis to be tested, default value is the 0 vector.
#' @param prob_alternative_prior an umempty list with each field a numeric vector of length 3 with sum 1.
#'                          Default value is list(c(1,0,0), c(0,1,0), c(0,0,1)).
#'                          Each field of the list represents a probabilistic prior models specified by weights
#'                          of $\Sigma^0$, $\Sigma^1$, $\Sigma^2$, etc.
#' @param Type1error_calib the method to calibrate Type 1 error rate of the test. Four values are allowed,
#'                    cube_root (default, cube-root transformation),
#'                    sqrt (square-root transformation),
#'                    chi_sq (chi-square approximation),
#'                    none (no calibration).
#' @param lambda_range Optional user-supplied lambda sequence; default is NULL, and ARHT chooses
#'    its own sequence
#' @param nlambda  Optional user-supplied number of lambda's in grid search; Default is 2000; Must be postive numeric. The grid is progressively coarser.
#'
#' @param bs_size a numeric scalar, default value 1e6, only available when more than one prior models are
#'           specified in prob_alternative_prior; control the size of bootstrap sample used to approximate p-values.
#' @return Tstat: a vector of the same length of prob_alternative_prior.

ARHT = function(X,
                Y = NULL,
                mu_0 = NULL,
                prob_alternative_prior = list(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1)),
                Type1error_calib = c("cube_root", "sqrt", "chi_sq", "none"),
                lambda_range = NULL,
                nlambda = 2000,
                #bootstrap_sample = NULL,
                bs_size = 1e6){

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

        if(!is.list(prob_alternative_prior))
                stop("prob_alternative_prior must be a list of numeric vectors")

        if(!all(sapply(prob_alternative_prior, is.vector, mode = "numeric")))
                stop("prob_alternative_prior must be a list of numeric vectors")

        if(any(sapply(prob_alternative_prior, function(a){sum(a) != 1})))
                stop("prior weights must have sum 1")

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

        if(!(Type1error_calib[1] %in% c("cube_root", "sqrt", "chi_sq", "original"))){
                Type1error_calib = "cube_root"
                warning('Unknown value for Type1error_calib; default value "cube_root" is chosen instead')
        }
        if( (length(prob_alternative_prior) >3L) && (Type1error_calib[1] == "chi_sq")){
                stop("Chi-square calibration of Type 1 error is not available when the number of prior models
                     is larger than 3")
        }
        if(length(prob_alternative_prior)>1L ){
                if( (!is.numeric(bs_size)) || (length(bs_size)!= 1) )
                        stop("bs_size must be numeric of length 1")
                if(bs_size <=0 )
                        stop("bs_size must be postive")
                if(bs_size < 1e3)
                        warning("Bootstrap sample size is too small; estimated p-value is not reliable.")
        }
        bs_size = ceiling(bs_size)

        if((Type1error_calib[1] != "chi_sq") && (length(prob_alternative_prior)>1L)){
                bootstrap_sample = matrix(rnorm(length(prob_alternative_prior)*bs_size),
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
        prior_max_order = max(sapply(prob_alternative_prior,length))
        unified_prob_alternative_prior = lapply(prob_alternative_prior, function(i) c(i, rep(0, times = max(prior_max_order,2) - length(i))))
        matrix_prob_alternative_prior = do.call(rbind, unified_prob_alternative_prior)
        if(prior_max_order <= 2L){
                rhos = rbind(mF, Theta1)
        }else{
                pop_moments = moments_PSD(emp_eig, n, prior_max_order-2)
                rhos = matrix(NA, nrow = prior_max_order, ncol = length(mF))
                rhos[1,] = mF
                rhos[2,] = Theta1
                # recursive formulas; cannot be parallel
                for(ii in 3:prior_max_order){
                        rhos[ii,] = (1 + gamma * Theta1) * (pop_moments[ii-2] - lambda * rhos[ii-1])
                }
        }
        powers = t(matrix_prob_alternative_prior %*% rhos) / sqrt(2*gamma*Theta2) # Column: prior model; Row: lambda
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
                if(length(prob_alternative_prior) == 1){
                        p_value = 1 - pnorm(RHT_std)
                        composite_p_value = p_value
                }else{
                        p_value = 1 - pnorm(RHT_std)
                        Tmax = apply(G_sqrt %*% bootstrap_sample,2,max)
                        composite_p_value = 1 - mean(max(RHT_std)>Tmax)
                }
        }

        if(Type1error_calib[1] == "chi_sq"){
                if(length(prob_alternative_prior) == 1L){
                        # when one prior model is specified, no need for bootstrap
                        constant_coef = Theta2[opt_lambda_index] / Theta1[opt_lambda_index]
                        degree_freedom =  p * (Theta1[opt_lambda_index])^2 / Theta2[opt_lambda_index]
                        p_value = 1 - pchisq( p * RHT / constant_coef, df = degree_freedom)
                        composite_p_value = p_value
                }else{
                        if( length(prob_alternative_prior) == 2L){
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
                        G_expand = G_expand * (G_expand >=0)

                        #constant * chisq(degree_freedom); round covariances down
                        constant_coef = Theta2[length3_opt_lambda_index] / Theta1[length3_opt_lambda_index]
                        degree_freedom = ceiling( p * (Theta1[length3_opt_lambda_index])^2 / Theta2[length3_opt_lambda_index])
                        covariance3 = floor(c(  (degree_freedom[1])^(1/2) * (degree_freedom[2])^(1/2) * G_expand[1,2],
                                                (degree_freedom[1])^(1/2) * (degree_freedom[3])^(1/2) * G_expand[1,3],
                                                (degree_freedom[2])^(1/2) * (degree_freedom[3])^(1/2) * G_expand[2,3]))

                        # Build three diag matrices Q1, Q2, Q3.  z^TQ_1z, z^TQ_2z, z^TQ_3z would be the chi-square vector
                        # See supplementary material for algorithm explanation
                        mincov = min(covariance3)
                        part1 = c(rep(1,mincov), rep(1,covariance3[1]-mincov), rep(1,covariance3[2]-mincov), rep(0,covariance3[3]-mincov))
                        part2 = c(rep(1,mincov), rep(1,covariance3[1]-mincov), rep(0,covariance3[2]-mincov), rep(1,covariance3[3]-mincov))
                        part3 = c(rep(1,mincov), rep(0,covariance3[1]-mincov), rep(1,covariance3[2]-mincov), rep(1,covariance3[3]-mincov))
                        ramainder1 = max(degree_freedom[1] - sum(part1), 0)
                        ramainder2 = max(degree_freedom[2] - sum(part2), 0)
                        ramainder3 = max(degree_freedom[3] - sum(part3), 0)
                        Q1 = c(part1, rep(1, ramainder1), rep(0, ramainder2), rep(0, ramainder3))
                        Q2 = c(part2, rep(0, ramainder1), rep(1, ramainder2), rep(0, ramainder3))
                        Q3 = c(part3, rep(0, ramainder1), rep(0, ramainder2), rep(1, ramainder3))

                        #chi-squared(1) sample
                        bootstrap_sample = matrix((rnorm(length(Q1) * bs_size))^2, nrow = length(Q1), ncol = bs_size)
                        chisq1 = colSums(Q1 * bootstrap_sample)
                        chisq2 = colSums(Q2 * bootstrap_sample)
                        chisq3 = colSums(Q3 * bootstrap_sample)
                        T1sample = (1/p*constant.coef[1] * chisq1 - Theta1[1])/sqrt(2*Theta2[1]/p)
                        T2sample = (1/p*constant.coef[2] * chisq2 - Theta1[2])/sqrt(2*Theta2[2]/p)
                        T3sample = (1/p*constant.coef[3] * chisq3 - Theta1[3])/sqrt(2*Theta2[3]/p)
                        Tmax = pmax(T1sample, T2sample, T3sample)

                        RHT_std = (RHT - Theta1[opt_lambda_index])/sqrt(2*Theta2[opt_lambda_index]/p)
                        p_value = 1 - mean(max(RHT_std)>Tmax)
                }
        }
}

