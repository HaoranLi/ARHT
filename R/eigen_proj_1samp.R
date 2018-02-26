#' An internal function
#' @name eigen_proj_1samp
#' @keywords internal

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


