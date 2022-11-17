# Statistical Programming Practical 4
# Group Member: Enyan WU (s2303128),
#               Huiyu WU (s2303136), 
#               Xuantao LI (s1822046) 

# github repo address: https://github.com/wenyan512/SP_Practical4_Group37.git

################################ Group contribution ############################

# Overall, the weights of contribution for each group member are fairly equal. 
# We wrote our own coding based on our personal understanding of these questions. 
# We shared our ideas when anyone stick in somewhere and all of us enjoyed the teamwork. 
# The final version upload now is the final work we all committed for.


############################### Project introduction ###########################

# This project aims to implement Newton’s method for minimization of functions by
# creating an optimization function that can operate broadly.

# Here is a brief description of Newton’s method:
# We start from an initial guess at the parameters, and find the value of the 
# objective function (D), the gradient vector and the Hessian matrix at 
# that guess. We then find an improved guess by adding the delta (a descent 
# direction that will the objective function) to the current guess, and repeat
# the process, until we reach at the convergence (which is judged by a certain
# criteria specified below).

# To avoid any possibility of divergence, we have to check that each new proposal 
# for the parameters has actually reduced D and if not, repeatedly halve delta  
# until it does. We also have to ensure that the Hessian matrix is positive definite,
# (tested by whether it can be decomposed by Cholesky), and perturbing it to be so, 
# if it isn’t.

# Cholesky decomposition:
# With positive definite matrices, tasks can be tackled especially easily and 
# efficiently, which is largely because a positive definite matrix can be decomposed 
# into the Cholesky factor (the crossproduct of an upper triangular matrix with itself).
# Thus, in this project, instead of using solve(), we use chol(), since the constant 
# of proportionality is smaller for Cholesky.

# Before generating the newt(), we created two additional function.
# 1) a function that approximates the Hessian matrix when it is not provided
# 2) a function that checks whether the Hessian matrix is positive definite or 
#    not at convergence 


################################## Functions ###################################

# Function purpose:
# This function aims to approximate the Hessian matrix when it is not provided by 
# finite differencing the gradient vector 

# Input:
# theta -> a vector of the values for the optimization parameters 
# grad -> gradient function of the objective
# eps -> the finite difference intervals

# Output:
# an approximation to the Hessian by finite differencing of the gradient vector

approx_H <- function(theta, grad, eps, ...) {
  
  # initialize Hessian matrix
  hess <- matrix(0, nrow = length(theta), ncol = length(theta))
  
  for (i in 1:length(theta)){
    # theta vector after adding the finite difference intervals (eps)
    theta_eps <- theta
    theta_eps[i] <- theta[i] + eps
    
    # generate the values of Hessian matrix for each column
    hess[,i] <- (grad(theta_eps,...) - grad(theta,...))/eps
  }
  # create a symmetric version of the generated Hessian matrix since a matrix must  
  # be first symmetric, then it can be positive definite under certain conditions
  hess <- 0.5 * (t(hess) + hess)
  
  return(hess)
}


# Function purpose:
# This function aims to check whether the Hessian matrix is positive definite at 
# convergence

# Input:
# H -> Hessian matrix which need to be judged 

# Output:
# When the Hessian matrix is positive definite, return the inverse of the Hessian 
# matrix, otherwise, stop the function and issue the error information.

check_positive <- function(H){
  # Cholesky decomposition fails for matrices that are not positive definite
  check <- try(H_chol <- chol(H), silent = TRUE)
  if ('try-error' %in% class(check)){
    stop('The Hessian is not positive definite at convergence.')
  } else{
    # calculate the inverse of the Hessian matrix
    Diag <- diag(dim(H_chol)[1])
    H_inverse <- backsolve(H_chol,forwardsolve(t(H_chol), Diag))
  }
  return(H_inverse)

}


# Function purpose:
# This function aims to minimize the objective function by Newton’s method

# Input:
# theta -> a vector of initial values for the optimization parameters
# func -> the objective function need to be minimized
# grad -> gradient function of the objective
# hess -> the Hessian matrix function of the objective
# ... -> any arguments of func, grad and hess after the 'theta' are passed using this
# tol -> the convergence tolerance
# fscale -> a rough estimate of the magnitude of func near the optimum
# maxit -> maximum number of Newton iterations to try before giving up
# max.half -> maximum number of times a step should be halved before concluding 
#            that the step has failed to improve the objective.
# eps -> the finite difference intervals to use when a Hessian function is not provided

# Output:
# if the function can find the minimization within the maximum number of Newton 
# iterations successfully, it will return a list contains the following:

# 1) f -> the value of the objective function at the minimum
# 2) theta -> the value of the parameters at the minimum
# 3) iter -> the number of iterations taken to reach the minimum
# 4) g -> the gradient vector at the minimum
# 5) Hi -> the inverse of the Hessian matrix at the minimum 

# Otherwise, the function will stop with the error or warnings under the following 
# circumstances:

# 1) If the objective or derivatives are not finite at the initial *theta*
# 2) If the step fails to reduce the objective despite trying *max.half* step halvings
# 3) If *maxit* is reached without convergence
# 4) If the Hessian is not positive definite at convergence

newt <- function(theta,func,grad,hess=NULL,...,tol=1e-8,fscale=1,maxit=100,
                 max.half=20,eps=1e-6,k=2){
  
  # obtain the value of the objective function at the initial values of the parameters
  ob_func <- func(theta,...)
  # obtain the gradient vector at initial point
  gradient <- grad(theta,...) 
  # obtain the Hessian matrix at initial point
  
  # if the Hessian matrix is not supplied, then calculate the approximation to it 
  # by finite differencing of the gradient vector
  if (is.null(hess)) {
    # the approximation is generated from the function approx_H()
    H <- approx_H(theta, grad, eps,...)
  } 
  else{
    H <- hess(theta,...) 
  }
  
  # check if the value of the objective function are not finite at the initial theta
  if (is.infinite(ob_func)){
    # if infinite, stop the function and issue an error
    stop("The value of the objective function is not finite at the initial theta.")
  }
  # check if the derivatives are not finite at the initial theta
  if (any(is.infinite(gradient)) | any(is.infinite(H))){
    # if infinite, stop the function and issue an error
    stop("The derivatives are not finite at the initial theta.")
  }
  
  
##--------------------------- Start Minimization -----------------------------##
  
  # keep repeating the process until we arrive at the optimal convergence
  # initialize the number of Newton iterations 
  iter = 0
  # loop until the maximum number of Newton iterations to try before giving up
  while (iter < maxit) {
    
    # convergence reached before getting to the maximum number of Newton iterations
    # convergence is judged by seeing whether all elements of the gradient vector 
    # have absolute value less than *tol* times the absolute value of the 
    # objective function plus *fscale*
    if (max(abs(gradient)) < (abs(ob_func)+fscale)*tol){
      
      # check if the Hessian matrix is positive definite at convergence
      # if Hessian is positive definite, return the following list,
      # Otherwise this function (newt) will be stopped
      H_inverse <- check_positive(H)
      return(list('f'=ob_func, 'theta'=theta, 'iter'=iter, 
                  'g'=gradient, 'Hi'=H_inverse))
    }
    
    else{
      
      # if the Hessian matrix is not positive definite, perturb it to be so.
      # obtain the result of try
      # (tested by Cholesky decomposition -> only positive definite matrix can pass)
      result <- try({chol(H)}, silent = TRUE)
      while ("try-error" %in% class(result)) {
        
        # add an identity matrix to the original Hessian
        # and for each iteration, add one more identity matrix
        H <- H + diag(dim(H)[1])
        # obtain the result of try
        # (tested by Cholesky -> only positive definite matrix can pass)
        result <- try({chol(H)}, silent = TRUE)
        
      }
      # after obtaining the positive definite Hessian matrix, decomposed it by 
      # Cholesky for later calculation
      H_chol <- chol(H)
      
      # calculate delta -> a descent direction that will reduce the objective function
      delta <- backsolve(H_chol, forwardsolve(t(H_chol), -gradient))
      
      # delta might overshoot and increase the objective function
      # If so, repeatedly halve delta until the value of the objective is decreased
      
      # initialize a counter for the number of times a step halved 
      half_iter <- 0
      # only halve delta when the value of the objective at the current theta plus 
      # delta is greater than the current value of the objective 
      while (func(theta + delta, ...) > ob_func) {
        # check if the step fails to reduce the objective despite trying *max.half*
        # step halvings
        if (half_iter < max.half) {
          delta <- delta / 2
          half_iter <- half_iter + 1
        } 
        else {
          # if the step fails to reduce the objective despite trying *max.half*
          # step halvings
          stop(paste("The step fails to reduce the objective despite trying", 
                     as.character(max.half), "halvings"))
        }
      }
      
      # update theta
      theta <- theta + delta
      # update function values
      ob_func <- func(theta, ...)
      # update gradient vector
      gradient <- grad(theta,...) 
      # update Hessian matrix
      # if the Hessian matrix is not supplied, then calculate the approximation to it 
      # by finite differencing of the gradient vector
      if (is.null(hess)) {
        # the approximation is generated from the function approx_H()
        H <- approx_H(theta, grad, eps,...)
      } 
      else{
        H <- hess(theta,...) 
      }
      
      # check if the value of the objective function are not finite
      if (is.infinite(ob_func)){
        # if infinite, stop the function and issue an error
        stop(paste("The value of the objective function is not finite at",
                   as.character(iter), "iteration(s)"))
      }
      # check if the derivatives are not finite
      if (any(is.infinite(gradient)) | any(is.infinite(H))){
        # if infinite, stop the function and issue an error
        stop(paste("The derivatives are not finite at",as.character(iter), 
                   "iteration(s)"))
      }
      
      # Iteration update
      iter <- iter + 1
      
    }
  }
  
  # check whether convergence occurs when reaching the maximum number of iterations
  if (max(abs(gradient)) < (abs(ob_func)+fscale)*tol){
    # check if the Hessian matrix is positive definite at convergence
    H_inverse <- check_positive(H)
    return(list('f'=ob_func, 'theta'=theta, 'iter'=iter, 
                'g'=gradient, 'Hi'=H_inverse))
    
  } 
  
  # if the maximum number of iterations is reached without convergence, issue the warnings
  else {
    warning(paste("Newton optimizer failed to converage after",
                  as.character(maxit), " iterations"))
  }
}


