#' @title Multiple State Model Class
#' @exportClass msm
#' @examples udd <- new("msm", name="Uniform Distribution Deaths", 
#'     states = c("Alive", "Dead"),
#'     Qxt = Qxt.diag( matrix(c(
#'        uxt00 = function(t=0, x=0)
#'            0
#'        , uxt01 = function(t=0, x=0, w=120)
#'            1/(w-x-t)
#'        , uxt10 = function(t=0, x=0)
#'            0 * t
#'        , uxt11 = function(t=0, x=0)
#'            0
#'     ), 2, 2)))
#' print(udd)
#' @param representation Each multiple state model contains a name, state names and a transition matrix
setClass("msm", representation(name = "character",
                               states = "character",
                               Qxt = "matrix"))

#' @title Diagonal of transition matrix
#' @description Generates functions for the diagonal of the transition matrix, Q
#' @param Qxt.fns A matrix of transition functions
#' @param zero_diag Indicate whether the diagonal should be all functions which return 0
#' @param joint Indicate whether the tranisition matrix is for joint lives (i.e. functions of x and y)
#' @return A matrix of transition functions with updated diagonal entries
#' @export 
Qxt.diag <- function(Qxt.fns, zero_diag=FALSE, joint=FALSE) {
    N <- length(Qxt.fns[1,])
    tmp <- Qxt.fns[[1]] 
    
    for(i in seq(1,N*N,N+1)) {
        if(zero_diag) {
            if(joint) {
                Qxt.fns[[i]] <- function(t=0,x=0,y=0,muxt) {
                    0 * t
                }
            } else {
                Qxt.fns[[i]] <- function(t=0,x=0,muxt) {
                    0 * t
                }
            }
        } else {
            if(joint) {
                Qxt.fns[[i]] <- function(t=0,x=0,y=0,muxt) {
                    -sum(unlist(lapply(muxt,function(f) f(t,x,y))))
                }
            } else {
                Qxt.fns[[i]] <- function(t=0,x=0,muxt) {
                    -sum(unlist(lapply(muxt,function(f) f(t,x))))
                }
            }
        }
        
        if(N==1) { # i.e. a mis-specified single decrement model
            formals(Qxt.fns[[1]])$muxt <- list(tmp)
        } else {
            row <- ceiling(i/N)
            mu <- seq((row-1)*N + 1, row*N) ; mu <- mu[which(mu != i)]
            formals(Qxt.fns[[i]])$muxt <- Qxt.fns[mu]
        }
    }
    Qxt.fns
}

#' @title Transition Probabilities
#' @description Probability of transitioning from one state to another within t years. Approximation based on the Kolmogorov forward equations.
#' @param object An object of class msm
#' @param t The number of years forward
#' @param x The age of the individual.
#' @param h The time step. By default 1/12 is used for the approximation
#' @param all Indicates whether the probabilities at each time point should be returned or whether only the final probability matrix should be returned (at time t)
#' @param joint Indicate whether there are two lives, x and y that should should be used in the transition matrix functions
#' @param y The age of the second individual if the model is joint
#' @return A matrix of transition probabilities
#' @export 
tPx <- function(object, t=0, x=0, h=1/12, all=FALSE, joint=FALSE, y=0) {
    Qxt.fns <- object@Qxt
    N <- length(Qxt.fns[1,])
    
    if(t==0) {
        return(diag(N))
    } 
    
    tPx.value <- diag(N)
    it <- t/h
    tPx.all <- vector("list", it + 1)
    
    for(i in 1:it) {
        tPx.all[[i]] <- tPx.value
        
        # KFE
        if(joint) {
            Qxt.vals <- matrix(unlist(
                lapply(Qxt.fns, function(f) f(t=(i-1)*h,x,y))), 
                N, N, byrow=TRUE)
        } else {
            Qxt.vals <- matrix(unlist(
            lapply(Qxt.fns, function(f) f(t=(i-1)*h,x))), 
            N, N, byrow=TRUE)
        }
        tPx.value <- tPx.value + h * tPx.value %*% Qxt.vals
    }
    
    tPx.all[[it+1]] <- tPx.value
    
    if(all) 
        tPx.all
    else 
        tail(tPx.all,1)[[1]]
}

#' @title Show msm object
#' @description A "show" method for objects of class msm
#' @export
setMethod(f = "show", signature = "msm",
          definition = function(object ) {
              N <- length(object@states)
              x <- seq(0,N-1,1)
              states <- paste0("[", x, "] ", object@states)
              cat("Name     = ", object@name,
                  "\nStates   = ", states, 
                  "\nQx+t     = ")
              for(i in 1:N) {
                  qxt_msg <- ifelse(i==1, "", "           ")
                  for(j in 1:N) {
                      qxt_msg <- ifelse(j==i, paste0(qxt_msg, "-vx+t", i, "  "),
                                        paste0(qxt_msg, " ux+t", i-1, j-1, " "))
                  }
                  cat(qxt_msg, "\n")
              }
              
              for(i in 1:N) {
                  for(j in 1:N) {
                      qxt_msg <- ifelse(j==i, paste0("-vx+t", i, "  "),
                                        paste0(" ux+t", i-1, j-1, " "))
                      cat(qxt_msg)
                      if(j!=i) {
                          cat(" = ")
                          msg <- body(object@Qxt[[(i-1)*N+j]])
                          print(msg)
                      } else if (j==i) {
                          cat("\n")
                      }
                  }
              }
          })

#' @title Expected value of annuity (matrix form)
#' @description Calculates the expected values of annuities whose values depend on transitioning from one state to another
#' @param object An object of class msm 
#' @param n The number of years for the annuity
#' @param x The age of the individual
#' @param h The time step used in the approximation
#' @param delta The force of interest (continuously compounded interest rate)
#' @param dis Indicates that the annuity is discrete (with period h). By default, continuous annuities are calculated
#' @param joint Indicates whether the annuity depends on two lives
#' @param y The age of the second individual if the model is joint
#' @return A matrix of expected values for annuities 
#' @export 
ax <- function(object, n=10, x=0, h=1/12, delta=log(1.05), dis=0,
               joint=FALSE, y=0) {

    tPx.all <- tPx(object, t = n, x, h, all = TRUE, joint, y)
    
    it <- n/h
    for(i in 1:it) {
        t <- i*h
        tPx.all[[i+1]] <- tPx.all[[i+1]] * exp(-delta*t)
    }
    
    # simpsons rule
    a <- 0; b <- n
    h = (b-a)/it
    x1 = round(seq(a+h,b,h*2)/h + 1, 1)
    x2 = round(seq(a+2*h,b-h,h*2)/h + 1, 1)
    y1 = tPx.all[x1]
    y2 = tPx.all[x2]
    fa = tPx.all[[1]]
    fb = tPx.all[[it+1]]
    
    if(dis==1) {
        return(h*(Reduce("+", tPx.all)))
    }
    h/3 * (fa + fb + 4*Reduce("+", y1) + 2*Reduce("+", y2))
}

#' @title Expected value of insurance (matrix form)
#' @description Calculates the expected values of annuities whose values depend on transitioning from one state to another
#' @param object An object of class msm 
#' @param n The number of years for the insurance
#' @param x The age of the individual
#' @param h The time step used in the approximation
#' @param delta The force of interest (continuously compounded interest rate)
#' @param dis Indicates that the insurance is discrete (with period h). By default, continuous insurance are calculated
#' @param joint Indicates whether the insurance depends on two lives
#' @param y The age of the second individual if the model is joint
#' @return A matrix of expected values for insurances
#' @export 
Ax <- function(object, n=10, x=0, h=1/12, delta=log(1.05),
               joint=FALSE, y=0) {
    # we need a matrix of functions Qx+t' which is Qx+t with a zero diagonal
    tPx.all <- tPx(object, t = n, x, h, all = TRUE, joint, y)
    Qxt.prime <- Qxt.diag(object@Qxt, zero_diag = T)
    N <- length(Qxt.prime[1,])
    
    if(joint) {
        Qxt.vals <- matrix(unlist(
            lapply(Qxt.prime, function(f) f(t=0,x,y))), 
            N, N, byrow=TRUE)
    } else { 
        Qxt.vals <- matrix(unlist(
        lapply(Qxt.prime, function(f) f(t=0,x))), 
        N, N, byrow=TRUE)
    }
    tPx.all[[1]] <- tPx.all[[1]] %*% Qxt.vals
    
    it <- n/h
    for(i in 1:it) {
        t <- i*h
        
        # evaluate the Qxt.prime component
        if(joint) {
            Qxt.vals <- matrix(unlist(
            lapply(Qxt.prime, function(f) f(t=(i)*h,x,y))), 
            N, N, byrow=TRUE)
        } else { 
            Qxt.vals <- matrix(unlist(
            lapply(Qxt.prime, function(f) f(t=(i)*h,x))), 
            N, N, byrow=TRUE)
        }
        tPx.all[[i+1]] <- tPx.all[[i+1]] %*% Qxt.vals * exp(-delta*t)
    }
    
    # simpsons rule
    a <- 0; b <- n
    h = (b-a)/it
    x1 = round(seq(a+h,b,h*2)/h + 1, 1)
    x2 = round(seq(a+2*h,b-h,h*2)/h + 1, 1)
    y1 = tPx.all[x1]
    y2 = tPx.all[x2]
    fa = tPx.all[[1]]
    fb = tPx.all[[it+1]]
    h/3 * (fa + fb + 4*Reduce("+", y1) + 2*Reduce("+", y2))
}