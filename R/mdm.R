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

#' @title Non-diagonal of transition matrix
#' @description Sets functions for the non-diagonal of the transition matrix, Q, to 0
#' @param Qxt.fns A matrix of transition functions 
#' @param zero_diag Indicate whether the diagonal should be all functions which return 0
#' @param joint Indicate whether the tranisition matrix is for joint lives (i.e. functions of x and y)
#' @return A matrix of transition functions with updated non-diagonal entries
#' @details Should be called after Qxt.diag has been used on Qxt.fns. Used by tV() function
#' @export 
Qxt.nondiag <- function(Qxt.fns, joint=FALSE) {
    N <- length(Qxt.fns[1,])
    tmp <- Qxt.fns[[1]] 
    
    for(i in rep(seq(2:(N+1)), N-1) + rep(0:(N-2),each=N)*(N+1) + 1) {
        if(joint) {
            Qxt.fns[[i]] <- function(t=0,x=0,y=0,muxt) 0 * t
        } else {
            Qxt.fns[[i]] <- function(t=0,x=0,muxt) 0 * t
        }
    }
    Qxt.fns
}

#' @title Zero function matrix
#' @description Creates a matrix of functions. Each element is a function of t, which returns 0 * t
#' @param nrow The number of rows for the matrix
#' @param ncol The number of columns for the matrix
#' @return A matrix of functions (of t) which each return 0 
#' @export 
zfm <- function(nrow=1, ncol=1) {
    zfm.matrix <- list()
    for(i in (1:(nrow*ncol))) {
        zfm.matrix[[i]] <- function(t) 0 * t
    }
    matrix(zfm.matrix, nrow=nrow, ncol=ncol)
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
annx <- function(object, n=10, x=0, h=1/12, delta=log(1.05), dis=0,
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
    # we need a matrix of fuecho $DESKTOP_SESSIONnctions Qx+t' which is Qx+t with a zero diagonal
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

#' @title Reserve at time-t (matrix form)
#' @description Calculates the value of time-t Reserves by using Euler's method 
#'  with Thiele's differential equation
#' @param object An object of class msm 
#' @param t The time to calculate the reserve at
#' @param x The age of the policyholder
#' @param h The time step to use
#' @param V0 The initial reserve values (a forward recursion is used)
#' @param B A matrix of benefit payment functions to denote the rate of payment of net benefit at time t 
#'  while the policyholder is in state i. B should be a vector of length N where N is the number of states (i = 0... n)
#' @param S A matrix of lump sum benefit functions to denote the rate of lump sum benefit payments, payable at time t,
#'  on transition from state i to state j.
#' @param delta The force of interest (assumed constant)
#' @param all indicate whether all the tV values should be returned (by default only the last tV is returned)
#' @return The t-time reserve for each state
#' @export 
tV <- function(object, t=0, x=0, h=1/12, V0=NA, B=NA, S=NA, delta=log(1.05),
               all=FALSE) {
    Qxt.prime <- Qxt.nondiag(object@Qxt)
    N <- length(Qxt.prime[1,])

    # initial value of reserve
    if(is.na(V0[1])) tV.value <- matrix(0,nrow=N)
    else tV.value <- V0
    
    # benefit payment functions
    if(is.na(B[1])) B <- zfm(nrow=N, ncol=1)
    else B <- B

    # lump sum benefit functions
    S <- if(is.na(S[1])) S <- zfm(nrow=N,ncol=N)
    else S <- S
    
    it <- t/h
    tV.all <- vector("list", it + 1) # store all the reserves
    
    # i.e. if t = 1, h=1/12, then i = 12 is the number of iterations
    it <- t/h 
    for(i in 1:it) {
        tV.all[[i]] <- tV.value
        
        # evaluate transition intensities
        Qxt.vals <- matrix(unlist(
            lapply(Qxt.prime, function(f) f(t=(i-1)*h,x))), 
            N, N, byrow=TRUE)
        
        # evaluate benefit payment functions
        B.vals <- matrix(unlist(
            lapply(B, function(f) f(t=(i-1)*h))),
            N, 1)
        
        # evaluate lump sum benefit functions
        S.vals <- matrix(unlist(
            lapply(S, function(f) f(t=(i-1)*h))),
            N, N, byrow=TRUE)
    
        # V*
        V.vals <- matrix(rep(tV.value,N), N, byrow=TRUE)
        
        # tPx * = St + V* + t(V*)
        tPx.vals <- S.vals + V.vals + t(V.vals)
        
        tV.value <- tV.value + h * (delta * tV.value - B.vals - apply(Qxt.vals %*% tPx.vals, 1, sum))
    }
    
    tV.all[[it+1]] <- tV.value
    
    if(all) 
        tV.all
    else 
        tail(tV.all,1)[[1]]
}