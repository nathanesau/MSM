#' @title Disability Model [msm]
#' @description A multiple state model with states Alive, Disabled and Dead
#' @export 
dm <- new("msm", name="Disability model",
          states = c("Alive", "Disabled", "Dead"),
          Qxt = Qxt.diag( matrix(c(
              uxt00 = function(t=0, x=0) 
                  0
              , uxt01 = function(t=0, x=0, a1=4e-04, b1=3.4674e-06, c1=0.138155) 
                  a1 + b1 * exp(c1*(x+t))
              , uxt02 = function(t=0, x=0, a2=5e-04, b2=7.5858e-05, c2=0.087498)
                  a2 + b2*exp(c2*(x+t))
              , uxt10 = function(t=0, x=0, a1=4e-04, b1=3.4674e-06, c1=0.138155)
                  0.1  * (a1 + b1  * exp(c1*(x+t)))
              , uxt11 = function(t=0, x=0) 
                  0
              , uxt12 = function(t=0, x=0, a2=5e-04, b2=7.5858e-05, c2=0.087498)
                  a2 + b2*exp(c2*(x+t))
              , uxt20 = function(t=0, x=0)
                  0 * t
              , uxt21 = function(t=0, x=0)
                  0 * t
              , uxt22 = function(t=0, x=0)
                  0
          ), 3, 3)))

#' @title Accidental Death Model [msm]
#' @description A multiple state model with states Alive, Disabled and Dead
#' @details This model is shown in Exercise 8.1 of Actuarial Mathematics for Life Contingent Risks (1st edition)
#' @export 
adm <- new("msm", name="Accidental Death Model",
           states=c("Alive", "Accidental", "Other"), 
           Qxt = Qxt.diag( matrix(c(
               uxt00 = function(t=0, x=0) 
                   0
               , uxt01 = function(t=0, x=0, A=5e-04, B=7.6e-05, c=1.09) 
                   A + B*c^(x+t)
               , uxt02 = function(t=0, x=0)
                    10^(-5)
               , uxt10 = function(t=0, x=0, a1=4e-04, b1=3.4674e-06, c1=0.138155)
                   0 * t
               , uxt11 = function(t=0, x=0) 
                   0
               , uxt12 = function(t=0, x=0, a2=5e-04, b2=7.5858e-05, c2=0.087498)
                   0 * t
               , uxt20 = function(t=0, x=0)
                   0 * t
               , uxt21 = function(t=0, x=0)
                   0 * t
               , uxt22 = function(t=0, x=0)
                   0
           ), 3, 3)))

#' @title Critical Illness Model [msm]
#' @description A multiple state model with states Healthy, Sick, Dead, Critical
#' @details This model is shown in Exercise 8.2 of Actuarial Mathematics for Life Contingent Risks (1st edition)
#' @export 
crit <- new("msm", name="Critical Illness Model",
           states=c("Healthy", "Sick", "Dead", "Critical"), 
           Qxt = Qxt.diag( matrix(c(
               uxt00 = function(t=0, x=0) 
                   0
               , uxt01 = function(t=0, x=0, a1=4e-04, b1=3.5e-06, c1=0.14)
                    a1 + b1 * exp(c1*(x+t))
               , uxt02 = function(t=0, x=0, a2=5e-04, b2=7.6e-05, c2=0.09)
                    a2 + b2 * exp(c2*(x+t))
               , uxt03 = function(t=0, x=0, a1=4e-04, b1=3.5e-06, c1=0.14)
                    0.05 * (a1 + b1 * exp(c1*(x+t)))
               , uxt10 = function(t=0, x=0, a1=4e-04, b1=3.5e-06, c1=0.14)
                    0.1 * (a1 + b1 * exp(c1*(x+t)))
               , uxt11 = function(t=0, x=0) 
                   0
               , uxt12 = function(t=0, x=0, a2=5e-04, b2=7.6e-05, c2=0.09)
                    a2 + b2 * exp(c2*(x+t))
               , uxt13 = function(t=0, x=0, a1=4e-04, b1=3.5e-06, c1=0.14)
                   0.05 * (a1 + b1 * exp(c1*(x+t)))
               , uxt20 = function(t=0, x=0)
                   0 * t
               , uxt21 = function(t=0, x=0)
                   0 * t
               , uxt22 = function(t=0, x=0)
                   0
               , uxt23 = function(t=0, x=0) 
                   0 * t
               , uxt30 = function(t=0, x=0)
                   0 * t
               , uxt31 = function(t=0, x=0)
                   0 * t
               , uxt32 = function(t=0, x=0, a2=5e-04, b2=7.6e-05, c2=0.09)
                   1.2 * (a2 + b2 * exp(c2*(x+t)))
               , uxt33 = function(t=0, x=0)
                   0 
           ), 4, 4)))

#' @title Common Shock Model [msm]
#' @description A multiple state model with states Both Alive, Husband Alive Only, Wife Alive Only, Both Dead
#' @details This model is shown in Exercise 8.5 of Actuarial Mathematics for Life Contingent Risks (1st edition)
#' @export 
common <- new("msm", name="Common Shock Model",
            states=c("Both", "Husband", "Wife", "Dead"), 
            Qxt = Qxt.diag( matrix(c(
                uxt00 = function(t=0, x=0, y=0) 
                    0
                , uxt01 = function(t=0, x=0, y=0, A=0.0001, B=0.0003, c=1.075)
                    A + B*c^(y+t)
                , uxt02 = function(t=0, x=0, y=0, A=0.0001, D=0.00035, c=1.075)
                    A + D*c^(x+t)
                , uxt03 = function(t=0, x=0, y=0)
                    5 * 10^(-5) * t^0
                , uxt10 = function(t=0, x=0, y=0)
                    0 * t
                , uxt11 = function(t=0, x=0, y=0) 
                    0
                , uxt12 = function(t=0, x=0, y=0)
                    0 * t
                , uxt13 = function(t=0, x=0, y=0)
                    0 * t
                , uxt20 = function(t=0, x=0, y=0)
                    0 * t
                , uxt21 = function(t=0, x=0, y=0)
                    0 * t
                , uxt22 = function(t=0, x=0, y=0)
                    0
                , uxt23 = function(t=0, x=0, y=0) 
                    0 * t
                , uxt30 = function(t=0, x=0, y=0)
                    0 * t
                , uxt31 = function(t=0, x=0, y=0)
                    0 * t
                , uxt32 = function(t=0, x=0, y=0)
                    0 * t
                , uxt33 = function(t=0, x=0, y=0)
                    0 
            ), 4, 4), joint=TRUE))

#' @title Constant Force of Mortality [msm]
#' @description An alive dead model with a constant force of mortality (i.e. exponential distribution).
#' @export 
cfm <- new("msm", name="Constant Force Mortality",
           states = c("Alive", "Dead"),
           Qxt = Qxt.diag( matrix(c(
               uxt00 = function(t=0, x=0) 
                   0
               , uxt01 = function(t=0, x=0, mu=0.05) 
                   mu
               , uxt10 = function(t=0, x=0)
                   0 * t
               , uxt11 = function(t=0, x=0) 
                   0
           ), 2, 2)))

#' @title Uniform Distribution of Deaths [msm]
#' @description An alive dead model based on De Moivre's Law (i.e. uniform distribution)
#' @export 
udd <- new("msm", name="Uniform Distribution Deaths",
           states = c("Alive", "Dead"),
           Qxt = Qxt.diag( matrix(c(
               uxt00 = function(t=0, x=0)
                   0
               , uxt01 = function(t=0, x=0, w=120)
                   1/(w-x-t)
               , uxt10 = function(t=0, x=0)
                   0 * t
               , uxt11 = function(t=0, x=0)
                   0
           ), 2, 2)))
