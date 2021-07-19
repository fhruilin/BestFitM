#' @title  Fitting model
#' 
#' @description 
#' \code{FitM} calculate the parameters of the selected model.
#'  
#' @details 
#' This function calculate the parameters of the selected model. 
#' There are seven model parameters that can be calculated using this function.
#' The seven basic models include one simple linear fitting model(line2P),
#' one quadratic curve model(line3P),
#' one logarithmic model(log2P), 
#' two exponential models(exp2P and exp3P), 
#' and two power-law models(power2P and power3P).
#' 
#' @param data a dataframe containing the variables x and y.
#' @param x  name of the independent variable x.
#' @param y name of the dependent variable y.
#' @param model the best fitting model.
#' @param summary summarizing the model fits. Default is TRUE.
#' @param Pvalue.corrected if P-value corrected or not, the value is one of c("TRUE", "FALSE").
#' @export
#' @return 
#' @author Ruilin Huang <fhruilin@@163.com>
#' 
#' @examples 
#' data("mtcars")
#' FitM(data= mtcars, x= "mpg", y = "disp",model = "line2P")
FitM <- function(data,x= "x",y = "y", 
                  model = "line2P",
                  summary = TRUE,
                  Pvalue.corrected = TRUE){
  df <- data.frame(data)
  x <- df[,x]
  y <- df[,y]
  model <- model
  if(model == "line2P"){
    fitline2P <- lm(y ~ x) #line2P model
    if(summary == TRUE){
      summary(fitline2P)
    } else{
      print(fitline2P)
    }
  } else {
    if(model == "line3P") {
      fitline3P <- lm(y ~ I(x^2)+ x) # line3P model
      if(summary == TRUE) {
        summary(fitline3P)
      } else{
        print(fitline3P)
      }
    } else {
      if(model == "exp3P") {
        adjy_exp3P <- y-min(y)+1
        xadjy_exp3P <- data.frame(x,adjy_exp3P)
        lmFit_exp3P <- lm(log(adjy_exp3P) ~ x)
        coefs_exp3P <- coef(lmFit_exp3P)
        get.b_exp3P <- coefs_exp3P[2]   #slope
        nlsFit_exp3P <- nls(adjy_exp3P ~ cbind(1+exp(b*x),exp(b*x)),
                            start = list(b=get.b_exp3P),
                            data = xadjy_exp3P,algorithm = "plinear",
                            nls.control(maxiter = 5000000,minFactor = 10^(-10)))
        coef_exp3P <- coef(nlsFit_exp3P)
        b_exp3P <- coef_exp3P[1]
        c_exp3P <- coef_exp3P[2]+min(y)-1
        a_exp3P <- coef_exp3P[3]+coef_exp3P[2]
        fitexp3P <- nls(y ~ a*exp(b*x)+ c, 
                        start = list(a = a_exp3P, b = b_exp3P, c=c_exp3P))
        sfitexp3P <- summary(fitexp3P)
        #exp3P model
        #calculate p,r2 and adj.r2
        n=length(x)
        k = 3 # k means the count numbers of parameters(i.e., 'a', 'b', 'c'  in this case)
        ss.res<-sum((residuals(fitexp3P))^2) # Residual Sum of Squares, DF= n-k
        ss.total.uncor<-sum(y^2)             # Uncorrected Total Sum of Squares, DF=n
        ss.total.cor<-sum((y-mean(y))^2)     # Corrected Total Sum of Squares, DF=n-1
        if (Pvalue.corrected == TRUE){
          ss.reg <- ss.total.cor - ss.res  # Regression Sum of Squares, DF= (n-1)-(n-k) = k-1 in this case
          dfR= k-1
        }else{
          ss.reg <- ss.total.uncor - ss.res  # Regression Sum of Squares, DF= n-(n-k) = k in this case
          dfR= k
        }
        dfE= n-k
        F.value=(ss.reg/dfR)/(ss.res/dfE)
        p.value=pf(F.value,dfR,dfE,lower.tail = F)
        p.value<-unname(p.value)#p value
        
        RSE<- sfitexp3P$sigma  # Residual standard error, 
        SSE<- (RSE^2)*(n-1)      # Sum of Squares for Error, not equal to 'ss.res'.
        
        adjr2 <- 1-SSE/((var(y))*(n-1))
        r2<- 1-(1-adjr2)*((n-k)/(n-1))
        pAndr2 <- data.frame(F.value,p.value,r2,adjr2)
        if(summary== TRUE) {
          list(sfitexp3P,pAndr2)
        } else{
          print(fitexp3P)
        }
      } else {
        if(model == "log2P") {
          log2P <- function(x) {
            dt = vector()
            if(min(x) >0) {
              for (i in seq_along(x)) {
                dt[[i]] = log(x[[i]])
              }
            } else {
              stop("x values greater than 0 are required")
            }
            return(dt)
          }
          fitlog2P <- lm( y ~ log2P(x))#log2P model
          if(summary == TRUE) {
            summary(fitlog2P)
          } else{
            print(fitlog2P)
          }
        } else {
          if(model == "exp2P") {
            if(min(y) >0) {
              lmFit_exp2P <- lm(log(y) ~ x)
              coefs_exp2P <- coef(lmFit_exp2P)
              a_exp2P <- exp(coefs_exp2P[1])  #intercept
              b_exp2P <- coefs_exp2P[2]   #slope
              fitexp2P <- nls(y ~ a*exp(b*x), 
                              start = list(a = a_exp2P, b = b_exp2P))
              sfitexp2P <- summary(fitexp2P)
              #exp2P model
              #calculate p,r2 and adj.r2
              n=length(x)
              k = 2 # k means the count numbers of parameters(i.e., 'a', 'b' in this case)
              ss.res<-sum((residuals(fitexp2P))^2) # Residual Sum of Squares, DF= n-k
              ss.total.uncor<-sum(y^2)             # Uncorrected Total Sum of Squares, DF=n
              ss.total.cor<-sum((y-mean(y))^2)     # Corrected Total Sum of Squares, DF=n-1
              if (Pvalue.corrected == TRUE){
                ss.reg <- ss.total.cor - ss.res  # Regression Sum of Squares, DF= (n-1)-(n-k) = k-1 in this case
                dfR= k-1
              }else{
                ss.reg <- ss.total.uncor - ss.res  # Regression Sum of Squares, DF= n-(n-k) = k in this case
                dfR= k
              }
              dfE= n-k
              F.value=(ss.reg/dfR)/(ss.res/dfE)
              p.value=pf(F.value,dfR,dfE,lower.tail = F)
              p.value<-unname(p.value)#p value
              
              RSE<- sfitexp2P$sigma  # Residual standard error, 
              SSE<- (RSE^2)*(n-1)      # Sum of Squares for Error, not equal to 'ss.res'.
              
              adjr2 <- 1-SSE/((var(y))*(n-1))
              r2<- 1-(1-adjr2)*((n-k)/(n-1))
              pAndr2 <- data.frame(F.value,p.value,r2,adjr2)
              if(summary == TRUE) {
                list(sfitexp2P,pAndr2)
              } else{
                print(fitexp2P)
              }
            } else {
              stop("y values greater than 0 are required")
            }
          } else {
            if(model == "power2P") {
              if(min(x) > 0) {
                if(min(y) > 0) {
                  lmFit_power2P <- lm(log(y) ~ log(x))
                  coef_power2P <- coef(lmFit_power2P)
                  a_power2P <- exp(coef_power2P[1])#intercept
                  b_power2P <- coef_power2P[2]#slope
                  fitpower2P <- nls(y ~ a*(x^b),
                                    start = list(a = a_power2P, b = b_power2P)) 
                  sfitpower2P <- summary(fitpower2P)
                  #power2P model
                  #calculate p,r2 and adj.r2
                  n=length(x)
                  k = 2 # k means the count numbers of parameters(i.e., 'a', 'b'  in this case)
                  ss.res<-sum((residuals(fitpower2P))^2) # Residual Sum of Squares, DF= n-k
                  ss.total.uncor<-sum(y^2)             # Uncorrected Total Sum of Squares, DF=n
                  ss.total.cor<-sum((y-mean(y))^2)     # Corrected Total Sum of Squares, DF=n-1
                  if (Pvalue.corrected == TRUE){
                    ss.reg <- ss.total.cor - ss.res  # Regression Sum of Squares, DF= (n-1)-(n-k) = k-1 in this case
                    dfR= k-1
                  }else{
                    ss.reg <- ss.total.uncor - ss.res  # Regression Sum of Squares, DF= n-(n-k) = k in this case
                    dfR= k
                  }
                  dfE= n-k
                  F.value=(ss.reg/dfR)/(ss.res/dfE)
                  p.value=pf(F.value,dfR,dfE,lower.tail = F)
                  p.value<-unname(p.value)#p value
                  
                  RSE<-sfitpower2P$sigma  # Residual standard error, 
                  SSE<-(RSE^2)*(n-1)      # Sum of Squares for Error, not equal to 'ss.res'.
                  
                  adjr2 <- 1-SSE/((var(y))*(n-1))
                  r2<- 1-(1-adjr2)*((n-k)/(n-1))
                  pAndr2 <- data.frame(F.value,p.value,r2,adjr2)
                  
                  if(summary == TRUE) {
                    list(sfitpower2P,pAndr2)
                  } else{
                    print(fitpower2P)
                  }
                } else {
                  stop("y values greater than 0 are required")
                }
              } else {
                stop("x values greater than 0 are required")
              }
            } else {
              if(model == "power3P") {
                if(min(x) > 0) {
                  if(min(y) > 0) {
                    adjy_power3P = y-min(y)+1
                    xadjy_power3P <- data.frame(x,adjy_power3P)
                    lmFit_power3P <- lm(log(adjy_power3P) ~ log(x)) 
                    coefs_power3P <- coef(lmFit_power3P)
                    get.b_power3P <- coefs_power3P[2]   #slope
                    nlsFit_power3P  <- nls(adjy_power3P ~ cbind(1+x^b,x^b),
                                           start = list(b = get.b_power3P),
                                           data = xadjy_power3P,algorithm = "plinear",
                                           nls.control(maxiter = 5000000,minFactor = 10^(-10)))
                    coef_power3P <- coef(nlsFit_power3P)
                    b_power3P <- coef_power3P[1]
                    c_power3P <- coef_power3P[2]+min(y)-1
                    a_power3P <- coef_power3P[3]+coef_power3P[2]
                    fitpower3P <- nls(y ~ a*(x^b)+c, 
                                      start = list(a = a_power3P, b = b_power3P, c = c_power3P))
                    #power3P model
                    sfitpower3P <- summary(fitpower3P)
                    #calculate p,r2 and adj.r2
                    n=length(x)
                    k = 3 # k means the count numbers of parameters(i.e., 'a', 'b' and 'c' in this case)
                    ss.res<-sum((residuals(fitpower3P))^2) # Residual Sum of Squares, DF= n-k
                    ss.total.uncor<-sum(y^2)             # Uncorrected Total Sum of Squares, DF=n
                    ss.total.cor<-sum((y-mean(y))^2)     # Corrected Total Sum of Squares, DF=n-1
                    if (Pvalue.corrected == TRUE){
                      ss.reg <- ss.total.cor - ss.res  # Regression Sum of Squares, DF= (n-1)-(n-k) = k-1 in this case
                      dfR= k-1
                    }else{
                      ss.reg <- ss.total.uncor - ss.res  # Regression Sum of Squares, DF= n-(n-k) = k in this case
                      dfR= k
                    }
                    dfE= n-k
                    F.value=(ss.reg/dfR)/(ss.res/dfE)
                    p.value=pf(F.value,dfR,dfE,lower.tail = F)
                    p.value<-unname(p.value)#p value
                    
                    RSE<-sfitpower3P$sigma  # Residual standard error, 
                    SSE<-(RSE^2)*(n-1)      # Sum of Squares for Error, not equal to 'ss.res'.
                    
                    adjr2 <- 1-SSE/((var(y))*(n-1))
                    r2<-1-(1-adjr2)*((n-k)/(n-1))
                    pAndr2 <- data.frame(F.value,p.value,r2,adjr2)
                    
                    if(summary == TRUE){
                      list(sfitpower3P,pAndr2)
                    } else{
                      print(fitpower3P)
                    }
                  } else {
                    stop("y values greater than 0 are required")
                  }
                } else {
                  stop("x values greater than 0 are required") 
                } 
              } else{
                stop("Without this model")
              }
            }
          }
        }
      }
    }
  }
}



