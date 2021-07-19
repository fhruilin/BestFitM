#' @title Best fitting model
#' 
#' @description 
#' \code{bestFitM} selected the best fitting model from seven 
#' basic fitting models, according to AIC and BIC values.
#' 
#' @details 
#' This function selected the best fitting model from seven 
#' basic fitting models. The seven basic models include 
#' one simple linear fitting model(line2P), one quadratic curve model(line3P),
#' one logarithmic model(log2P), two exponential models(exp2P and exp3P),
#' and two power-law models(power2P and power3P).
#' 
#' @param data a dataframe containing the variables x and y.
#' @param x  name of the independent variable x.
#' @param y name of the dependent variable y.
#' @export
#' @return a list with 2 elements:
#' \item{BestFitM}{the best fitting model}
#' \item{AllModel}{all seven basic models, NA represents the model is not applicable}
#' @author Ruilin Huang <fhruilin@@163.com>
#' @examples 
#' data("mtcars")
#' bestFitM(mtcars, x= "mpg", y = "disp")
bestFitM <- function(data,x= "x",y = "y") {
  dat <- data.frame(data)
  x <- dat[,x]
  y <- dat[,y]
  if(is.na(x) || is.na(y))
  {
    stop("There may be missing values in the variable x or y")
  } else{
    fitline2P <- lm(y ~ x)#line2P model
    fitline3P <- lm(y ~ I(x^2)+ x)#line3P model
    
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
    #exp3P model
    AIC <- c(AIC(fitline2P),AIC(fitline3P),AIC(fitexp3P))
    BIC <- c(BIC(fitline2P),BIC(fitline3P),BIC(fitexp3P))
    if(min(x) > 0) {
      log2P <- function(x) {
        dt = vector()
        for (i in seq_along(x)) {
          dt[[i]] = log(x[[i]])
        }
        return(dt) 
      }
      fitlog2P <- lm( y ~ log2P(x))#log2P model
      AIC <- c(AIC(fitline2P),AIC(fitline3P),AIC(fitexp3P),
               AIC(fitlog2P))
      BIC <- c(BIC(fitline2P),BIC(fitline3P),BIC(fitexp3P),
               BIC(fitlog2P))
    } else{
      AIC <- c(AIC(fitline2P),AIC(fitline3P),AIC(fitexp3P),
               NA)
      BIC <- c(BIC(fitline2P),BIC(fitline3P),BIC(fitexp3P),
               NA)
    }
    
    if(min(y) >0) {
      lmFit_exp2P <- lm(log(y) ~ x)
      coefs_exp2P <- coef(lmFit_exp2P)
      a_exp2P <- exp(coefs_exp2P[1])  #intercept
      b_exp2P <- coefs_exp2P[2]   #slope
      fitexp2P <- nls(y ~ a*exp(b*x), 
                      start = list(a = a_exp2P, b = b_exp2P))
      #exp2P model
      if(min(x) > 0) {
        AIC <- c(AIC(fitline2P),AIC(fitline3P),AIC(fitexp3P),
                 AIC(fitlog2P),AIC(fitexp2P))
        BIC <- c(BIC(fitline2P),BIC(fitline3P),BIC(fitexp3P),
                 BIC(fitlog2P),BIC(fitexp2P))
      } else{
        AIC <- c(AIC(fitline2P),AIC(fitline3P),AIC(fitexp3P),
                 NA,AIC(fitexp2P))
        BIC <- c(BIC(fitline2P),BIC(fitline3P),BIC(fitexp3P),
                 NA,BIC(fitexp2P))
      }
    } else{
      if(min(x) > 0) {
        AIC <- c(AIC(fitline2P),AIC(fitline3P),AIC(fitexp3P),
                 AIC(fitlog2P),NA)
        BIC <- c(BIC(fitline2P),BIC(fitline3P),BIC(fitexp3P),
                 BIC(fitlog2P),NA)
      } else{
        AIC <- c(AIC(fitline2P),AIC(fitline3P),AIC(fitexp3P),
                 NA,NA)
        BIC <- c(BIC(fitline2P),BIC(fitline3P),BIC(fitexp3P),
                 NA,NA)
      }
    }
    
    if(min(x) > 0) {
      if(min(y) > 0) {
        lmFit_power2P <- lm(log(y) ~ log(x))
        coef_power2P <- coef(lmFit_power2P)
        a_power2P <- exp(coef_power2P[1])#intercept
        b_power2P <- coef_power2P[2]#slope
        fitpower2P <- nls(y ~ a*(x^b), 
                          start = list(a = a_power2P, b = b_power2P))
        # power2P model
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
        AIC <- c(AIC(fitline2P),AIC(fitline3P),AIC(fitexp3P),
                 AIC(fitlog2P),AIC(fitexp2P), AIC(fitpower2P),
                 AIC(fitpower3P))
        BIC <- c(BIC(fitline2P),BIC(fitline3P),BIC(fitexp3P),
                 BIC(fitlog2P),BIC(fitexp2P), BIC(fitpower2P),
                 BIC(fitpower3P))
      } else {
        AIC <- c(AIC(fitline2P),AIC(fitline3P),AIC(fitexp3P),
                 AIC(fitlog2P),NA, NA,NA)
        BIC <- c(BIC(fitline2P),BIC(fitline3P),BIC(fitexp3P),
                 BIC(fitlog2P),NA, NA,NA)
      }
    } else{
      AIC <- c(AIC(fitline2P),AIC(fitline3P),AIC(fitexp3P),
               NA,NA, NA,NA)
      BIC <- c(BIC(fitline2P),BIC(fitline3P),BIC(fitexp3P),
               NA,NA, NA,NA)
    }
    Fitting_model = data.frame(
      model = c("line2P","line3P", "exp3P","log2P",
                "exp2P", "power2P","power3P"),
      formula = c("y = ax+b","y = ax^2+bx+c","y = a*exp(b*x)+ c",
                  "y=a*ln(x)+b","y = a*exp(b*x)","y = a*x^b",
                  "y = a*x^b+ c"))
    df <- data.frame(AIC,BIC,Fitting_model)
    df2 <- df[order(AIC,BIC,decreasing=F),]
    selectModel = df2[1,]
    model <- list(selectModel,df2)
    names(model) <- c("BestFitM","AllModel")
    print(model)
  }
}

