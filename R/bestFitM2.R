#' @title Best fitting model2
#' 
#' @description
#' \code{bestFitM2} selected the best fitting model from five. 
#' basic fitting models, according to AIC and BIC values.
#' 
#' @details 
#' This function selected the best fitting model from five 
#' basic fitting models.The five basic models include 
#' one simple linear fitting model(line2P), one quadratic curve model(line3P),
#' one logarithmic model(log2P), one exponential models(exp2P), and one power-law models(power2P).
#' This function is selected, when "bestFitM" reported some error messages. 
#' 
#' @param data a dataframe containing the variables x and y.
#' @param x  name of the independent variable x.
#' @param y name of the dependent variable y.
#' @export
#' @return a list with 2 elements:
#' \item{BestFitM2}{the best fitting model}
#' \item{AllModel}{all seven basic models, NA represents the model is not applicable.}
#' @author Ruilin Huang <fhruilin@@163.com>
#' @examples  
#' data("mtcars")
#' bestFitM2(data= mtcars, x= "mpg", y = "disp")
bestFitM2 <- function(data,x= "x",y = "y") {
  dat <- data.frame(data)
  x <- dat[,x]
  y <- dat[,y]
  if(is.na(x) || is.na(y))
  {
    stop("There may be missing values in the variable x or y")
  } else{
    fitline2P <- lm(y ~ x)#line2P model
    fitline3P <- lm(y ~ I(x^2)+ x)#line3P model

    AIC <- c(AIC(fitline2P),AIC(fitline3P))
    BIC <- c(BIC(fitline2P),BIC(fitline3P))
    if(min(x) > 0) {
      log2P <- function(x) {
        dt = vector()
        for (i in seq_along(x)) {
          dt[[i]] = log(x[[i]])
        }
        return(dt) 
      }
      fitlog2P <- lm( y ~ log2P(x))#log2P model
      AIC <- c(AIC(fitline2P),AIC(fitline3P),
               AIC(fitlog2P))
      BIC <- c(BIC(fitline2P),BIC(fitline3P),
               BIC(fitlog2P))
    } else{
      AIC <- c(AIC(fitline2P),AIC(fitline3P),
               NA)
      BIC <- c(BIC(fitline2P),BIC(fitline3P),
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
        AIC <- c(AIC(fitline2P),AIC(fitline3P),
                 AIC(fitlog2P),AIC(fitexp2P))
        BIC <- c(BIC(fitline2P),BIC(fitline3P),
                 BIC(fitlog2P),BIC(fitexp2P))
      } else{
        AIC <- c(AIC(fitline2P),AIC(fitline3P),
                 NA,AIC(fitexp2P))
        BIC <- c(BIC(fitline2P),BIC(fitline3P),
                 NA,BIC(fitexp2P))
      }
    } else{
      if(min(x) > 0) {
        AIC <- c(AIC(fitline2P),AIC(fitline3P),
                 AIC(fitlog2P),NA)
        BIC <- c(BIC(fitline2P),BIC(fitline3P),
                 BIC(fitlog2P),NA)
      } else{
        AIC <- c(AIC(fitline2P),AIC(fitline3P),
                 NA,NA)
        BIC <- c(BIC(fitline2P),BIC(fitline3P),
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
        AIC <- c(AIC(fitline2P),AIC(fitline3P),
                 AIC(fitlog2P),AIC(fitexp2P), AIC(fitpower2P))
        BIC <- c(BIC(fitline2P),BIC(fitline3P),
                 BIC(fitlog2P),BIC(fitexp2P), BIC(fitpower2P))
      } else {
        AIC <- c(AIC(fitline2P),AIC(fitline3P),
                 AIC(fitlog2P),NA, NA)
        BIC <- c(BIC(fitline2P),BIC(fitline3P),
                 BIC(fitlog2P),NA, NA)
      }
    } else{
      AIC <- c(AIC(fitline2P),AIC(fitline3P),
               NA,NA, NA)
      BIC <- c(BIC(fitline2P),BIC(fitline3P),
               NA,NA, NA)
    }
    Fitting_model = data.frame(
      model = c("line2P","line3P", "log2P",
                "exp2P", "power2P"),
      formula = c("y = ax+b","y = ax^2+bx+c",
                  "y=a*ln(x)+b","y = a*exp(b*x)","y = a*x^b"))
    df <- data.frame(AIC,BIC,Fitting_model)
    df2 <- df[order(AIC,BIC,decreasing=F),]
    selectModel = df2[1,]
    model <- list(selectModel,df2)
    names(model) <- c("BestFitM2","AllModel")
    print(model)
  }
}



