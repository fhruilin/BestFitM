# BestFitM
This package selects the best-fit model by comparing AIC and BIC. Currently, seven basic fitting models are supported. The seven basic models include one simple linear fitting model (line2P), one quadratic curve model (line3P), one logarithmic model (log2P), two exponential models (exp2P and exp3P), and two power-law models (power2P and power3P).
The formulas of these seven basic fitting models are y = ax+b (line2P), y = ax^2+bx+c (line3P), y = a*exp(b*x)+ c (exp3P), y=a*ln(x)+b (log2P), y = a*exp(b*x) (exp2P), y = a*x^b (power2P), y = a*x^b+ c (power3P), respectively.

# Installation
##require "devtools"

devtools::install_github("fhruilin/BestFitM")

# example
library(BestFitM)
data("mtcars")
bestFitM(data= mtcars, x= "mpg", y = "disp")

![image](https://user-images.githubusercontent.com/50893444/126153624-7b74ff97-08b0-4ea0-8c1e-10938a0d98ea.png)


#Considering that not all data is applicable to the exp3P and power3P models, I give another function, BestFitM2, which contains only five basic fitting models.

# example
bestFitM2(data= mtcars, x= "mpg", y = "disp")
$BestFitM2
       AIC      BIC model        formula
4 354.6943 359.0915 exp2P y = a*exp(b*x)

$AllModel
       AIC      BIC   model        formula
4 354.6943 359.0915   exp2P y = a*exp(b*x)
2 355.6720 361.5349  line3P  y = ax^2+bx+c
3 356.5136 360.9108   log2P    y=a*ln(x)+b
5 359.6060 364.0032 power2P      y = a*x^b
1 363.7164 368.1136  line2P       y = ax+b

#Once the best-fitting model has been selected, you can use the FitM function to look at the other parameters of the model.

# example
FitM(data= mtcars, x= "mpg", y = "disp",model = "line2P")
Call:
lm(formula = y ~ x)

Residuals:
    Min      1Q  Median      3Q     Max 
-103.05  -45.74   -8.17   46.65  153.75 

Coefficients:
            Estimate Std. Error t value Pr(>|t|)    
(Intercept)  580.884     41.740  13.917 1.26e-14 ***
x            -17.429      1.993  -8.747 9.38e-10 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 66.86 on 30 degrees of freedom
Multiple R-squared:  0.7183,	Adjusted R-squared:  0.709 
F-statistic: 76.51 on 1 and 30 DF,  p-value: 9.38e-10
