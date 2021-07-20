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


#Considering that not all data is applicable to the exp3P and power3P models, we give another function, BestFitM2, which contains only five basic fitting models.

# example
bestFitM2(data= mtcars, x= "mpg", y = "disp")

![image](https://user-images.githubusercontent.com/50893444/126153865-b7bc7c49-0b7c-4d6a-8392-31602aa3cfdd.png)


#Once the best-fitting model has been selected, you can use the FitM function to look at the other parameters of the model.

# example
FitM(data= mtcars, x= "mpg", y = "disp",model = "line2P")

![image](https://user-images.githubusercontent.com/50893444/126153914-95642b8c-b347-48c0-acbf-9e5cf40c6b08.png)

