
* Import the data
clear all
set more off
import delimited "C:\Users\Samuel\PycharmProjects\SR-GWAS\Example\Data\CovariantSnp.csv"


* Model 1 **********************************************************************

* regress BMI on G, sex, YoB, PCs
regress bmi rs012 gender age pc1-pc10

* Model 2 **********************************************************************

* Residualise BMI
reg bmi gender age pc1-pc10
predict BMIres, residuals

* regress residualised BMI on G
regress BMIres rs012 

* Model 3 **********************************************************************
* Residualised G
reg rs012 gender age pc1-pc10
predict Gres, residuals

* regress BMI on residualised G
regress bmi Gres 

* Model 4 **********************************************************************
* regress residualised BMI on residualised G
regress BMIres Gres
