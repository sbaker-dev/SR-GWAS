/*
Author : Stephanie von Hinke, Samuel Baker
Version: 0.2.0
Date   : 17/05/2021
Paper  : SR-GWAS
Purpose: This script tests a single snp on the four models to compare the stata 
         code principle to python and show they are identifical as well showing 
		 the full process for a single snp.
*/

* Import the data
clear all
set more off
import delimited "C:\Users\Samuel\PycharmProjects\SR-GWAS\Example\Data\CovariantSnp.csv"

log using "C:\Users\Samuel\PycharmProjects\SR-GWAS\Example\Data\Output.log", replace

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

log close
