/*
Author : Frank Windmeijer, Samuel Baker
Version: 0.1.0
Date   : 28/05/2021
Paper  : SR-GWAS
Purpose: This script tests a single snp on the four models to compare the stata
         code principle to python and show they are identifical as well showing
		 the full process for a single snp.
*/
clear
use "I:\Work\Genetics\Residuals\AT3\Data\ExampleSR.dta"

* Generate phenotypic residuals
qui reg bmi gender age pc1 pc2 pc3 pc4 pc5 pc6 pc7 pc8 pc9 pc10
predict BmiRes, residuals

* Reg Model 1
qui reg bmi rs61747664 gender age pc1 pc2 pc3 pc4 pc5 pc6 pc7 pc8 pc9 pc10

* Extract the OLS estimate and the estimated variance of the snp
matrix snp_est_N = e(b)
matrix var_est_N = e(V)

* Squared T statistics, using the sample size N as the denominator in the variance estimator (RSS/N) instead of RSS/(N-k_2)
scalar w1 = (snp_est_N[1,1]^2) / (var_est_N[1,1]) * e(N) / e(df_r)
di w1

* Reg Model 2
qui reg BmiRes rs61747664

* Extract the OLS estimate and the estimated variance of the snp
matrix snp_est_R = e(b)
matrix var_est_R = e(V)

* Squared T statistics, using the sample size N as the denominator in the variance estimator (RSS/N) instead of RSS/(N-k_2)
scalar w2 = (snp_est_R[1,1]^2) / (var_est_R[1,1]) * e(N) / e(df_r)
scalar pw2 = chi2tail(1, w2)

di w2
di pw2

* Genetic Residuals
qui reg rs61747664 gender age pc1 pc2 pc3 pc4 pc5 pc6 pc7 pc8 pc9 pc10
scalar r2 = e(r2)
predict genres, r

* Reg Model 4
qui reg BmiRes genres

* Extract the OLS estimate of the snp
matrix snp_est_G = e(b)

* Extract the estimated variance of the snp
matrix var_est_G = e(V)

scalar w4 = (snp_est_G[1,1]^2) / (var_est_G[1,1]) * e(N) / e(df_r)
scalar pw4 = chi2tail(1,w4)

di r2
di w4
di pw4


* M1
di w1

* M2
di w2
di pw2

* M4
di r2
di w4
di pw4




