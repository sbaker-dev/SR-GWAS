/*
Author : Frank Windmeijer, Samuel Baker
Version: 0.1.0
Date   : 28/05/2021
Paper  : SR-GWAS
Purpose: This script tests a runs a looping version of the single example
*/


clear
use "I:\Work\Genetics\Residuals\AT3\Data\ExampleSR.dta"


* Adjust Standard errors Program
program sqt_stat
	args est_out
	* Extract the OLS estimate and the estimated variance of the snp
	matrix snp_est = `est_out'(b)
	matrix snp_var = `est_out'(V)
	
	* Squared T statistics, using the sample size N as the denominator in the variance estimator (RSS/N) instead of RSS/(N-k_2)
	scalar cw = (snp_est[1,1]^2) / (snp_var[1,1]) * `est_out'(N) / `est_out'(df_r)

end


* Generate phenotypic residuals
qui reg bmi gender age pc1 pc2 pc3 pc4 pc5 pc6 pc7 pc8 pc9 pc10
predict BmiRes, residuals

* Set snps to loop through
global snps_list rs61747664 rs6677965 rs79578966 rs76744859 rs57498207 rs1539695 rs12091143 rs1473487 rs1171 rs10864532 rs320462 rs34623281 rs8880 rs115448612 rs6426320 rs116396279 rs2247893 rs113336870 rs857825 rs12136312 rs111776064 rs10799828 rs1132185 rs10914250 rs2257763 rs116000650 rs73127292 rs12125026 rs114787228 rs75484590 rs9332660 rs114340315 rs115670157 rs4636400 rs76606177 rs850762 rs841859 rs77953585 rs75369822 rs585273 rs34349246 rs2230678 rs116731477 rs1466423 rs112563471 rs1180967 rs6692377 rs259464 rs114649431 rs17118223 rs1458847 rs72715196 rs698961 rs2873296 rs35877876 rs17038472 rs17377142 rs12087085 rs16855673 rs6676197 rs67018530 rs61769940 rs59666456 rs76438629 rs10789143 rs34938215 rs6683443 rs28653656 rs749917 rs12137480 rs116454739 rs79211218 rs72695340 rs12090635 rs712025 rs215773 rs17479518 rs10923867 rs75752646 rs10921007 rs61775265 rs339569 rs115360270 rs11122383 rs2788887 rs67142165 rs1049434 rs78260708 rs72652923 rs76024980 rs17561193 rs79450364 rs41370845 rs6699951 rs10802390 rs6659226 rs11811071 rs34114668 rs75395157 rs429328

* Create IO Stream
file open results using "C:/Users/Samuel/PycharmProjects/SR-GWAS/Example/Data/Adjusted_output.txt", write replace
file write results "snp" _tab "w1" _tab  "w2" _tab "w2_chi2tail" _tab "w4" _tab "w4_chi2tail" _n

foreach snp of varlist $snps_list{

	* Reg Model 1
	qui reg bmi `snp' gender age pc1 pc2 pc3 pc4 pc5 pc6 pc7 pc8 pc9 pc10
	sqt_stat e
	scalar w1 = cw

	* Reg Model 2
	qui reg BmiRes `snp'
	sqt_stat e
	scalar w2 = cw
	scalar pw2 = chi2tail(1, w2)

	* Genetic Residuals
	qui reg `snp' gender age pc1 pc2 pc3 pc4 pc5 pc6 pc7 pc8 pc9 pc10
	scalar r2 = e(r2)
	predict genres, r

	* Reg Model 4
	qui reg BmiRes genres
	sqt_stat e
	scalar w4 = cw
	scalar pw4 = chi2tail(1,w4)
	drop genres
	
	di w1
	di w4
	di r2
	di

	* Note: This isn't working, not sure why.
	file write results %9s "`snp'" _tab ///
					   %9s "`w1'" _tab ///
					   %9s "`w2'" _tab ///
					   %9s "`pw2'" _tab ///
					   %9s "`w4'" _tab ///
					   %9s "`r2'" _tab ///
					   %9s "`pw4'" _tab _n


}

file close results
program drop sqt_stat
