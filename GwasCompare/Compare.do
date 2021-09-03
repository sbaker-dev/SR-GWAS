import delimited "I:\Work\Genetics\Residuals\Compare two\SkinCompare.csv"

* Coefficent compare
scatter beta_no_res beta_res, xtitle("Residualised beta") ytitle("Non residualised beta")
graph export "I:\Work\Genetics\Residuals\Compare two\Beta comparisions.png", as(png) replace

* P value count
gen no_res_sig = 0
replace no_res_sig = 1 if p_no_res < 0.000000005

gen res_sig = 0
replace res_sig = 1 if p_res < 0.000000005

* No res 	23,597
* Res 		23,443
* No res found an additional 154 more hits at p value 
di 23597 - 23443
