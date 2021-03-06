import delimited "Z:\UKB\GeographicID\Paper Data Extraction\SB_Papers\SW_GWAS\GWAScompare\NoRes.csv"

global phenotypes skincolour haircolour milkconsumption
foreach var in $phenotypes{
	destring `var', ignore("nan") replace
}


*Create residuals
foreach var in $phenotypes{
	reg `var' pc1 pc2 pc3 pc4 pc5 pc6 pc7 pc8 pc9 pc10 pc11 pc12 pc13 pc14 pc15 pc16 pc17 pc18 pc19 pc20 pc21 pc22 pc23 pc24 pc25 pc26 pc27 pc28 pc29 pc30 pc31 pc32 pc33 pc34 pc35 pc36 pc37 pc38 pc39 pc40
	predict `var'_res, residuals

}

* Order variables to assist python sorting
order skincolour_res, after(milkconsumption)
order haircolour_res, after(milkconsumption)
order milkconsumption_res, after(milkconsumption)
order gender, before(pc1)

* Export the csv
export delimited using "Z:\UKB\GeographicID\Paper Data Extraction\SB_Papers\SW_GWAS\GWAScompare\Res.csv", replace
