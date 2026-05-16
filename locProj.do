clear all
set more off

import delimited "analysis_data.csv", clear
gen date2 = date(date, "YMD")
format date2 %td
tsset date2, delta(30)

* Define as globals so they work everywhere
global shock   target
global vol     centered_vol
global xvars   lag1_effective_exchange_rate          ///
               lag1_import_price_index_rolling        ///
               lag1_path lag2_path                    ///
               lag1_vix_mps                           ///
               dummy_r1 dummy_r2 dummy_r3 dummy_r4

* Test OLS for one variable, short horizon
locproj log_gdp_sa,             ///
    lcs($shock $vol)             ///
    slags(2)                     ///
    ylags(2)                     ///
    controls($xvars)             ///
    horizon(5)                   ///
    longdiff                     ///
    saving(test_ols, replace)

* Check what locproj saved
use test_ols, clear
describe
list
