* =============================================================
* local_projections_stata.do
* State-dependent Local Projections
* Three specs: OLS | Rigobon IV | Lewbel robustness
* =============================================================

clear all
set more off

* ssc install locproj
* ssc install ivreg2
* ssc install ivreg2h

import delimited "analysis_data.csv", clear
gen date2 = date(date, "YMD")
format date2 %td
tsset date2, delta(30)

* =============================================================
* PARAMETERS
* =============================================================

local shock   target
local vol     centered_vol
local H       24
local z90     1.645

local controls lag1_target lag2_target                   ///
               lag1_centered_vol lag2_centered_vol        ///
               lag1_effective_exchange_rate               ///
               lag1_import_price_index_rolling            ///
               lag1_vol_mps lag2_vol_mps                  ///
               lag1_path lag2_path                        ///
               lag1_vix_mps                               ///
               dummy_r1 dummy_r2 dummy_r3 dummy_r4

local ronia_iv iv1_dummy_r1 iv1_dummy_r2 iv1_dummy_r3 iv1_dummy_r4 ///
               iv2_dummy_r1 iv2_dummy_r2 iv2_dummy_r3 iv2_dummy_r4

local type1_vars log_gdp_sa log_core_cpi
local type2_vars gb2yt gb5yt gb10yt household_mortgage_rates_sa
local all_vars `type1_vars' `type2_vars'

quietly sum `vol'
local vol_mean = r(mean)
local vol_sd   = r(sd)
local vol_hi   = `vol_mean' + `vol_sd'

* =============================================================
* HELPER PROGRAM: build IRF matrix and plot
* =============================================================

cap program drop irf_plot
program define irf_plot
    args spec var H vol_mean vol_hi z90

    matrix IRF_`spec'_`var' = J(`H'+1, 7, .)

    forval h = 0/`H' {
        local b  = beta_`spec'_`var'_`h'
        local d  = delta_`spec'_`var'_`h'
        local sb = se_b_`spec'_`var'_`h'
        local sd = se_d_`spec'_`var'_`h'

        local est_m = `b' + `d' * `vol_mean'
        local se_m  = sqrt(`sb'^2 + `vol_mean'^2 * `sd'^2)

        local est_h = `b' + `d' * `vol_hi'
        local se_h  = sqrt(`sb'^2 + `vol_hi'^2 * `sd'^2)

        matrix IRF_`spec'_`var'[`=`h'+1', 1] = `h'
        matrix IRF_`spec'_`var'[`=`h'+1', 2] = `est_m'
        matrix IRF_`spec'_`var'[`=`h'+1', 3] = `est_m' - `z90'*`se_m'
        matrix IRF_`spec'_`var'[`=`h'+1', 4] = `est_m' + `z90'*`se_m'
        matrix IRF_`spec'_`var'[`=`h'+1', 5] = `est_h'
        matrix IRF_`spec'_`var'[`=`h'+1', 6] = `est_h' - `z90'*`se_h'
        matrix IRF_`spec'_`var'[`=`h'+1', 7] = `est_h' + `z90'*`se_h'
    }

    matrix colnames IRF_`spec'_`var' = h est_m lo_m hi_m est_h lo_h hi_h
    svmat IRF_`spec'_`var', names(col)

    twoway  (rarea lo_m hi_m h,  color(blue%20))              ///
            (rarea lo_h hi_h h,  color(orange%20))            ///
            (line est_m h,       lcolor(blue) lwidth(medium)) ///
            (line est_h h,       lcolor(orange) lwidth(medium) ///
                                 lpattern(dash)),              ///
            yline(0, lcolor(black) lwidth(thin))               ///
            title("`var' — `spec'")                            ///
            xtitle("Months since shock")                       ///
            ytitle("Response")                                 ///
            legend(order(3 "Mean vol" 4 "High vol (+1SD)"))    ///
            name(irf_`spec'_`var', replace)

    cap drop h est_m lo_m hi_m est_h lo_h hi_h
end

* =============================================================
* 1. OLS — locproj + extract coefficients
* =============================================================

foreach var of local type1_vars {
    locproj `var', shock(`shock') interact(`vol')   ///
        controls(`controls') lags(2)                ///
        horizon(`H') longdiff                       ///
        saving(lp_ols_`var', replace)
}

foreach var of local type2_vars {
    locproj `var', shock(`shock') interact(`vol')   ///
        controls(`controls') lags(2)                ///
        horizon(`H') level                          ///
        saving(lp_ols_`var', replace)
}

foreach var of local all_vars {
    preserve
    use lp_ols_`var', clear
    forval h = 0/`H' {
        quietly sum b_`shock' if horizon == `h'
        scalar beta_ols_`var'_`h'  = r(mean)
        quietly sum b_inter if horizon == `h'
        scalar delta_ols_`var'_`h' = r(mean)
        quietly sum se_`shock' if horizon == `h'
        scalar se_b_ols_`var'_`h'  = r(mean)
        quietly sum se_inter if horizon == `h'
        scalar se_d_ols_`var'_`h'  = r(mean)
    }
    restore
    irf_plot ols `var' `H' `vol_mean' `vol_hi' `z90'
}

* =============================================================
* 2. RIGOBON IV — ivreg2
* =============================================================

foreach var of local type1_vars {
    forval h = 0/`H' {
        cap drop lhs inter_term
        gen lhs        = F`h'.`var' - L.`var'
        gen inter_term = `vol' * `shock'
        ivreg2 lhs `shock' `controls'              ///
            (inter_term = `ronia_iv'),              ///
            gmm2s bw(`=`h'+1') robust
        scalar beta_rig_`var'_`h'  = _b[`shock']
        scalar delta_rig_`var'_`h' = _b[inter_term]
        scalar se_b_rig_`var'_`h'  = _se[`shock']
        scalar se_d_rig_`var'_`h'  = _se[inter_term]
        drop lhs inter_term
    }
    irf_plot rig `var' `H' `vol_mean' `vol_hi' `z90'
}

foreach var of local type2_vars {
    forval h = 0/`H' {
        cap drop lhs inter_term
        gen lhs        = F`h'.`var'
        gen inter_term = `vol' * `shock'
        ivreg2 lhs `shock' `controls'              ///
            (inter_term = `ronia_iv'),              ///
            gmm2s bw(`=`h'+1') robust
        scalar beta_rig_`var'_`h'  = _b[`shock']
        scalar delta_rig_`var'_`h' = _b[inter_term]
        scalar se_b_rig_`var'_`h'  = _se[`shock']
        scalar se_d_rig_`var'_`h'  = _se[inter_term]
        drop lhs inter_term
    }
    irf_plot rig `var' `H' `vol_mean' `vol_hi' `z90'
}

* =============================================================
* 3. LEWBEL — ivreg2h
* =============================================================

foreach var of local type1_vars {
    forval h = 0/`H' {
        cap drop lhs inter_term
        gen lhs        = F`h'.`var' - L.`var'
        gen inter_term = `vol' * `shock'
        ivreg2h lhs `shock' `controls'             ///
            (inter_term = ),                        ///
            gmm2s bw(`=`h'+1') robust
        scalar beta_lew_`var'_`h'  = _b[`shock']
        scalar delta_lew_`var'_`h' = _b[inter_term]
        scalar se_b_lew_`var'_`h'  = _se[`shock']
        scalar se_d_lew_`var'_`h'  = _se[inter_term]
        drop lhs inter_term
    }
    irf_plot lew `var' `H' `vol_mean' `vol_hi' `z90'
}

foreach var of local type2_vars {
    forval h = 0/`H' {
        cap drop lhs inter_term
        gen lhs        = F`h'.`var'
        gen inter_term = `vol' * `shock'
        ivreg2h lhs `shock' `controls'             ///
            (inter_term = ),                        ///
            gmm2s bw(`=`h'+1') robust
        scalar beta_lew_`var'_`h'  = _b[`shock']
        scalar delta_lew_`var'_`h' = _b[inter_term]
        scalar se_b_lew_`var'_`h'  = _se[`shock']
        scalar se_d_lew_`var'_`h'  = _se[inter_term]
        drop lhs inter_term
    }
    irf_plot lew `var' `H' `vol_mean' `vol_hi' `z90'
}

* =============================================================
* 4. COMPARISON PLOTS — OLS vs Rigobon vs Lewbel
* =============================================================

foreach var of local all_vars {

    matrix CMP_`var' = J(`H'+1, 6, .)

    forval h = 0/`H' {
        local b_ols = beta_ols_`var'_`h' + delta_ols_`var'_`h' * `vol_mean'
        local b_rig = beta_rig_`var'_`h' + delta_rig_`var'_`h' * `vol_mean'
        local b_lew = beta_lew_`var'_`h' + delta_lew_`var'_`h' * `vol_mean'

        local sb = se_b_rig_`var'_`h'
        local sd = se_d_rig_`var'_`h'
        local se = sqrt(`sb'^2 + `vol_mean'^2 * `sd'^2)

        matrix CMP_`var'[`=`h'+1', 1] = `h'
        matrix CMP_`var'[`=`h'+1', 2] = `b_ols'
        matrix CMP_`var'[`=`h'+1', 3] = `b_rig'
        matrix CMP_`var'[`=`h'+1', 4] = `b_lew'
        matrix CMP_`var'[`=`h'+1', 5] = `b_rig' - `z90' * `se'
        matrix CMP_`var'[`=`h'+1', 6] = `b_rig' + `z90' * `se'
    }

    matrix colnames CMP_`var' = h ols rig lew lo_rig hi_rig
    svmat CMP_`var', names(col)

    twoway  (rarea lo_rig hi_rig h, color(blue%15))              ///
            (line ols h, lcolor(black) lwidth(medium))           ///
            (line rig h, lcolor(blue)  lwidth(medium))           ///
            (line lew h, lcolor(red)   lwidth(medium)            ///
                         lpattern(dash)),                         ///
            yline(0, lcolor(gray) lwidth(thin))                  ///
            title("`var': OLS vs Rigobon vs Lewbel")             ///
            subtitle("Mean vol state | 90% CI for Rigobon")      ///
            xtitle("Months since shock")                         ///
            ytitle("Response")                                   ///
            legend(order(2 "OLS" 3 "Rigobon IV" 4 "Lewbel"))    ///
            name(cmp_`var', replace)

    cap drop h ols rig lew lo_rig hi_rig
}

* =============================================================
* 5. MARGINAL EFFECTS — ECB Figure 3 style
* =============================================================

local trough_gdp = 0
local min_est    = 0

forval h = 0/`H' {
    local est_m = beta_rig_log_gdp_sa_`h' +             ///
                  delta_rig_log_gdp_sa_`h' * `vol_mean'
    if `est_m' < `min_est' {
        local min_est    = `est_m'
        local trough_gdp = `h'
    }
}

di "GDP trough horizon: `trough_gdp'"

local b_t  = beta_rig_log_gdp_sa_`trough_gdp'
local d_t  = delta_rig_log_gdp_sa_`trough_gdp'
local sb_t = se_b_rig_log_gdp_sa_`trough_gdp'
local sd_t = se_d_rig_log_gdp_sa_`trough_gdp'

matrix ME = J(41, 4, .)
forval i = 1/41 {
    local z     = -2 + (`i'-1) * 0.1
    local xval  = `vol_mean' + `z' * `vol_sd'
    local est   = `b_t' + `d_t' * `xval'
    local se_me = sqrt(`sb_t'^2 + `xval'^2 * `sd_t'^2)
    matrix ME[`i', 1] = `z'
    matrix ME[`i', 2] = `est'
    matrix ME[`i', 3] = `est' - `z90' * `se_me'
    matrix ME[`i', 4] = `est' + `z90' * `se_me'
}

matrix colnames ME = zvol est lo90 hi90
svmat ME, names(col)

twoway  (rarea lo90 hi90 zvol, color(blue%20))              ///
        (line  est  zvol,      lcolor(blue) lwidth(medium)), ///
        yline(0, lcolor(black) lwidth(thin))                 ///
        xline(-0.5, lcolor(blue)  lpattern(dash))            ///
        xline(0,    lcolor(gray)  lpattern(dash))            ///
        xline(1,    lcolor(red)   lpattern(dash))            ///
        title("Marginal effect of vol on GDP transmission")  ///
        subtitle("At trough h=`trough_gdp' | Rigobon IV")   ///
        xtitle("Short-rate volatility (z-score)")            ///
        ytitle("GDP response to 1pp shock")                  ///
        legend(off)                                          ///
        name(marginal_gdp, replace)

cap drop zvol est lo90 hi90

* =============================================================
* 6. COMBINE ALL COMPARISON PLOTS
* =============================================================

graph combine cmp_log_gdp_sa cmp_log_core_cpi              ///
              cmp_gb2yt cmp_gb5yt cmp_gb10yt               ///
              cmp_household_mortgage_rates_sa,              ///
    cols(3)                                                 ///
    title("OLS vs Rigobon IV vs Lewbel — Mean vol state")
local type1_vars log_gdp_sa log_core_cpi
local type2_vars gb2yt gb5yt gb10yt household_mortgage_rates_sa

* =============================================================
* OLS — locproj
* =============================================================

foreach var of local type1_vars {
    locproj `var', shock(`shock') interact(`vol')   ///
        controls(`controls') lags(2)                ///
        horizon(`H') longdiff                       ///
        saving(lp_ols_`var', replace)
}

foreach var of local type2_vars {
    locproj `var', shock(`shock') interact(`vol')   ///
        controls(`controls') lags(2)                ///
        horizon(`H') level                          ///
        saving(lp_ols_`var', replace)
}

* =============================================================
* RIGOBON IV — ivreg2 with pre-constructed instruments
* =============================================================

foreach var of local type1_vars {
    forval h = 0/`H' {
        cap drop lhs inter_term
        gen lhs        = F`h'.`var' - L.`var'
        gen inter_term = `vol' * `shock'

        ivreg2 lhs `shock' `controls'              ///
            (inter_term = `ronia_iv'),              ///
            gmm2s bw(`=`h'+1') robust

        scalar beta_`var'_`h'  = _b[`shock']
        scalar delta_`var'_`h' = _b[inter_term]
        scalar se_b_`var'_`h'  = _se[`shock']
        scalar se_d_`var'_`h'  = _se[inter_term]
        drop lhs inter_term
    }
}

foreach var of local type2_vars {
    forval h = 0/`H' {
        cap drop lhs inter_term
        gen lhs        = F`h'.`var'
        gen inter_term = `vol' * `shock'

        ivreg2 lhs `shock' `controls'              ///
            (inter_term = `ronia_iv'),              ///
            gmm2s bw(`=`h'+1') robust

        scalar beta_`var'_`h'  = _b[`shock']
        scalar delta_`var'_`h' = _b[inter_term]
        scalar se_b_`var'_`h'  = _se[`shock']
        scalar se_d_`var'_`h'  = _se[inter_term]
        drop lhs inter_term
    }
}

* =============================================================
* LEWBEL ROBUSTNESS — ivreg2h
* =============================================================

foreach var of local type1_vars {
    forval h = 0/`H' {
        cap drop lhs inter_term
        gen lhs        = F`h'.`var' - L.`var'
        gen inter_term = `vol' * `shock'

        ivreg2h lhs `shock' `controls'             ///
            (inter_term = ),                        ///
            gmm2s bw(`=`h'+1') robust

        scalar beta_lew_`var'_`h'  = _b[`shock']
        scalar delta_lew_`var'_`h' = _b[inter_term]
        scalar se_b_lew_`var'_`h'  = _se[`shock']
        drop lhs inter_term
    }
}

foreach var of local type2_vars {
    forval h = 0/`H' {
        cap drop lhs inter_term
        gen lhs        = F`h'.`var'
        gen inter_term = `vol' * `shock'

        ivreg2h lhs `shock' `controls'             ///
            (inter_term = ),                        ///
            gmm2s bw(`=`h'+1') robust

        scalar beta_lew_`var'_`h'  = _b[`shock']
        scalar delta_lew_`var'_`h' = _b[inter_term]
        scalar se_b_lew_`var'_`h'  = _se[`shock']
        drop lhs inter_term
    }
}

* =============================================================
* IRF PLOTS — Rigobon IV
* =============================================================

quietly sum `vol'
local vol_mean = r(mean)
local vol_sd   = r(sd)
local vol_hi   = `vol_mean' + `vol_sd'
local z90      = 1.645

foreach var in `type1_vars' `type2_vars' {

    matrix IRF_`var' = J(`H'+1, 7, .)

    forval h = 0/`H' {
        local b  = beta_`var'_`h'
        local d  = delta_`var'_`h'
        local sb = se_b_`var'_`h'
        local sd = se_d_`var'_`h'

        local est_m = `b' + `d' * `vol_mean'
        local se_m  = sqrt(`sb'^2 + `vol_mean'^2 * `sd'^2)
        local lo_m  = `est_m' - `z90' * `se_m'
        local hi_m  = `est_m' + `z90' * `se_m'

        local est_h = `b' + `d' * `vol_hi'
        local se_h  = sqrt(`sb'^2 + `vol_hi'^2 * `sd'^2)
        local lo_h  = `est_h' - `z90' * `se_h'
        local hi_h  = `est_h' + `z90' * `se_h'

        matrix IRF_`var'[`=`h'+1', 1] = `h'
        matrix IRF_`var'[`=`h'+1', 2] = `est_m'
        matrix IRF_`var'[`=`h'+1', 3] = `lo_m'
        matrix IRF_`var'[`=`h'+1', 4] = `hi_m'
        matrix IRF_`var'[`=`h'+1', 5] = `est_h'
        matrix IRF_`var'[`=`h'+1', 6] = `lo_h'
        matrix IRF_`var'[`=`h'+1', 7] = `hi_h'
    }

    matrix colnames IRF_`var' = h est_m lo_m hi_m est_h lo_h hi_h
    svmat IRF_`var', names(col)

    twoway  (rarea lo_m hi_m h,  color(blue%20))              ///
            (rarea lo_h hi_h h,  color(orange%20))            ///
            (line est_m h,       lcolor(blue) lwidth(medium)) ///
            (line est_h h,       lcolor(orange) lwidth(medium) ///
                                 lpattern(dash)),              ///
            yline(0, lcolor(black) lwidth(thin))               ///
            title("`var' response to `shock'")                 ///
            xtitle("Months since shock")                       ///
            ytitle("Response")                                 ///
            legend(order(3 "Mean vol" 4 "High vol (+1SD)"))    ///
            name(irf_`var', replace)

    cap drop h est_m lo_m hi_m est_h lo_h hi_h
}

* =============================================================
* MARGINAL EFFECTS — replicates ECB Figure 3
* =============================================================

* Find trough horizon for GDP
local trough_gdp = 0
local min_est    = 0

forval h = 0/`H' {
    local est_m = beta_log_gdp_sa_`h' +              ///
                  delta_log_gdp_sa_`h' * `vol_mean'
    if `est_m' < `min_est' {
        local min_est    = `est_m'
        local trough_gdp = `h'
    }
}

di "GDP trough horizon: `trough_gdp'"

local b_t  = beta_log_gdp_sa_`trough_gdp'
local d_t  = delta_log_gdp_sa_`trough_gdp'
local sb_t = se_b_log_gdp_sa_`trough_gdp'
local sd_t = se_d_log_gdp_sa_`trough_gdp'

local npts = 41
matrix ME = J(`npts', 4, .)

forval i = 1/`npts' {
    local z     = -2 + (`i'-1) * 0.1
    local xval  = `vol_mean' + `z' * `vol_sd'
    local est   = `b_t' + `d_t' * `xval'
    local se_me = sqrt(`sb_t'^2 + `xval'^2 * `sd_t'^2)
    local lo    = `est' - `z90' * `se_me'
    local hi    = `est' + `z90' * `se_me'

    matrix ME[`i', 1] = `z'
    matrix ME[`i', 2] = `est'
    matrix ME[`i', 3] = `lo'
    matrix ME[`i', 4] = `hi'
}

matrix colnames ME = zvol est lo90 hi90
svmat ME, names(col)

twoway  (rarea lo90 hi90 zvol, color(blue%20))             ///
        (line  est  zvol,      lcolor(blue) lwidth(medium)), ///
        yline(0, lcolor(black) lwidth(thin))                ///
        xline(-0.5, lcolor(blue)   lpattern(dash))          ///
        xline(0,    lcolor(gray)   lpattern(dash))          ///
        xline(1,    lcolor(red)    lpattern(dash))          ///
        title("Marginal effect of vol on GDP transmission") ///
        subtitle("At trough horizon h=`trough_gdp'")        ///
        xtitle("Short-rate volatility (z-score)")           ///
        ytitle("GDP response to 1pp shock")                 ///
        legend(off)                                         ///
        name(marginal_gdp, replace)

cap drop zvol est lo90 hi90
