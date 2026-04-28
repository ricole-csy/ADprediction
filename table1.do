
* 1. 导入数据
import delimited "alzheimers_disease_data_clean.csv", clear

* 2. 关键步骤：将所有原始变量名转换为小写
rename _all, lower

* 3. 定义变量列表 (此时全部使用小写)
* 连续变量
local cont_vars age bmi alcoholconsumption physicalactivity dietquality sleepquality systolicbp diastolicbp cholesteroltotal cholesterolldl cholesterolhdl cholesteroltriglycerides mmse functionalassessment adl

* 分类变量 (0/1 或 类别变量)
local cat_vars gender ethnicity educationlevel smoking familyhistoryalzheimers cardiovasculardisease diabetes depression headinjury hypertension memorycomplaints behavioralproblems confusion disorientation personalitychanges difficultycompletingtasks forgetfulness

* 3.1 生成你要求的分组变量（用于 Table 1 的 categorical 展示）
*
* Gender: 假设 gender==0 为 Male, gender==1 为 Female
* 如果你的数据相反，把下面两行 label 的顺序互换即可
label define gender_lbl 0 "Male" 1 "Female", replace
label values gender gender_lbl

* Ethnicity: 0=Caucasian, 1=African American, 2=Asian, 3=Other
label define ethnicity_lbl 0 "Caucasian" 1 "African American" 2 "Asian" 3 "Other", replace
label values ethnicity ethnicity_lbl

* Age: 分成 3 个 subgroup（默认 <65, 65–74, >=75；需要的话你可以改 cutpoint）
gen age_grp = .
replace age_grp = 1 if age < 65
replace age_grp = 2 if inrange(age, 65, 74)
replace age_grp = 3 if age >= 75
label define age_grp_lbl 1 "<65" 2 "65–74" 3 ">=75", replace
label values age_grp age_grp_lbl

* EducationLevel: low (0–1) / high (2–3)
gen edu_high = .
replace edu_high = 0 if inlist(educationlevel, 0, 1)
replace edu_high = 1 if inlist(educationlevel, 2, 3)
label define edu_high_lbl 0 "Low (0–1)" 1 "High (2–3)", replace
label values edu_high edu_high_lbl

* 4. 设置结果存储结构
tempname memhold
postfile `memhold' str80 variable str40 non_ad str40 ad str40 or_95ci str40 rr_95ci str20 p_value using "final_lowercase_results.dta", replace

* 5. 连续变量：一行一个变量（mean ± sd）
foreach var in `cont_vars' {
    
    * 检查变量是否存在
    capture confirm variable `var'
    if _rc != 0 continue

    *--- A. 描述性统计 (依据是否在连续变量列表中切换格式) ---
    quietly summarize `var' if diagnosis == 0
    local m0 : display %9.2f r(mean)
    local s0 : display %9.2f r(sd)
    
    quietly summarize `var' if diagnosis == 1
    local m1 : display %9.2f r(mean)
    local s1 : display %9.2f r(sd)
    
    local stats0 "`m0' ± `s0'"
    local stats1 "`m1' ± `s1'"

    *--- B. 计算 OR (Logistic) ---
    quietly logistic diagnosis `var'
    local or_val : display %9.2f exp(_b[`var'])
    local p_val : display %9.3f (2 * normal(-abs(_b[`var']/_se[`var'])))
    local or_ci : display "(" %9.2f exp(_b[`var']-1.96*_se[`var']) "-" %9.2f exp(_b[`var']+1.96*_se[`var']) ")"

    *--- C. 计算 RR (稳健 Poisson 回归) ---
    quietly glm diagnosis `var', family(poisson) link(log) vce(robust)
    local rr_val : display %9.2f exp(_b[`var'])
    local rr_ci : display "(" %9.2f exp(_b[`var']-1.96*_se[`var']) "-" %9.2f exp(_b[`var']+1.96*_se[`var']) ")"

    *--- D. 写入 ---
    post `memhold' ("`var'") ("`stats0'") ("`stats1'") ("`or_val' `or_ci'") ("`rr_val' `rr_ci'") ("`p_val'")
}

* 6. 分类变量：按照要求展开 gender(2 行)、age_grp(3 行)、edu_high(2 行)；
* 其他二分类变量仍按 “>=1 作为 Yes” 输出一行（与原脚本一致）。

* 6.1 Gender（Male/Female 两行；以 Male 为 reference）
    quietly count if diagnosis == 0
    local n0_tot = r(N)
    quietly count if diagnosis == 1
    local n1_tot = r(N)

    * 模型：diagnosis ~ i.gender
    quietly logistic diagnosis i.gender
    local overall_p : display %9.3f chi2tail(e(df_m), e(chi2))
    quietly glm diagnosis i.gender, family(poisson) link(log) vce(robust)

    * Male (ref)
    quietly count if gender == 0 & diagnosis == 0
    local stats0 : display r(N) " (" %4.1f (r(N)/`n0_tot'*100) "%)"
    quietly count if gender == 0 & diagnosis == 1
    local stats1 : display r(N) " (" %4.1f (r(N)/`n1_tot'*100) "%)"
    post `memhold' ("Gender (Male)") ("`stats0'") ("`stats1'") ("1.00 (ref)") ("1.00 (ref)") ("`overall_p'")

    * Female (vs Male)
    quietly count if gender == 1 & diagnosis == 0
    local stats0 : display r(N) " (" %4.1f (r(N)/`n0_tot'*100) "%)"
    quietly count if gender == 1 & diagnosis == 1
    local stats1 : display r(N) " (" %4.1f (r(N)/`n1_tot'*100) "%)"

    quietly logistic diagnosis i.gender
    local or_val : display %9.2f exp(_b[1.gender])
    local p_val  : display %9.3f (2 * normal(-abs(_b[1.gender]/_se[1.gender])))
    local or_ci  : display "(" %9.2f exp(_b[1.gender]-1.96*_se[1.gender]) "-" %9.2f exp(_b[1.gender]+1.96*_se[1.gender]) ")"

    quietly glm diagnosis i.gender, family(poisson) link(log) vce(robust)
    local rr_val : display %9.2f exp(_b[1.gender])
    local rr_ci  : display "(" %9.2f exp(_b[1.gender]-1.96*_se[1.gender]) "-" %9.2f exp(_b[1.gender]+1.96*_se[1.gender]) ")"

    post `memhold' ("Gender (Female)") ("`stats0'") ("`stats1'") ("`or_val' `or_ci'") ("`rr_val' `rr_ci'") ("`p_val'")

* 6.2 Age subgroup（三组；以 <65 为 reference）
    quietly count if diagnosis == 0
    local n0_tot = r(N)
    quietly count if diagnosis == 1
    local n1_tot = r(N)

    quietly logistic diagnosis i.age_grp
    local overall_p : display %9.3f chi2tail(e(df_m), e(chi2))
    quietly glm diagnosis i.age_grp, family(poisson) link(log) vce(robust)

    foreach lvl in 1 2 3 {
        local lbl : label age_grp_lbl `lvl'
        quietly count if age_grp == `lvl' & diagnosis == 0
        local stats0 : display r(N) " (" %4.1f (r(N)/`n0_tot'*100) "%)"
        quietly count if age_grp == `lvl' & diagnosis == 1
        local stats1 : display r(N) " (" %4.1f (r(N)/`n1_tot'*100) "%)"

        if `lvl' == 1 {
            post `memhold' ("Age group (`lbl')") ("`stats0'") ("`stats1'") ("1.00 (ref)") ("1.00 (ref)") ("`overall_p'")
        }
        else {
            quietly logistic diagnosis i.age_grp
            local or_val : display %9.2f exp(_b[`lvl'.age_grp])
            local p_val  : display %9.3f (2 * normal(-abs(_b[`lvl'.age_grp]/_se[`lvl'.age_grp])))
            local or_ci  : display "(" %9.2f exp(_b[`lvl'.age_grp]-1.96*_se[`lvl'.age_grp]) "-" %9.2f exp(_b[`lvl'.age_grp]+1.96*_se[`lvl'.age_grp]) ")"

            quietly glm diagnosis i.age_grp, family(poisson) link(log) vce(robust)
            local rr_val : display %9.2f exp(_b[`lvl'.age_grp])
            local rr_ci  : display "(" %9.2f exp(_b[`lvl'.age_grp]-1.96*_se[`lvl'.age_grp]) "-" %9.2f exp(_b[`lvl'.age_grp]+1.96*_se[`lvl'.age_grp]) ")"

            post `memhold' ("Age group (`lbl')") ("`stats0'") ("`stats1'") ("`or_val' `or_ci'") ("`rr_val' `rr_ci'") ("`p_val'")
        }
    }

* 6.3 Education level（low/high；以 low 为 reference）
    quietly count if diagnosis == 0
    local n0_tot = r(N)
    quietly count if diagnosis == 1
    local n1_tot = r(N)

    quietly logistic diagnosis i.edu_high
    local overall_p : display %9.3f chi2tail(e(df_m), e(chi2))
    quietly glm diagnosis i.edu_high, family(poisson) link(log) vce(robust)

    * Low (ref)
    quietly count if edu_high == 0 & diagnosis == 0
    local stats0 : display r(N) " (" %4.1f (r(N)/`n0_tot'*100) "%)"
    quietly count if edu_high == 0 & diagnosis == 1
    local stats1 : display r(N) " (" %4.1f (r(N)/`n1_tot'*100) "%)"
    post `memhold' ("EducationLevel (Low 0–1)") ("`stats0'") ("`stats1'") ("1.00 (ref)") ("1.00 (ref)") ("`overall_p'")

    * High (vs Low)
    quietly count if edu_high == 1 & diagnosis == 0
    local stats0 : display r(N) " (" %4.1f (r(N)/`n0_tot'*100) "%)"
    quietly count if edu_high == 1 & diagnosis == 1
    local stats1 : display r(N) " (" %4.1f (r(N)/`n1_tot'*100) "%)"

    quietly logistic diagnosis i.edu_high
    local or_val : display %9.2f exp(_b[1.edu_high])
    local p_val  : display %9.3f (2 * normal(-abs(_b[1.edu_high]/_se[1.edu_high])))
    local or_ci  : display "(" %9.2f exp(_b[1.edu_high]-1.96*_se[1.edu_high]) "-" %9.2f exp(_b[1.edu_high]+1.96*_se[1.edu_high]) ")"

    quietly glm diagnosis i.edu_high, family(poisson) link(log) vce(robust)
    local rr_val : display %9.2f exp(_b[1.edu_high])
    local rr_ci  : display "(" %9.2f exp(_b[1.edu_high]-1.96*_se[1.edu_high]) "-" %9.2f exp(_b[1.edu_high]+1.96*_se[1.edu_high]) ")"

    post `memhold' ("EducationLevel (High 2–3)") ("`stats0'") ("`stats1'") ("`or_val' `or_ci'") ("`rr_val' `rr_ci'") ("`p_val'")

* 6.4 Ethnicity subgroup（四组；以 Caucasian(0) 为 reference）
    quietly count if diagnosis == 0
    local n0_tot = r(N)
    quietly count if diagnosis == 1
    local n1_tot = r(N)

    quietly logistic diagnosis i.ethnicity
    local overall_p : display %9.3f chi2tail(e(df_m), e(chi2))
    quietly glm diagnosis i.ethnicity, family(poisson) link(log) vce(robust)

    foreach lvl in 0 1 2 3 {
        local lbl : label ethnicity_lbl `lvl'
        quietly count if ethnicity == `lvl' & diagnosis == 0
        local stats0 : display r(N) " (" %4.1f (r(N)/`n0_tot'*100) "%)"
        quietly count if ethnicity == `lvl' & diagnosis == 1
        local stats1 : display r(N) " (" %4.1f (r(N)/`n1_tot'*100) "%)"

        if `lvl' == 0 {
            post `memhold' ("Ethnicity (`lbl')") ("`stats0'") ("`stats1'") ("1.00 (ref)") ("1.00 (ref)") ("`overall_p'")
        }
        else {
            quietly logistic diagnosis i.ethnicity
            local or_val : display %9.2f exp(_b[`lvl'.ethnicity])
            local p_val  : display %9.3f (2 * normal(-abs(_b[`lvl'.ethnicity]/_se[`lvl'.ethnicity])))
            local or_ci  : display "(" %9.2f exp(_b[`lvl'.ethnicity]-1.96*_se[`lvl'.ethnicity]) "-" %9.2f exp(_b[`lvl'.ethnicity]+1.96*_se[`lvl'.ethnicity]) ")"

            quietly glm diagnosis i.ethnicity, family(poisson) link(log) vce(robust)
            local rr_val : display %9.2f exp(_b[`lvl'.ethnicity])
            local rr_ci  : display "(" %9.2f exp(_b[`lvl'.ethnicity]-1.96*_se[`lvl'.ethnicity]) "-" %9.2f exp(_b[`lvl'.ethnicity]+1.96*_se[`lvl'.ethnicity]) ")"

            post `memhold' ("Ethnicity (`lbl')") ("`stats0'") ("`stats1'") ("`or_val' `or_ci'") ("`rr_val' `rr_ci'") ("`p_val'")
        }
    }

* 6.5 其他二分类变量：输出 “Yes(>=1)” 这一行
foreach var in `cat_vars' {
    if inlist("`var'", "gender", "age", "educationlevel", "ethnicity") continue

    capture confirm variable `var'
    if _rc != 0 continue

    quietly count if diagnosis == 0
    local n0_tot = r(N)
    quietly count if diagnosis == 1
    local n1_tot = r(N)

    quietly count if `var' >= 1 & diagnosis == 0
    local stats0 : display r(N) " (" %4.1f (r(N)/`n0_tot'*100) "%)"
    quietly count if `var' >= 1 & diagnosis == 1
    local stats1 : display r(N) " (" %4.1f (r(N)/`n1_tot'*100) "%)"

    quietly logistic diagnosis `var'
    local or_val : display %9.2f exp(_b[`var'])
    local p_val : display %9.3f (2 * normal(-abs(_b[`var']/_se[`var'])))
    local or_ci : display "(" %9.2f exp(_b[`var']-1.96*_se[`var']) "-" %9.2f exp(_b[`var']+1.96*_se[`var']) ")"

    quietly glm diagnosis `var', family(poisson) link(log) vce(robust)
    local rr_val : display %9.2f exp(_b[`var'])
    local rr_ci : display "(" %9.2f exp(_b[`var']-1.96*_se[`var']) "-" %9.2f exp(_b[`var']+1.96*_se[`var']) ")"

    post `memhold' ("`var'") ("`stats0'") ("`stats1'") ("`or_val' `or_ci'") ("`rr_val' `rr_ci'") ("`p_val'")
}

postclose `memhold'

* 7. 生成并导出 CSV
use "final_lowercase_results.dta", clear
export delimited using "AD_Full_Analysis_Table_Lower_clean.csv", replace
list, clean noobs
