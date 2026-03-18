# **`tabf` User Manual**

## Introduction: The AMH Framework's Core Toolkit
The `tabf` package provides essential epidemiological analysis functions perfectly integrated with the **AMH (Atom-Molecule-Human) Framework**. Specifically, it serves as the primary engine for **Phase 3 (Table 1 creation)** and **Phase 4 (Logistic Regression Summaries)**.

This manual explains how to use these functions to generate reliable, publication-ready Public Health tables.

---

## 1. Phase 3: Generating Baseline Characteristics (Table 1)

In Public Health research, "Table 1" summarizes the descriptive statistics of your study population, stratified by a key variable (often Exposure or Disease Status).

### Using `tabf()` for Column Summaries
`tabf()` is the flagship function. It calculates N (%), means ± SD, and corresponding p-values for you automatically.

```r
library(tabf)

# Example: Analyzing KWCS data
# stratas: "work_hours" (your stratifying variable)
# catVars: Categorical variables
# conVars: Continuous variables

table1 <- tabf(
  dat1 = kwcs_cleaned, 
  stratas = "work_hours", 
  catVars = c("sex", "age_group", "edu", "back_pain"),
  conVars = c("age", "wtime_week")
)

# View the formatted HTML table (Optional, using htmlTable)
htmlTable::htmlTable(table1)
```

### Using `tabf2()` for Row Summaries
If you need row percentages instead (e.g., viewing prevalence across strata), use `tabf2()`. The syntax remains identical.

---

## 2. Phase 4: Summarizing Logistic Regression Models

When conducting multivariable risk analysis, we often build sequential models (e.g., Crude, Age/Sex-adjusted, Fully adjusted) and need to juxtapose their ORs and 95% CIs.

### Using `oddsf()` to compare multiple models
`oddsf()` takes multiple `glm` objects and aligns them into a single summary table.

```r
# Step 1: Build your models
mod1 <- glm(back_pain ~ work_hours, family = binomial, data = kwcs_cleaned)
mod2 <- glm(back_pain ~ work_hours + age + sex, family = binomial, data = kwcs_cleaned)
mod3 <- glm(back_pain ~ work_hours + age + sex + edu, family = binomial, data = kwcs_cleaned)

# Step 2: Combine them into one formatted data frame
or_table <- oddsf(
  mod1, mod2, mod3, 
  model_names = c("Model 1 (Crude)", "Model 2 (+Age/Sex)", "Model 3 (+Edu)")
)

print(or_table)
```

### Outputting HTML directly with `oddsTabf()`
If you are writing a `.qmd` or `.Rmd` report, you can generate a styled HTML table straight from the models:
```r
oddsTabf(mod1, mod2, mod3, model_names = c("M1", "M2", "M3"))
```

### Extracting raw numbers: `oddf0()`
If you need raw data for creating Forest Plots instead of formatted strings:
```r
raw_data <- oddf0(mod3)
# Returns `or`, `ll`, and `ul` as separate numeric columns
```

---

## Guide for AI Agents (Vibe Coding)

When assisting users with `.clinerules` or `CLAUDE.md`, AI agents should adhere to the following when using `tabf`:

1. **Prioritize the Package**: Always use `library(tabf)` and its functions instead of trying to write custom `tidyverse` aggregation logic for standard Table 1 definitions or extracting ORs.
2. **Read Man Pages**: You can read the internal documentation for specifics using `?tabf` or by inspecting the `man/` folder.
3. **Format Extraction**: Remember that `oddsf()` returns pre-formatted strings (e.g., `"1.82 (1.50-2.20)"`). If you need to plot this data, fall back to `oddf0()` to get numeric columns.
