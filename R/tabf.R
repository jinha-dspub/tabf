# Global variables to avoid R CMD check NOTEs
utils::globalVariables(c(
  "variables", "values", "n", "prob", "smry", "avg", "std", "p.value",
  "ranks", "order", "tidied", "fit", "dat", ".stratas_tmp", "data",

  "or", "ll", "ul", "pvalue", "keys", "OR95CI", "Variables", "rank", "name",
  "term", "estimate", "conf.low", "conf.high", "Significance", "Label", "Missing N",
  "beta", "se", "Beta95CI", "Subgroup", "OR_CI", "p_int", "std.error"
))

#' Calculate Chi-Square Test Results for Table 1
#'
#' This function computes the p-values from a Pearson's Chi-squared test for multiple categorical variables against a stratifying variable. It is a core utility for generating Phase 3 Table 1 characteristics in epidemiological research.
#'
#' @param dat1 A data frame containing the variables for analysis.
#' @param stratas A character string specifying the name of the stratifying variable (e.g., exposure or disease status).
#' @param catVars A character vector specifying the names of the categorical variables to test against `stratas`.
#'
#' @return A nested tibble containing the `variables` names and the formatted `p.value` (as "<0.001" or a 3-decimal string).
#' @import dplyr tidyr purrr broom labelled stringr
#' @importFrom stats confint.default na.omit sd setNames t.test chisq.test aov
#' @importFrom utils as.roman head
#' @export
#'
#' @examples
#' \dontrun{
#' tab.Chisq(mydata, "exposure_group", c("sex", "age_group"))
#' }
tab.Chisq = function(dat1, stratas, catVars){
  suppressWarnings({
    labelled::var_label(dat1) <- NULL
    dat1 = labelled::remove_labels(dat1)
    dat1 %>%
      select(all_of(stratas), all_of(catVars)) %>%
      pivot_longer(-all_of(stratas), names_to = "variables", values_to ="values")%>%
      group_by(variables, values) %>%
      count(!!sym(stratas)) %>%
      pivot_wider(names_from = all_of(stratas), values_from = n,
                  values_fill = 0) %>%
      ungroup() %>%
      select(-values) %>%
      nest(dat = -variables) %>%
      mutate(
        fit = map(dat, ~{
          m <- as.matrix(.x)
          m <- m[rowSums(m) > 0, colSums(m) > 0, drop = FALSE]
          if (nrow(m) < 2 || ncol(m) < 2) return(list(p.value = NA_real_))
          exp_ok <- tryCatch(
            suppressWarnings(all(chisq.test(m)$expected >= 5)),
            error = function(e) FALSE
          )
          if (exp_ok) chisq.test(m)
          else fisher.test(m, simulate.p.value = TRUE, B = 10000)
        }),
        tidied = map(fit, ~{
          if (is.list(.x) && !inherits(.x, "htest")) {
            tibble::tibble(p.value = .x$p.value)
          } else {
            broom::tidy(.x)
          }
        })
      ) %>%
      unnest(tidied) %>%
      select(variables, p.value) %>%
      mutate(p.value = case_when(
        is.na(p.value)  ~ "",
        p.value < 0.001 ~ "<0.001",
        TRUE            ~ sprintf("%.3f", p.value)
      ))
  })
}

#' Calculate T-Test / ANOVA Results for Table 1
#'
#' This function computes p-values directly for comparing the means of multiple continuous variables across a stratifying variable. It performs Student's t-test for binary stratas and ANOVA for stratas with 3 or more levels. Used in Phase 3 Table 1 construction.
#'
#' @param dat1 A data frame containing the variables for analysis.
#' @param stratas A character string specifying the name of the binary stratifying variable.
#' @param conVars A character vector specifying the names of the continuous variables to test against `stratas`.
#'
#' @return A nested tibble containing the `variables` names and the formatted `p.value` (as "<0.001" or a 3-decimal string).
#' @import dplyr tidyr purrr broom labelled
#' @export
#'
#' @examples
#' \dontrun{
#' tab.Ttest(mydata, "exposure_binary", c("age", "bmi"))
#' }
tab.Ttest =function(dat1, stratas, conVars){
  suppressWarnings({
    labelled::var_label(dat1) <- NULL
    dat1 = labelled::remove_labels(dat1)
    dat1 %>%
      mutate(.stratas_tmp = factor(!!sym(stratas))) %>%
      select(.stratas_tmp, all_of(conVars)) %>%
      pivot_longer(-.stratas_tmp, names_to = "variables", values_to ="values") %>%
      group_by(variables) %>%
      nest() %>%
      mutate(
        p.value = map_dbl(data, function(df) {
          # Extract non-NA levels
          df <- na.omit(df)
          levels <- as.character(unique(df$.stratas_tmp))
          if(length(levels) == 2) {
            g1 <- df$values[df$.stratas_tmp == levels[1]]
            g2 <- df$values[df$.stratas_tmp == levels[2]]
            if(length(g1) == 0 || length(g2) == 0) return(as.numeric(NA))
            return(t.test(g1, g2)$p.value)
          } else if(length(levels) > 2) {
            # Perform ANOVA for 3 or more levels
            res <- aov(values ~ .stratas_tmp, data = df)
            return(summary(res)[[1]][["Pr(>F)"]][1])
          } else {
            return(as.numeric(NA))
          }
        })
      ) %>%
      select(variables, p.value) %>%
      mutate(p.value = ifelse(is.na(p.value), "", ifelse(p.value <0.001, "<0.001", sprintf("%.3f", p.value))))
  })
}

#' Generate Epidemiological Table 1 (Column Percentages)
#'
#' Creates a standardized 'Table 1' (baseline characteristics) commonly used in public health and epidemiological research.
#' For categorical variables, it calculates counts and column percentages. For continuous variables, it calculates means and standard deviations. It automatically appends p-values using chi-square tests and t-tests.
#'
#' @param dat1 A data frame.
#' @param stratas A character string specifying the stratifying variable (e.g., "sex"). If NULL, generates a univariate Table 1 (Total column only).
#' @param catVars A character vector of categorical variable names. Can be NULL.
#' @param conVars A character vector of continuous variable names. Can be NULL.
#'
#' @return A formatted data frame ready for HTML/Markdown rendering with `Variables`, `Values`, strata levels, and `p.value` (if stratas is provided).
#' @import dplyr tidyr purrr stringr broom labelled
#' @export
#'
#' @examples
#' \dontrun{
#' tabf(mydata, "disease", c("sex", "smoking"), c("age", "bmi"))
#' tabf(mydata, NULL, c("sex", "smoking"), c("age", "bmi")) # Univariate
#' }
tabf = function(dat1, stratas = NULL, catVars = NULL, conVars = NULL){
  suppressWarnings({
    labelled::var_label(dat1) <- NULL
    dat1 = labelled::remove_labels(dat1)
    
    is_univariate <- is.null(stratas)
    if(is_univariate) {
      dat1$Total <- "Total"
      stratas <- "Total"
    }
  
    if(!is.null(catVars)) {
      varOrdercat = tibble("variables"=c(catVars)) %>%
        mutate(order = row_number())
    } else {varOrdercat = tibble("variables"=NULL, "order"=NULL)}
    
    if (!is.null(conVars)) {
      varOrdercon = tibble("variables"=c(conVars)) %>%
        mutate(order = row_number())
    } else {varOrdercon = tibble("variables"=NULL, "order"=NULL)}
  
    varOrder = rbind(varOrdercat, varOrdercon) %>% na.omit()
  
    if(!is.null(catVars)){
      catTab = dat1 %>%
        select(all_of(stratas), all_of(catVars)) %>%
        mutate(across(everything(), as.character)) %>%
        pivot_longer(-all_of(stratas), names_to = "variables", values_to ="values")%>%
        group_by( variables, values) %>%
        count(!!sym(stratas)) %>%
        mutate(prob = n/sum(n),
               smry= sprintf("%.0f (%.1f%%)", n, prob*100)
        ) %>%
        select(-n, -prob) %>%
        ungroup() %>%
        pivot_wider(names_from = all_of(stratas), values_from =smry)
    } else {
      catTab = dat1 %>%
        select(all_of(stratas)) %>%
        mutate(nothing =1) %>%
        pivot_longer(-all_of(stratas), names_to = "variables", values_to ="values")%>%
        group_by( variables, values) %>%
        count(!!sym(stratas)) %>%
        mutate(prob = n/sum(n),
               smry= sprintf("%.0f (%.1f%%)", n, prob*100)
        ) %>%
        select(-n, -prob) %>%
        ungroup() %>%
        pivot_wider(names_from = all_of(stratas), values_from =smry) %>%
        filter(variables !="nothing")
    }
  
    if(!is.null(conVars)){
      conTab =
        dat1 %>%
        select(all_of(stratas), all_of(conVars)) %>%
        pivot_longer(-all_of(stratas), names_to = "variables", values_to ="values")%>%
        group_by( !!sym(stratas), variables) %>%
        summarise(avg = mean(values, na.rm =TRUE),
                  std = sd(values, na.rm =TRUE), .groups="drop"
        ) %>%
        mutate(smry  = sprintf("%.1f\u00b1%.1f", avg, std)) %>%
        select(all_of(stratas), variables, smry)%>%
        ungroup() %>%
        pivot_wider(names_from = all_of(stratas), values_from =smry) %>%
        mutate(values ="")
    } else {
      conTab =
        dat1 %>%
        select(all_of(stratas)) %>%
        mutate(nothing =1) %>%
        pivot_longer(-all_of(stratas), names_to = "variables", values_to ="values")%>%
        group_by( !!sym(stratas), variables) %>%
        summarise(avg = mean(values, na.rm =TRUE),
                  std = sd(values, na.rm =TRUE), .groups="drop"
        ) %>%
        mutate(smry  = sprintf("%.1f\u00b1%.1f", avg, std)) %>%
        select(all_of(stratas), variables, smry)%>%
        ungroup() %>%
        pivot_wider(names_from = all_of(stratas), values_from =smry) %>%
        mutate(values ="") %>%
        filter(variables !="nothing")
  
    }
  
    tabDat = rbind(catTab, conTab)
    
    if(is_univariate) {
      tab1 = tabDat %>%
        left_join(varOrder, by = c("variables")) %>%
        arrange(order, values) %>%
        group_by(variables) %>%
        mutate(ranks = row_number()) %>%
        mutate(variables = ifelse(ranks==min(ranks), variables, "")) %>%
        ungroup() %>%
        select(-order, -ranks)%>%
        mutate(values = str_replace(values, "[:digit:]\\.", ""))
    } else {
      if(!is.null(catVars)){
        catPvalue = tab.Chisq(dat1, stratas, catVars)
      } else {
        catPvalue = data.frame(variables=c(""), p.value=c(""))
      }
      
      if(!is.null(conVars)){
        conPvalue = tab.Ttest(dat1, stratas, conVars)
      } else {
        conPvalue = data.frame(variables=c(""), p.value=c(""))
      }
      
      tabPvalue = rbind(catPvalue, conPvalue) %>% na.omit()
    
      tab1 = tabDat %>%
        left_join(tabPvalue, by=c("variables")) %>%
        left_join(varOrder, by = c("variables")) %>%
        arrange(order, values) %>%
        group_by(variables) %>%
        mutate(ranks = row_number()) %>%
        mutate(p.value   = ifelse(ranks==min(ranks), p.value,   "")) %>%
        mutate(variables = ifelse(ranks==min(ranks), variables, "")) %>%
        ungroup() %>%
        select(-order, -ranks)%>%
        mutate(values = str_replace(values, "[:digit:]\\.", ""))
    }
    return(tab1)
  })
}

#' Generate Epidemiological Table 1 (Row Percentages)
#'
#' Similar to `tabf()`, but for categorical variables, it calculates row percentages instead of column percentages.
#' Useful when the strata variable represents an outcome and you want to see the prevalence/incidence across different category levels (e.g., prevalence of disease by age group).
#'
#' @param dat1 A data frame.
#' @param stratas A character string specifying the stratifying variable (e.g., "disease_status"). If NULL, generates univariate summary.
#' @param catVars A character vector of categorical variable names. Can be NULL.
#' @param conVars A character vector of continuous variable names. Can be NULL.
#'
#' @return A formatted data frame ready for HTML/Markdown rendering.
#' @import dplyr tidyr purrr stringr broom labelled
#' @export
#'
#' @examples
#' \dontrun{
#' tabf2(mydata, "disease_outcome", c("age_group", "occupation"))
#' }
tabf2 = function(dat1, stratas = NULL, catVars = NULL, conVars = NULL){
  suppressWarnings({
    labelled::var_label(dat1) <- NULL
    dat1 = labelled::remove_labels(dat1)
    
    is_univariate <- is.null(stratas)
    if(is_univariate) {
      dat1$Total <- "Total"
      stratas <- "Total"
    }
    if(!is.null(catVars)) {
      varOrdercat = tibble("variables"=c(catVars)) %>%
        mutate(order = row_number())
    } else {varOrdercat = tibble("variables"=NULL, "order"=NULL)}
    if (!is.null(conVars)) {
      varOrdercon = tibble("variables"=c(conVars)) %>%
        mutate(order = row_number())
    } else {varOrdercon = tibble("variables"=NULL, "order"=NULL)}
  
    varOrder = rbind(varOrdercat, varOrdercon) %>% na.omit()
  
    if(!is.null(catVars)){
      catTab = dat1 %>%
        select(all_of(stratas), all_of(catVars)) %>%
        mutate(across(everything(), as.character)) %>%
        pivot_longer(-all_of(stratas), names_to = "variables", values_to ="values")%>%
        group_by(!!sym(stratas), variables) %>%
        count(values) %>%
        mutate(prob = n/sum(n),
               smry= sprintf("%.0f (%.1f%%)", n, prob*100)
        ) %>%
        select(-n, -prob) %>%
        ungroup() %>%
        pivot_wider(names_from = all_of(stratas), values_from =smry)
    } else {
      catTab = dat1 %>%
        select(all_of(stratas)) %>%
        mutate(nothing =1) %>%
        pivot_longer(-all_of(stratas), names_to = "variables", values_to ="values")%>%
        group_by(!!sym(stratas), variables) %>%
        count(values) %>%
        mutate(prob = n/sum(n),
               smry= sprintf("%.0f (%.1f%%)", n, prob*100)
        ) %>%
        select(-n, -prob) %>%
        ungroup() %>%
        pivot_wider(names_from = all_of(stratas), values_from =smry) %>%
        filter(variables !="nothing")
    }
  
  
    if(!is.null(conVars)){
      conTab =
        dat1 %>%
        select(all_of(stratas), all_of(conVars)) %>%
        pivot_longer(-all_of(stratas), names_to = "variables", values_to ="values")%>%
        group_by( !!sym(stratas), variables) %>%
        summarise(avg = mean(values, na.rm =TRUE),
                  std = sd(values, na.rm =TRUE), .groups="drop"
        ) %>%
        mutate(smry  = sprintf("%.1f\u00b1%.1f", avg, std)) %>%
        select(all_of(stratas), variables, smry)%>%
        ungroup() %>%
        pivot_wider(names_from = all_of(stratas), values_from =smry) %>%
        mutate(values ="")
    } else {
      conTab =
        dat1 %>%
        select(all_of(stratas)) %>%
        mutate(nothing =1) %>%
        pivot_longer(-all_of(stratas), names_to = "variables", values_to ="values")%>%
        group_by( !!sym(stratas), variables) %>%
        summarise(avg = mean(values, na.rm =TRUE),
                  std = sd(values, na.rm =TRUE), .groups="drop"
        ) %>%
        mutate(smry  = sprintf("%.1f\u00b1%.1f", avg, std)) %>%
        select(all_of(stratas), variables, smry)%>%
        ungroup() %>%
        pivot_wider(names_from = all_of(stratas), values_from =smry) %>%
        mutate(values ="") %>%
        filter(variables !="nothing")
  
    }
  
    tabDat = rbind(catTab, conTab)
  
    if(is_univariate) {
      tab1 = tabDat %>%
        left_join(varOrder, by = c("variables")) %>%
        arrange(order, values) %>%
        group_by(variables) %>%
        mutate(ranks = row_number()) %>%
        mutate(variables = ifelse(ranks==min(ranks), variables, "")) %>%
        ungroup() %>%
        select(-order, -ranks)%>%
        mutate(values = str_replace(values, "[:digit:]\\.", ""))
    } else {
      if(!is.null(catVars)){
        catPvalue = tab.Chisq(dat1, stratas, catVars)
      } else {
        catPvalue = data.frame(variables=c(""), p.value=c(""))
      }
    
      if(!is.null(conVars)){
        conPvalue = tab.Ttest(dat1, stratas, conVars)
      } else {
        conPvalue = data.frame(variables=c(""), p.value=c(""))
      }
    
      tabPvalue = rbind(catPvalue, conPvalue)
    
      tab1 = tabDat %>%
        left_join(tabPvalue, by=c("variables")) %>%
        left_join(varOrder, by = c("variables")) %>%
        arrange(order, values) %>%
        group_by(variables) %>%
        mutate(ranks = row_number()) %>%
        mutate(p.value   = ifelse(ranks==min(ranks), p.value,   "")) %>%
        mutate(variables = ifelse(ranks==min(ranks), variables, "")) %>%
        ungroup() %>%
        select(-order, -ranks)%>%
        mutate(values = str_replace(values, "[:digit:]\\.", ""))
    }
    return(tab1)
  })
}

#' Summarize Logistic Regression Model 
#'
#' Extracts and exponentiates coefficients to calculate Odds Ratios (OR) and 95% Confidence Intervals (CI), along with p-values from a binary logistic regression model.
#'
#' @param mod A fitted logistic regression model object (e.g., from `glm(..., family = binomial)`).
#'
#' @return A matrix containing the exponentiated coefficients (OR), lower and upper 95% confidence bounds, and p-values.
#' @import dplyr broom
#' @export
#'
#' @examples
#' \dontrun{
#' modsmryf(glm(outcome ~ exposure, family = binomial, data = mydata))
#' }
modsmryf=function(mod) {
  cbind(mod$coefficients %>% exp(.), confint.default(mod)%>% exp(.), mod %>% tidy() %>% select(p.value))
}

# Helper function for xlevels
get_xlevels_df <- function(a) {
  if(!any(is.na(a$xlevels))){
    if(length(a$xlevels) == 0){
      tts =  a$model %>% colnames() %>% .[-1]
      t1 = list()
      for (i in 1:length(tts)){
        t1[[tts[i]]] = c("")
      }
    } else{
      t1 = a$xlevels
    }
    bm1 = purrr::map(1:length(t1), function(x){tibble(variables= names(t1)[x], values = t1[[x]])}) %>% do.call(rbind, .)
    return(bm1)
  } else {
    return(data.frame())
  }
}

#' Formatting Logistic Regression Summary (Single Model)
#'
#' Parses a single logistic regression model into a formatted data frame showing variables, levels (values), and a combined `OR (95% CI)` column. It automatically highlights significant findings (p < 0.05) using Markdown/HTML bold tags and marks reference groups.
#'
#' @param a A fitted logistic regression model object.
#'
#' @return A formatted data frame with columns `variables`, `values`, and `OR95CI`.
#' @import dplyr tidyr purrr
#' @export
#'
#' @examples
#' \dontrun{
#' oddf(glm(outcome ~ age + sex, family = binomial))
#' }
oddf=function(a){
  suppressWarnings({
    if(!missing(a)){
      mm = modsmryf(a)
      mm1 = mm%>%
        data.frame() %>%
        setNames(c("or", "ll", "ul", "pvalue")) %>%
        mutate(keys=rownames(mm))

      bm1 = get_xlevels_df(a)

      # 종속변수 제외
      outcome_var <- names(a$model)[1]
      numeric_cols = a$model %>% select(-all_of(outcome_var)) %>% select(where(is.numeric))
      if(ncol(numeric_cols) > 0){
        bm2 = numeric_cols %>% slice(1:2) %>% pivot_longer(everything()) %>% select(variables = name) %>% mutate(values="") %>% unique()
      } else {
        bm2 = data.frame()
      }
      if(nrow(bm1) == 0 & nrow(bm2) == 0) {
          bm0 = data.frame()
      } else if (nrow(bm1) == 0) {
          bm0 = bm2 %>% mutate(keys= paste0(variables, values)) %>% unique()
      } else if (nrow(bm2) == 0) {
          bm0 = bm1 %>% mutate(keys= paste0(variables, values)) %>% unique()
      } else {
          bm0 = rbind(bm1, bm2) %>% mutate(keys= paste0(variables, values)) %>% unique()
      }

      atab= bm0 %>%
        left_join(mm1, by=c("keys")) %>%
        mutate(OR95CI = case_when(
          is.na(or) ~ "<i>1.00 (reference)</i>",
          pvalue < 0.05 ~ sprintf("<b>%.2f (%.2f-%.2f)</b>", round(or, 2), round(ll, 2), round(ul, 2)),
          TRUE ~ sprintf("%.2f (%.2f-%.2f)", round(or, 2), round(ll, 2), round(ul, 2))
        )) %>%
        select(variables, values, OR95CI)
      return(atab)
    } else {
      atab = data.frame("variables"=c(NA), "values"=c(NA), "OR95CI"=c(NA))
      return(atab)
    }
  })
}


#' Mutiple Logistic Regression Models Summary Table
#'
#' Takes multiple logistic regression models (e.g., sequentially adjusted models M1, M2, M3) and combines their `oddf()` outputs into a single, comprehensive comparison table. Essential for Phase 4 of epidemiological research.
#'
#' @param ... Multiple fitted logistic regression model objects.
#' @param model_names An optional character vector of names for the models (e.g., c("Model 1", "Model 2")). If NULL, defaults to "Model.I", "Model.II", etc.
#'
#' @return A data frame juxtaposing the OR (95% CI) of the supplied models side-by-side.
#' @import dplyr purrr stringr
#' @export
#'
#' @examples
#' \dontrun{
#' oddsf(mod1, mod2, mod3, model_names = c("Crude", "+ Age/Sex", "Fully Adjusted"))
#' }
oddsf = function(..., model_names = NULL){
  arglist = list(...)
  tt = map(arglist, oddf) %>%
    reduce(full_join, by=c("variables", "values"))
  vl = c(length(tt)-2)

  if (is.null(model_names)) {
    model_names = paste0("Model.", as.roman(1:vl))
  }

  model_names = head(c(model_names, rep("", max(0, vl - length(model_names)))), vl)

  ys =  arglist[[1]]$formula[2] %>% as.character() %>% str_replace(., "\\=\\=", "of") %>%
    str_replace_all(., '\\"', "")

  tt = tt %>% setNames(c("Variables", "Values", model_names))

  return(tt)
}


#' Multiple Logistic Regression Models Summary Table (HTML)
#'
#' Similar to `oddsf()`, but immediately renders the output as a styled HTML table using the `htmlTable` package.
#'
#' @param ... Multiple fitted logistic regression model objects.
#' @param model_names An optional character vector of names for the models.
#'
#' @return An `htmlTable` object ready for viewing or knitting in an RMarkdown/Quarto document.
#' @import dplyr purrr stringr htmlTable
#' @export
#'
#' @examples
#' \dontrun{
#' oddsTabf(mod1, mod2)
#' }
oddsTabf = function(..., model_names = NULL){
  arglist = list(...)
  mod1 = arglist[[1]]
  tt = map(arglist, oddf) %>%
    reduce(full_join, by=c("variables", "values"))
  vl = c(length(tt)-2)
  if (is.null(model_names)) {
    model_names = paste0("Model.", as.roman(1:vl))
  }
  
  model_names = head(c(model_names, rep("", max(0, vl - length(model_names)))), vl)

  ys =  mod1$formula[2] %>% as.character() %>% str_replace(., "\\=\\=", "of") %>%
    str_replace_all(., '\\"', "")
  tt = tt %>% setNames(c("Variables", "Values", model_names))
  tt %>%  `rownames<-`(NULL) %>%
    group_by(Variables) %>%
    mutate(rank = row_number()) %>%
    mutate(Variables = ifelse(rank == min(rank), Variables, "")) %>%
    mutate(across(starts_with("Model"), ~replace(., is.na(.), ""))) %>%
    ungroup() %>% select(-rank) %>%
    htmlTable::addHtmlTableStyle(align = 'll') %>%
    htmlTable::htmlTable(
      rnames = FALSE,
      caption = sprintf("Table. OR(95%%CI) for %s", ys)

    )

}

#' Raw Logistic Regression Summary Data Frame
#'
#' Similar to `oddf()`, but instead of formatting the output into a single `OR (95% CI)` string, it returns the raw numeric values (`or`, `ll`, `ul`) in separate columns. Useful for downstream plotting (e.g., Forest plots).
#'
#' @param a A fitted logistic regression model object.
#'
#' @return A data frame with columns `variables`, `values`, `or`, `ll`, and `ul`.
#' @import dplyr tidyr purrr
#' @export
#'
#' @examples
#' \dontrun{
#' oddf0(mod) 
#' }
oddf0=function(a){
  suppressWarnings({
    if(!missing(a)){
      mm = modsmryf(a)
      mm1 = mm%>%
        data.frame() %>%
        setNames(c("or", "ll", "ul", "pvalue")) %>%
        mutate(keys=rownames(mm))

      bm1 = get_xlevels_df(a)

      # 종속변수 제외
      outcome_var <- names(a$model)[1]
      numeric_cols = a$model %>% select(-all_of(outcome_var)) %>% select(where(is.numeric))
      if(ncol(numeric_cols) > 0){
        bm2 = numeric_cols %>% slice(1:2) %>% pivot_longer(everything()) %>% select(variables = name) %>% mutate(values="") %>% unique()
      } else {
        bm2 = data.frame()
      }

      if(nrow(bm1) == 0 & nrow(bm2) == 0) {
          bm0 = data.frame()
      } else if (nrow(bm1) == 0) {
          bm0 = bm2 %>% mutate(keys= paste0(variables, values)) %>% unique()
      } else if (nrow(bm2) == 0) {
          bm0 = bm1 %>% mutate(keys= paste0(variables, values)) %>% unique()
      } else {
          bm0 = rbind(bm1, bm2) %>% mutate(keys= paste0(variables, values)) %>% unique()
      }

      atab= bm0 %>%
        left_join(mm1, by=c("keys")) %>%
        mutate(or = ifelse(is.na(or), 1.00, or),
               ll = ifelse(is.na(ll), 1.00, ll),
               ul = ifelse(is.na(ul), 1.00, ul)
        ) %>%
        select(variables, values, or, ll, ul)

      return(atab)
    } else {
      atab = data.frame("variables"=c(NA), "values"=c(NA), "or"=c(NA), "ll"=c(NA), "ul"=c(NA))
      return(atab)
    }
  })
}

#' Summarize Missing Values
#'
#' This function calculates the count and percentage of missing values (NA) for each variable in a data frame.
#'
#' @param dat A data frame.
#'
#' @return A tibble with Variables, Total N, Missing N, and Missing %.
#' @import dplyr tidyr purrr
#' @export
#'
#' @examples
#' \dontrun{
#' missingf(mydata)
#' }
missingf = function(dat) {
  res = purrr::map_df(names(dat), function(var) {
    vec <- dat[[var]]
    n_miss <- sum(is.na(vec))
    n_total <- length(vec)
    pct <- (n_miss / n_total) * 100
    dplyr::tibble(
      Variables = var,
      `Total N` = n_total,
      `Missing N` = n_miss,
      `Missing %` = sprintf("%.1f%%", pct)
    )
  }) %>%
    dplyr::arrange(dplyr::desc(`Missing N`)) %>%
    dplyr::filter(`Missing N` > 0)
    
  if(nrow(res) == 0) {
    message("No missing values found in the dataset.")
  }
  
  return(res)
}

#' Generate Forest Plot from Logistic Regression
#'
#' This function takes a logistic regression model (`glm`) and generates a ggplot2 forest plot for the odds ratios.
#'
#' @param model A logistic regression model (`glm` object with `family = binomial`).
#' @param title A character string for the plot title. Default is "Odds Ratios (95% CI)".
#'
#' @return A ggplot2 object representing the forest plot.
#' @import broom dplyr ggplot2
#' @export
#'
#' @examples
#' \dontrun{
#' mod <- glm(disease ~ age + sex, data = df, family = binomial)
#' forestf(mod)
#' }
forestf = function(model, title = "Odds Ratios (95% CI)") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package \"ggplot2\" must be installed to use this function.", call. = FALSE)
  }
  
  suppressWarnings({
    suppressMessages({
      res = broom::tidy(model, exponentiate = TRUE, conf.int = TRUE) %>%
        dplyr::filter(term != "(Intercept)") %>%
        dplyr::mutate(
          Significance = ifelse(p.value < 0.05, "Significant", "Not Significant"),
          Label = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high)
        ) %>%
        dplyr::arrange(dplyr::desc(estimate)) # Sort by OR
    })
  })
    
  p <- ggplot2::ggplot(res, ggplot2::aes(y = stats::reorder(term, estimate), x = estimate, xmin = conf.low, xmax = conf.high, color = Significance)) +
    ggplot2::geom_pointrange(size = 0.8) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
    ggplot2::scale_color_manual(values = c("Significant" = "#d73027", "Not Significant" = "#4575b4")) +
    ggplot2::labs(title = title, x = "Odds Ratio", y = "Variables") +
    ggplot2::theme_minimal() +
    ggplot2::geom_text(ggplot2::aes(label = Label), vjust = -1, size = 3, show.legend = FALSE) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom"
    )
  
  return(p)
}

# =============================================================================
# LINEAR REGRESSION FUNCTIONS
# =============================================================================

#' Linear Regression Model Summary
#'
#' Extracts beta coefficients with 95% confidence intervals from a linear regression model.
#'
#' @param model A fitted linear regression model object (`lm`).
#'
#' @return A data frame with columns: term, beta, conf.low, conf.high, p.value.
#' @import broom dplyr
#' @export
#'
#' @examples
#' \dontrun{
#' mod <- lm(bmi ~ age + sex, data = df)
#' lmsmryf(mod)
#' }
lmsmryf = function(model) {
  suppressWarnings({
    result <- broom::tidy(model, conf.int = TRUE) %>%
      dplyr::filter(term != "(Intercept)") %>%
      dplyr::select(term, estimate, conf.low, conf.high, p.value) %>%
      dplyr::rename(beta = estimate)
    return(result)
  })
}

#' Format Linear Regression Results as Beta (95% CI)
#'
#' Formats linear regression results similar to `oddf()` but for beta coefficients.
#'
#' @param a A fitted linear regression model object (`lm`).
#'
#' @return A data frame with columns: variables, values, Beta95CI.
#' @import dplyr tidyr purrr
#' @export
#'
#' @examples
#' \dontrun{
#' mod <- lm(bmi ~ age + sex, data = df)
#' betaf(mod)
#' }
betaf = function(a) {
  suppressWarnings({
    if(!missing(a)) {
      mm <- broom::tidy(a, conf.int = TRUE) %>%
        dplyr::filter(term != "(Intercept)") %>%
        dplyr::mutate(
          Beta95CI = sprintf("%.2f (%.2f, %.2f)", estimate, conf.low, conf.high),
          Beta95CI = ifelse(p.value < 0.05,
                           paste0("<b>", Beta95CI, "</b>"),
                           Beta95CI)
        ) %>%
        dplyr::select(term, Beta95CI)

      # Parse term into variables and values
      xlevels <- a$xlevels
      if(length(xlevels) > 0) {
        xlevels_df <- purrr::imap_dfr(xlevels, ~tibble::tibble(
          variables = .y,
          values = .x
        )) %>%
          dplyr::mutate(term = paste0(variables, values))
      } else {
        xlevels_df <- tibble::tibble(variables = character(), values = character(), term = character())
      }

      # Handle numeric variables
      numeric_vars <- names(a$model)[sapply(a$model, is.numeric)]
      numeric_vars <- setdiff(numeric_vars, names(a$model)[1])  # exclude outcome

      if(length(numeric_vars) > 0) {
        numeric_df <- tibble::tibble(
          variables = numeric_vars,
          values = "",
          term = numeric_vars
        )
        xlevels_df <- dplyr::bind_rows(xlevels_df, numeric_df)
      }

      # Add reference categories
      if(length(xlevels) > 0) {
        ref_df <- purrr::imap_dfr(xlevels, ~tibble::tibble(
          variables = .y,
          values = .x[1],
          term = paste0(.y, .x[1]),
          Beta95CI = "<i>0.00 (reference)</i>"
        ))

        result <- xlevels_df %>%
          dplyr::left_join(mm, by = "term") %>%
          dplyr::bind_rows(ref_df %>% dplyr::filter(!term %in% xlevels_df$term)) %>%
          dplyr::arrange(variables, values) %>%
          dplyr::mutate(Beta95CI = ifelse(is.na(Beta95CI), "<i>0.00 (reference)</i>", Beta95CI)) %>%
          dplyr::select(variables, values, Beta95CI)
      } else {
        result <- xlevels_df %>%
          dplyr::left_join(mm, by = "term") %>%
          dplyr::select(variables, values, Beta95CI)
      }

      return(result)
    } else {
      return(data.frame(variables = NA, values = NA, Beta95CI = NA))
    }
  })
}

#' Combined Linear Regression Beta Table
#'
#' Combines multiple linear regression models into a single HTML table showing Beta (95% CI).
#' Similar to `oddsTabf()` but for linear regression models.
#'
#' @param ... Multiple fitted linear regression model objects (`lm`).
#' @param model_names An optional character vector of names for the models.
#'
#' @return An `htmlTable` object ready for viewing.
#' @import dplyr purrr stringr htmlTable
#' @export
#'
#' @examples
#' \dontrun{
#' mod1 <- lm(bmi ~ shift_work, data = df)
#' mod2 <- lm(bmi ~ shift_work + age + sex, data = df)
#' lmTabf(mod1, mod2)
#' }
lmTabf = function(..., model_names = NULL) {
  arglist <- list(...)
  mod1 <- arglist[[1]]

  tt <- purrr::map(arglist, betaf) %>%
    purrr::reduce(dplyr::full_join, by = c("variables", "values"))

  vl <- length(tt) - 2

  if (is.null(model_names)) {
    model_names <- paste0("Model.", utils::as.roman(1:vl))
  }
  model_names <- utils::head(c(model_names, rep("", max(0, vl - length(model_names)))), vl)

  ys <- as.character(mod1$terms[[2]])
  tt <- tt %>% setNames(c("Variables", "Values", model_names))

  tt %>%
    `rownames<-`(NULL) %>%
    dplyr::group_by(Variables) %>%
    dplyr::mutate(rank = dplyr::row_number()) %>%
    dplyr::mutate(Variables = ifelse(rank == min(rank), Variables, "")) %>%
    dplyr::mutate(dplyr::across(dplyr::starts_with("Model"), ~replace(., is.na(.), ""))) %>%
    dplyr::ungroup() %>%
    dplyr::select(-rank) %>%
    htmlTable::addHtmlTableStyle(align = 'll') %>%
    htmlTable::htmlTable(
      rnames = FALSE,
      caption = sprintf("Table. Beta (95%% CI) for %s", ys)
    )
}

# =============================================================================
# INTERACTION ANALYSIS FUNCTIONS
# =============================================================================

#' Run Interaction Model
#'
#' Fits a logistic or linear regression model with an interaction term and returns
#' stratified results along with the p-value for interaction.
#'
#' @param data A data frame containing the variables.
#' @param Y Character string specifying the outcome variable name.
#' @param X Character string specifying the main exposure variable name.
#' @param modifier Character string specifying the effect modifier variable name.
#' @param covars Character vector of covariate names to adjust for.
#' @param family Model family: "binomial" for logistic regression, "gaussian" for linear.
#'
#' @return A list containing: interaction model, stratified models, and p-interaction.
#' @import dplyr broom
#' @export
#'
#' @examples
#' \dontrun{
#' result <- interf(data = df, Y = "depression", X = "work_hours",
#'                  modifier = "sex", covars = c("age"), family = "binomial")
#' }
interf = function(data, Y, X, modifier, covars = NULL, family = "binomial") {
  suppressWarnings({
    # Build formula
    if(!is.null(covars) && length(covars) > 0) {
      covar_str <- paste(covars, collapse = " + ")
      formula_int <- stats::as.formula(paste(Y, "~", X, "*", modifier, "+", covar_str))
      formula_main <- stats::as.formula(paste(Y, "~", X, "+", covar_str))
    } else {
      formula_int <- stats::as.formula(paste(Y, "~", X, "*", modifier))
      formula_main <- stats::as.formula(paste(Y, "~", X))
    }

    # Fit interaction model
    if(family == "binomial") {
      mod_int <- stats::glm(formula_int, data = data, family = stats::binomial())
    } else {
      mod_int <- stats::lm(formula_int, data = data)
    }

    # Extract p-value for interaction term
    int_term <- paste0(X, ":", modifier)
    int_term_alt <- paste0(modifier, ":", X)

    tidy_int <- broom::tidy(mod_int)
    p_int <- tidy_int %>%
      dplyr::filter(grepl(paste0("^", X, ".*:", modifier), term) |
                   grepl(paste0("^", modifier, ".*:", X), term)) %>%
      dplyr::pull(p.value)

    if(length(p_int) == 0) p_int <- NA
    if(length(p_int) > 1) p_int <- min(p_int)  # take minimum if multiple interaction terms

    # Fit stratified models
    modifier_levels <- unique(stats::na.omit(data[[modifier]]))
    stratified_models <- list()

    for(level in modifier_levels) {
      subset_data <- data[data[[modifier]] == level, ]
      if(family == "binomial") {
        stratified_models[[as.character(level)]] <- stats::glm(formula_main,
                                                               data = subset_data,
                                                               family = stats::binomial())
      } else {
        stratified_models[[as.character(level)]] <- stats::lm(formula_main,
                                                              data = subset_data)
      }
    }

    result <- list(
      interaction_model = mod_int,
      stratified_models = stratified_models,
      p_interaction = p_int,
      modifier = modifier,
      modifier_levels = modifier_levels,
      family = family
    )

    class(result) <- c("interf_result", class(result))
    return(result)
  })
}

#' Interaction Analysis OR/Beta Table
#'
#' Creates an HTML table showing stratified OR (or Beta) with p-value for interaction.
#'
#' @param result An object returned by `interf()`.
#' @param digits Number of decimal places for estimates. Default is 2.
#'
#' @return An `htmlTable` object showing stratified results and p-interaction.
#' @import dplyr purrr htmlTable broom
#' @export
#'
#' @examples
#' \dontrun{
#' result <- interf(data = df, Y = "depression", X = "work_hours",
#'                  modifier = "sex", covars = c("age"), family = "binomial")
#' interTabf(result)
#' }
interTabf = function(result, digits = 2) {
  suppressWarnings({
    family <- result$family
    modifier <- result$modifier
    p_int <- result$p_interaction

    # Extract results from stratified models
    strat_results <- purrr::imap_dfr(result$stratified_models, function(mod, level) {
      if(family == "binomial") {
        tidy_res <- broom::tidy(mod, exponentiate = TRUE, conf.int = TRUE)
        tidy_res %>%
          dplyr::filter(term != "(Intercept)") %>%
          dplyr::mutate(
            OR_CI = sprintf(paste0("%.", digits, "f (%.", digits, "f-%.", digits, "f)"),
                           estimate, conf.low, conf.high),
            OR_CI = ifelse(p.value < 0.05, paste0("<b>", OR_CI, "</b>"), OR_CI),
            Subgroup = level
          ) %>%
          dplyr::select(term, Subgroup, OR_CI)
      } else {
        tidy_res <- broom::tidy(mod, conf.int = TRUE)
        tidy_res %>%
          dplyr::filter(term != "(Intercept)") %>%
          dplyr::mutate(
            OR_CI = sprintf(paste0("%.", digits, "f (%.", digits, "f, %.", digits, "f)"),
                           estimate, conf.low, conf.high),
            OR_CI = ifelse(p.value < 0.05, paste0("<b>", OR_CI, "</b>"), OR_CI),
            Subgroup = level
          ) %>%
          dplyr::select(term, Subgroup, OR_CI)
      }
    })

    # Pivot to wide format
    wide_results <- strat_results %>%
      tidyr::pivot_wider(names_from = Subgroup, values_from = OR_CI)

    # Add p-interaction column
    wide_results$p_int <- ""
    wide_results$p_int[1] <- ifelse(is.na(p_int), "-",
                                    ifelse(p_int < 0.001, "<0.001",
                                           sprintf("%.3f", p_int)))

    # Format table
    col_names <- c("Variable", result$modifier_levels, "p-interaction")
    names(wide_results) <- col_names

    estimate_label <- ifelse(family == "binomial", "OR (95% CI)", "Beta (95% CI)")

    wide_results %>%
      htmlTable::addHtmlTableStyle(align = paste0("l", paste(rep("c", ncol(wide_results)-1), collapse = ""))) %>%
      htmlTable::htmlTable(
        rnames = FALSE,
        caption = sprintf("Table. Stratified %s by %s", estimate_label, modifier)
      )
  })
}

#' Forest Plot for Interaction Analysis
#'
#' Creates a forest plot showing stratified estimates from interaction analysis.
#'
#' @param result An object returned by `interf()`.
#' @param title Plot title. Default is "Stratified Analysis".
#'
#' @return A ggplot2 object.
#' @import ggplot2 dplyr purrr broom
#' @export
#'
#' @examples
#' \dontrun{
#' result <- interf(data = df, Y = "depression", X = "work_hours",
#'                  modifier = "sex", covars = c("age"), family = "binomial")
#' interPlotf(result)
#' }
interPlotf = function(result, title = "Stratified Analysis") {
  suppressWarnings({
    family <- result$family
    modifier <- result$modifier
    p_int <- result$p_interaction

    # Extract results from stratified models
    plot_data <- purrr::imap_dfr(result$stratified_models, function(mod, level) {
      if(family == "binomial") {
        tidy_res <- broom::tidy(mod, exponentiate = TRUE, conf.int = TRUE)
      } else {
        tidy_res <- broom::tidy(mod, conf.int = TRUE)
      }
      tidy_res %>%
        dplyr::filter(term != "(Intercept)") %>%
        dplyr::mutate(Subgroup = level)
    })

    # Reference line
    ref_line <- ifelse(family == "binomial", 1, 0)
    x_label <- ifelse(family == "binomial", "Odds Ratio", "Beta Coefficient")

    p_label <- ifelse(is.na(p_int), "p-int: NA",
                     ifelse(p_int < 0.001, "p-int < 0.001",
                            sprintf("p-int = %.3f", p_int)))

    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = estimate, y = term,
                                                  xmin = conf.low, xmax = conf.high,
                                                  color = Subgroup)) +
      ggplot2::geom_pointrange(position = ggplot2::position_dodge(width = 0.5), size = 0.8) +
      ggplot2::geom_vline(xintercept = ref_line, linetype = "dashed", color = "gray50") +
      ggplot2::labs(title = title,
                   subtitle = p_label,
                   x = x_label,
                   y = "Variables",
                   color = modifier) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        legend.position = "bottom"
      )

    return(p)
  })
}
