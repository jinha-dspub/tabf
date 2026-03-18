#' Calculate Chi-Square Test Results for Table 1
#'
#' This function computes the p-values from a Pearson's Chi-squared test for multiple categorical variables against a stratifying variable. It is a core utility for generating Phase 3 Table 1 characteristics in epidemiological research.
#'
#' @param dat1 A data frame containing the variables for analysis.
#' @param stratas A character string specifying the name of the stratifying variable (e.g., exposure or disease status).
#' @param catVars A character vector specifying the names of the categorical variables to test against `stratas`.
#'
#' @return A nested tibble containing the `variables` names and the formatted `p.value` (as "<0.001" or a 3-decimal string).
#' @import dplyr tidyr purrr broom labelled
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
      pivot_wider(names_from = all_of(stratas), values_from =n) %>%
      ungroup() %>%
      select(-values) %>%
      nest(dat = -variables) %>%
      mutate(
        fit = map(dat,
                  ~chisq.test(.x)),
        tidied = map(fit, tidy)
      ) %>%
      unnest(tidied) %>%
      select(variables, p.value) %>%
      mutate(p.value = ifelse(p.value <0.001, "<0.001", sprintf("%.3f", p.value)))
  })
}

#' Calculate T-Test Results for Table 1
#'
#' This function computes p-values using Student's t-test for comparing the means of multiple continuous variables across a two-level stratifying variable. Used in Phase 3 Table 1 construction.
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
        fit   = map(data, function(df) {
          # Extract values for each level of the binary stratas
          levels <- unique(df$.stratas_tmp)
          if(length(levels) != 2) return(NULL) # Safety check
          g1 <- df$values[df$.stratas_tmp == levels[1]]
          g2 <- df$values[df$.stratas_tmp == levels[2]]
          if(length(g1) == 0 || length(g2) == 0) return(NULL)
          t.test(g1, g2)
        }),
        tidied = map(fit, ~if(is.null(.x)) tibble(p.value = as.numeric(NA)) else tidy(.x))
      ) %>%
      unnest(tidied) %>%
      select(variables, p.value) %>%
      mutate(p.value = ifelse(p.value <0.001, "<0.001", sprintf("%.3f", p.value)))
  })
}

#' Generate Epidemiological Table 1 (Column Percentages)
#'
#' Creates a standardized 'Table 1' (baseline characteristics) commonly used in public health and epidemiological research.
#' For categorical variables, it calculates counts and column percentages. For continuous variables, it calculates means and standard deviations. It automatically appends p-values using chi-square tests and t-tests.
#'
#' @param dat1 A data frame.
#' @param stratas A character string specifying the stratifying variable (e.g., "sex" or "disease_status"). If not stratifying, provide a dummy variable.
#' @param catVars A character vector of categorical variable names. Can be NULL.
#' @param conVars A character vector of continuous variable names. Can be NULL.
#'
#' @return A formatted data frame ready for HTML/Markdown rendering with `Variables`, `Values`, strata levels, and `p.value`.
#' @import dplyr tidyr purrr stringr broom labelled
#' @export
#'
#' @examples
#' \dontrun{
#' tabf(mydata, "disease", c("sex", "smoking"), c("age", "bmi"))
#' }
tabf = function(dat1, stratas, catVars = NULL, conVars = NULL){
  suppressWarnings({
    labelled::var_label(dat1) <- NULL
    dat1 = labelled::remove_labels(dat1)
  
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
    return(tab1)
  })
}

#' Generate Epidemiological Table 1 (Row Percentages)
#'
#' Similar to `tabf()`, but for categorical variables, it calculates row percentages instead of column percentages.
#' Useful when the strata variable represents an outcome and you want to see the prevalence/incidence across different category levels (e.g., prevalence of disease by age group).
#'
#' @param dat1 A data frame.
#' @param stratas A character string specifying the stratifying variable (e.g., "disease_status"). 
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
tabf2 = function(dat1, stratas, catVars = NULL, conVars = NULL){
  suppressWarnings({
    labelled::var_label(dat1) <- NULL
    dat1 = labelled::remove_labels(dat1)
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
      
      numeric_cols = a$model %>% select(where(is.numeric))
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
      
      numeric_cols = a$model %>% select(where(is.numeric))
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
