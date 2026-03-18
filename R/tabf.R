#' Chisq test
#'
#' @param dat1
#' @param stratas
#' @param catVars
#'
#' @return
#' @export
#'
#' @examples
tab.Chisq = function(dat1, stratas, catVars){
  var_label(dat1) <- NULL
  dat1 = remove_labels(dat1)
  dat1 %>%
    select(stratas, all_of(catVars)) %>%
    pivot_longer(-c(stratas), names_to = "variables", values_to ="values")%>%
    group_by(variables, values) %>%
    count(!!sym(stratas)) %>%
    pivot_wider(names_from = stratas, values_from =n) %>%
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
}

#' T-test
#'
#' @param dat1
#' @param stratas
#' @param conVars
#'
#' @return
#' @export
#'
#' @examples
tab.Ttest =function(dat1, stratas, conVars){
  var_label(dat1) <- NULL
  dat1 = remove_labels(dat1)
  dat1 %>%
    mutate(stratas = !!sym(stratas)) %>%
    select(stratas, conVars) %>%
    pivot_longer(-c(stratas), names_to = "variables", values_to ="values") %>%
    nest(dat = -variables) %>%
    mutate(
      fit   =map(dat, ~t.test(.$values ~ .$stratas)),
      tidied=map(fit, tidy)
    ) %>%
    unnest(tidied) %>%
    select(variables, p.value) %>%
    mutate(p.value = ifelse(p.value <0.001, "<0.001", sprintf("%.3f", p.value)))
}

#' table 1 functon for public health (col sum)
#'
#' @param dat1
#' @param stratas
#' @param catVars
#' @param conVars
#'
#' @return
#' @export
#'
#' @examples
tabf = function(dat1, stratas, catVars, conVars){
  var_label(dat1) <- NULL
  dat1 = remove_labels(dat1)

  if(!missing(catVars)) {
    varOrdercat = tibble("variables"=c(catVars)) %>%
      mutate(order = row_number())
  } else {varOrdercat = tibble("variables"=NULL, "order"=NULL)}
  if (!missing(conVars)) {
    varOrdercon = tibble("variables"=c(conVars)) %>%
      mutate(order = row_number())
  } else {varOrdercon = tibble("variables"=NULL, "order"=NULL)}

  varOrder = rbind(varOrdercat, varOrdercon) %>% na.omit()



  if(!missing(catVars)){
    catTab = dat1 %>%
      select(stratas, all_of(catVars)) %>%
      mutate(across(everything(), as.character)) %>%
      pivot_longer(-c(stratas), names_to = "variables", values_to ="values")%>%
      group_by( variables, values) %>%
      count(!!sym(stratas)) %>%
      mutate(prob = n/sum(n),
             smry= sprintf("%.0f (%.1f%%)", n, prob*100)
      ) %>%
      select(-n, -prob) %>%
      ungroup() %>%
      pivot_wider(names_from = stratas, values_from =smry)
  } else {
    catTab = dat1 %>%
      select(stratas) %>%
      mutate(nothing =1) %>%
      pivot_longer(-c(stratas), names_to = "variables", values_to ="values")%>%
      group_by( variables, values) %>%
      count(!!sym(stratas)) %>%
      mutate(prob = n/sum(n),
             smry= sprintf("%.0f (%.1f%%)", n, prob*100)
      ) %>%
      select(-n, -prob) %>%
      ungroup() %>%
      pivot_wider(names_from = stratas, values_from =smry) %>%
      filter(variables !="nothing")
  }

  if(!missing(conVars)){
    conTab =
      dat1 %>%
      select(stratas, all_of(conVars)) %>%
      pivot_longer(-c(stratas), names_to = "variables", values_to ="values")%>%
      group_by( !!sym(stratas), variables) %>%
      summarise(avg = mean(values, na.rm =TRUE),
                std = sd(values, na.rm =TRUE)
      ) %>%
      mutate(smry  = sprintf("%.1f\u00b1%.1f", avg, std)) %>%
      select(stratas, variables, smry)%>%
      ungroup() %>%
      pivot_wider(names_from = stratas, values_from =smry) %>%
      mutate(values ="")
  } else {
    conTab =
      dat1 %>%
      select(stratas) %>%
      mutate(nothing =1) %>%
      pivot_longer(-c(stratas), names_to = "variables", values_to ="values")%>%
      group_by( !!sym(stratas), variables) %>%
      summarise(avg = mean(values, na.rm =TRUE),
                std = sd(values, na.rm =TRUE)
      ) %>%
      mutate(smry  = sprintf("%.1f\u00b1%.1f", avg, std)) %>%
      select(stratas, variables, smry)%>%
      ungroup() %>%
      pivot_wider(names_from = stratas, values_from =smry) %>%
      mutate(values ="") %>%
      filter(variables !="nothing")

  }

  tabDat = rbind(catTab, conTab)
  if(!missing(catVars)){
    catPvalue =
      dat1 %>%
      select(stratas, catVars) %>%
      mutate(across(everything(), as.character)) %>%
      pivot_longer(-c(stratas), names_to = "variables", values_to ="values")%>%
      group_by(variables, values) %>%
      count(!!sym(stratas)) %>%
      pivot_wider(names_from = stratas, values_from =n) %>%
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
  } else {
    catPvalue = data.frame(variables=c(""), p.value=c(""))
  }
  if(!missing(conVars)){
    conPvalue=dat1 %>%
      mutate(stratas = !!sym(stratas)) %>%
      select(stratas, all_of(conVars)) %>%
      pivot_longer(-c(stratas), names_to = "variables", values_to ="values") %>%
      nest(dat = -variables) %>%
      mutate(
        fit   =map(dat, ~t.test(.$values ~ .$stratas)),
        tidied=map(fit, tidy)
      ) %>%
      unnest(tidied) %>%
      select(variables, p.value) %>%
      mutate(p.value = ifelse(p.value <0.001, "<0.001", sprintf("%.3f", p.value)))
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
}  %>% suppressWarnings()





#' Table 1 for public health (row sum)
#'
#' @param dat1
#' @param stratas
#' @param catVars
#' @param conVars
#'
#' @return
#' @export
#'
#' @examples
tabf2 = function(dat1, stratas, catVars, conVars){
  var_label(dat1) <- NULL
  dat1 = remove_labels(dat1)
  if(!missing(catVars)) {
    varOrdercat = tibble("variables"=c(catVars)) %>%
      mutate(order = row_number())
  } else {varOrdercat = tibble("variables"=NULL, "order"=NULL)}
  if (!missing(conVars)) {
    varOrdercon = tibble("variables"=c(conVars)) %>%
      mutate(order = row_number())
  } else {varOrdercon = tibble("variables"=NULL, "order"=NULL)}

  varOrder = rbind(varOrdercat, varOrdercon) %>% na.omit()

  if(!missing(catVars)){
    catTab = dat1 %>%
      select(stratas, all_of(catVars)) %>%
      mutate(across(everything(), as.character)) %>%
      pivot_longer(-c(stratas), names_to = "variables", values_to ="values")%>%
      group_by(!!sym(stratas), variables) %>%
      count(values) %>%
      mutate(prob = n/sum(n),
             smry= sprintf("%.0f (%.1f%%)", n, prob*100)
      ) %>%
      select(-n, -prob) %>%
      ungroup() %>%
      pivot_wider(names_from = stratas, values_from =smry)
  } else {
    catTab = dat1 %>%
      select(stratas) %>%
      mutate(nothing =1) %>%
      pivot_longer(-c(stratas), names_to = "variables", values_to ="values")%>%
      group_by(!!sym(stratas), variables) %>%
      count(values) %>%
      mutate(prob = n/sum(n),
             smry= sprintf("%.0f (%.1f%%)", n, prob*100)
      ) %>%
      select(-n, -prob) %>%
      ungroup() %>%
      pivot_wider(names_from = stratas, values_from =smry) %>%
      filter(variables !="nothing")
  }


  if(!missing(conVars)){
    conTab =
      dat1 %>%
      select(stratas, all_of(conVars)) %>%
      pivot_longer(-c(stratas), names_to = "variables", values_to ="values")%>%
      group_by( !!sym(stratas), variables) %>%
      summarise(avg = mean(values, na.rm =TRUE),
                std = sd(values, na.rm =TRUE)
      ) %>%
      mutate(smry  = sprintf("%.1f\u00b1%.1f", avg, std)) %>%
      select(stratas, variables, smry)%>%
      ungroup() %>%
      pivot_wider(names_from = stratas, values_from =smry) %>%
      mutate(values ="")
  } else {
    conTab =
      dat1 %>%
      select(stratas) %>%
      mutate(nothing =1) %>%
      pivot_longer(-c(stratas), names_to = "variables", values_to ="values")%>%
      group_by( !!sym(stratas), variables) %>%
      summarise(avg = mean(values, na.rm =TRUE),
                std = sd(values, na.rm =TRUE)
      ) %>%
      mutate(smry  = sprintf("%.1f\u00b1%.1f", avg, std)) %>%
      select(stratas, variables, smry)%>%
      ungroup() %>%
      pivot_wider(names_from = stratas, values_from =smry) %>%
      mutate(values ="") %>%
      filter(variables !="nothing")

  }

  tabDat = rbind(catTab, conTab)

  if(!missing(catVars)){
    catPvalue =
      dat1 %>%
      select(stratas, catVars) %>%
      mutate(across(everything(), as.character)) %>%
      pivot_longer(-c(stratas), names_to = "variables", values_to ="values")%>%
      group_by(variables, values) %>%
      count(!!sym(stratas)) %>%
      pivot_wider(names_from = stratas, values_from =n) %>%
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
  } else {
    catPvalue = data.frame(variables=c(""), p.value=c(""))
  }

  if(!missing(conVars)){
    conPvalue=dat1 %>%
      mutate(stratas = !!sym(stratas)) %>%
      select(stratas, all_of(conVars)) %>%
      pivot_longer(-c(stratas), names_to = "variables", values_to ="values") %>%
      nest(dat = -variables) %>%
      mutate(
        fit   =map(dat, ~t.test(.$values ~ .$stratas)),
        tidied=map(fit, tidy)
      ) %>%
      unnest(tidied) %>%
      select(variables, p.value) %>%
      mutate(p.value = ifelse(p.value <0.001, "<0.001", sprintf("%.3f", p.value)))
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
}  %>% suppressWarnings()

#' logistic regression model summary
#'
#' @param mod
#'
#' @return
#' @export
#'
#' @examples
modsmryf=function(mod) {
  cbind(mod$coefficients %>% exp(.), confint.default(mod)%>% exp(.), mod %>% tidy() %>% select(p.value))}

#' logistic regression model summary
#'
#' @param a
#'
#' @return
#' @export
#'
#' @examples
oddf=function(a){
  if(!missing(a)){
    mm = modsmryf(a)
    mm1 = mm%>%
      data.frame() %>%
      setNames(c("or", "ll", "ul", "pvalue")) %>%
      mutate(keys=rownames(mm))
    if(!any(is.na(a$xlevels))){
      if(length(a$xlevels) ==0){
        tts =  a$model %>% colnames() %>% .[-1]
        t1 = list()
        for (i in 1:c(length(tts))){
          t1[[tts[i]]] = c("")
        }

      } else{
        t1 = a$xlevels
      }

      bm1 = map(1:length(t1),function(x){tibble(variables= names(t1)[x], values = t1[[x]])}) %>% do.call(rbind, .)
    } else {
      t1 = data.frame();bm1=data.frame()
    }
    numeric_cols = a$model %>% select(where(is.numeric))
    if(ncol(numeric_cols) > 0){
      bm2 = numeric_cols %>% slice(1:2) %>% pivot_longer(everything()) %>% select(variables = name) %>% mutate(values="") %>% unique()
    } else {
      bm2 = data.frame()
    }
    bm0 = rbind(bm1, bm2) %>% mutate(keys= paste0(variables, values)) %>%
      unique()

    atab= bm0 %>%
      left_join(mm1, by=c("keys")) %>%
      mutate(OR95CI = case_when(
        is.na(or) ~ "<i>1.00 (reference)</i>",
        pvalue < 0.05 ~ sprintf("<b>%.2f (%.2f-%.2f)</b>", round(or, 2), round(ll, 2), round(ul, 2)),
        TRUE ~ sprintf("%.2f (%.2f-%.2f)", round(or, 2), round(ll, 2), round(ul, 2))
      )) %>%
      #mutate(values = case_when(
      #  pvalue <0.05 ~ sprintf("<b>%s</b>", values),
      #  TRUE ~ values
      #)) %>%
      select(variables, values, OR95CI)
    return(atab)
  } else {
    atab = data.frame("variables"=c(NA), "values"=c(NA), "OR95CI"=c(NA))
    return(atab)
  }
} %>% suppressWarnings()


#' Two or more logistic regression model summary table
#'
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
oddsf = function(..., model_names = NULL){
  arglist = list(...)
  tt = map(arglist, oddf) %>%
    reduce(full_join, by=c("variables", "values"))
  vl = c(length(tt)-2)

  if (is.null(model_names)) {
    model_names = paste0("Model.", as.roman(1:vl))
  }

  # Ensure that the number of models matches the number of columns
  model_names = head(c(model_names, rep("", vl - length(model_names))), vl)

  ys =  arglist[[1]]$formula[2] %>% as.character() %>% str_replace(., "\\=\\=", "of") %>%
    str_replace_all(., '\\"', "")

  tt = tt %>% setNames(c("Variables", "Values", model_names))

  return(tt)
}


#' Two or more logistic regression model summary table with html
#'
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
oddsTabf = function(..., model_names = NULL){
  arglist = list(...)
  mod1 = arglist[[1]]
  tt = map(arglist, oddf) %>%
    reduce(full_join, by=c("variables", "values"))
  vl = c(length(tt)-2)
  if (is.null(model_names)) {
    model_names = paste0("Model.", as.roman(1:vl))
  }
  ys =  mod1$formula[2] %>% as.character() %>% str_replace(., "\\=\\=", "of") %>%
    str_replace_all(., '\\"', "")
  tt = tt %>% setNames(c("Variables", "Values", model_names))
  tt %>%  `rownames<-`(NULL) %>%
    group_by(Variables) %>%
    mutate(rank = row_number()) %>%
    mutate(Variables = ifelse(rank == min(rank), Variables, "")) %>%
    mutate_at(., vars(starts_with("Model")), ~replace(., is.na(.), "")) %>%
    ungroup() %>% select(-rank) %>%
    addHtmlTableStyle(align = 'll') %>%
    htmlTable(
      rnames = FALSE,
      caption = sprintf("Table. OR(95%%CI) for %s", ys)

    )

}






#' Two or more logistic regression model data.frame
#'
#' @param a
#'
#' @return
#' @export
#'
#' @examples
oddf0=function(a){
  if(!missing(a)){
    mm = modsmryf(a)
    mm1 = mm%>%
      data.frame() %>%
      setNames(c("or", "ll", "ul", "pvalue")) %>%
      mutate(keys=rownames(mm))
    if(!any(is.na(a$xlevels))){
      t1 = a$xlevels
      bm1 = map(1:length(t1),function(x){tibble(variables= names(t1)[x], values = t1[[x]])}) %>% do.call(rbind, .)
    } else {
      t1 = data.frame();bm1=data.frame()
    }
    numeric_cols = a$model %>% select(where(is.numeric))
    if(ncol(numeric_cols) > 0){
      bm2 = numeric_cols %>% slice(1:2) %>% pivot_longer(everything()) %>% select(variables = name) %>% mutate(values="") %>% unique()
    } else {
      bm2 = data.frame()
    }
    bm0 = rbind(bm1, bm2) %>% mutate(keys= paste0(variables, values))

    atab= bm0 %>%
      left_join(mm1, by=c("keys")) %>%
      mutate(or = ifelse(is.na(or), 1.00, or),
             ll = ifelse(is.na(ll), 1.00, ll),
             ul = ifelse(is.na(ul), 1.00, ul)
      ) %>%
      select(variables, values, or, ll, ul)

    return(atab)
  } else {
    atab = data.frame("variables"=c(NA), "values"=c(NA), "OR95CI"=c(NA))
    return(atab)
  }
} %>% suppressWarnings()
