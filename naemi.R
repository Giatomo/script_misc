---
title: "Naemi"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
library(readxl)
library(tidyverse)
library(rstatix)
```

##

```{r}
df<-read_xlsx("~/Bureau/données_naemi.xlsx")
```

```{r}
df |> group_by(conditions, replicat) |>
  summarise(
    mean = mean(relative_val),
    sd = sd(relative_val),
    error = qt(0.975, df = n()-1)*sd/sqrt(n())) -> summarized

summarized |>
  ggplot(aes(y = mean, x = replicat)) +
  geom_crossbar(aes(ymin = mean-sd, ymax = mean+sd)) +
  facet_wrap(~conditions)

summarized |>
  mutate(all_sd = sd(mean)) |>
  group_by(conditions) |>
  summarise(
    g_mean = mean(mean),
    g_sd = sd(mean),
    error = qt(0.975, df = n()-1)*g_sd/sqrt(n()),
    all_sd) |>
  mutate(ctrl_mean = if_else(conditions == "pBamA",g_mean, NA_real_)) |>
  ungroup() |>
  fill(ctrl_mean, .direction = "downup") |>
  mutate(effect_size = (ctrl_mean - g_mean)/all_sd) -> grand_summ


grand_summ |>
  ggplot(aes(y = g_mean, x = conditions, color = conditions)) +
  geom_crossbar(aes(ymin = g_mean-g_sd, ymax = g_mean+g_sd))

```

```{r}
grand_summ |> kruskal.test(g_mean ~ conditions) |>
  broom::tidy()
```

```{r}
kw_comp. <- nparcomp::nparcomp(data = grand_summ, g_mean ~ conditions)

kw_comp. |> summary()
```
