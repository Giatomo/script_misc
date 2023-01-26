---
  title: "Curve Fitting"
output: html_notebook
---
  
  $$
  P(t) = p_{min} + \frac{p_{max}-p_{min}}{1 + e^{-r.(t-s)}}
$$
  
  {r}
richard <- function(t, p_max, p_min, r, s){
  return(p_min + (p_max - p_min)  / (1 + exp(-r * (t - s))))
}


$$ \begin{align*}
P(t) &= p_{min} + \frac{p_{max} - p_{min}}
{1 + e^{-r \cdot (t - s)}} \\

P'(t) &= -\frac{(p_{min} - p_{max}) \cdot r \cdot e^{r \cdot (t - s)}}
              {(e^{r \cdot (t - s)} + 1)^2} \\

P''(t) &= -\frac{(p_{min} - p_{max}) \cdot r^2 \cdot (e^{r \cdot (t - s)} - 1) \cdot e^{r \cdot (t - s)} }
               {(e^{r \cdot (t - s)} + 1)^3} \\

P'''(t) &= -\frac{(p_{min} - p_{max}) \cdot r^3 \cdot e^{r \cdot (t - s)} \cdot (e^{2r \cdot (t - s)} - 4e^{r \cdot (t - s}) + 1}
{(e^{r \cdot (t - s)} + 1)^4} \\
\end{align*} $$
  
  Replace $r$ in $P(t)$ by $r_{max}$ which is the extremum of $P'(t)$ i.e. the value of $P'(t)$ where $P''(t) = 0$
  
  $$
  -\frac{(p_{min} - p_{max}) . r^2 . (e^{r . (t - s)} - 1) . e^{r . (t - s)} } {(e^{r . (t - s)} + 1)^3} \stackrel{?}{=} 0
$$
  
  Case 1 : $p_{min} - p_{max} = 0$
  
  $$
  \begin{align*}
p_{min} - p_{max} &= 0 \\
p_{min} &= p_{max}
\end{align*}
$$
  
  Result in no growth

Case 2 : $r^2 = 0$
  
  $$
  \begin{align*}
r^2 &= 0 \\
r &= 0
\end{align*}
$$
  
  Result in no growth

Case 3 : $e^{r(t-s)} = 0$
  
  $$
  \begin{align*}
e^{r \cdot (t - s)} &= 0
\end{align*}
$$
  
  Impossible

Case 4 : $e^{r(t-s)} -1 = 0$
  
  $$
  \begin{align*}
e^{r \cdot (t - s)} - 1 &= 0 \\
e^{r \cdot (t - s)} &= 1 \\
r \cdot (t - s) &= 0 \\
\end{align*}
$$
  
  $$
  r = 0
$$
  
  Result in no growth

$$
  t = s
$$
  
  Ok

Determine whether the extremum is a minima, a maxima or 0

Replace $t$ by $s$ in $P'''(t)$

$$ \begin{align*}
P'''(t) &= -\frac{(p_{min} - p_{max}) \cdot r^3 \cdot e^{r \cdot (t - s)} \cdot (e^{2r \cdot (t - s)} - 4e^{r \cdot (t - s)} + 1)}
{(e^{r \cdot (t - s)} + 1)^4} \\

P'''(s) &= -\frac{(p_{min} - p_{max}) \cdot r^3 \cdot e^{r \cdot (s - s)} \cdot (e^{2r \cdot (s - s)} - 4e^{r \cdot (s - s)} + 1)}
                {(e^{r \cdot (s - s)} + 1)^4} \\

P'''(s) &= -\frac{(p_{min} - p_{max}) \cdot r^3 \cdot 1 \cdot (1 - 4 + 1)}
{(1 + 1)^4} \\

P'''(s) &= -\frac{-2r^3 \cdot (p_{min} - p_{max}) }
                {16} \\
\end{align*} $$

Case 1 : $r > 0$

$$
P'''(S) < 0 \\
P'(S) \rightarrow maximum
$$

Case 2 : $r < 0$

$$
P'''(S) > 0 \\
P'(S) \rightarrow minimum
$$

Case 3 : $r = 0$

$$
P'''(S) < 0 \\
P'(S) \rightarrow 0
$$

Let's find the value of the extremum i.e. $P'(s)$

$$ \begin{align*}

P'(t) &= -\frac{(p_{min} - p_{max}) \cdot r \cdot e^{r \cdot (t - s)}}
{(e^{r \cdot (t - s)} + 1)^2} \\

P'(s) &= -\frac{(p_{min} - p_{max}) \cdot r \cdot e^{r \cdot (s - s)}}
              {(e^{r \cdot (s - s)} + 1)^2} \\

\frac{-4P'(s)}{(p_{min} - p_{max})} &= r \\

\frac{-4r_{max}}{(p_{min} - p_{max})} &= r \\


\end{align*} $$
  
  Let's then replace $r$ in $P(t)$ by its value

$$
\begin{align*}
 P(t) &= p_{min} + \frac{p_{max}-p_{min}}{1 + e^{-r.(t-s)}} \\
 P(t) &= p_{min} + \frac{p_{max}-p_{min}}{1 + e^{4r_{max}.(t-s)/p_{min}- p_{max}}}
\end{align*}
$$

{r}
richard2 <- function(t, p_max, p_min, r_max, s){
  return(p_min + (p_max - p_min)  / (1 + exp(4 * r_max * (t - s)/(p_min - p_max))))
}

We can now use this new equation to fit our date using a genetic algorithm

{r}
fit_richard2 <- function(y, t) {
  
  max_y <- max(y)
  min_y <- min(y)
  r_boundary <- plyr::round_any(max_y-min_y,
                      10^(floor(log10(max_y-min_y))),
                      f = ceiling)
  s_boundary <- max(t)
  
  GA::ga(type = "real-valued",
         fitness = \(p) -sum((y - richard2(t = t,  p[1], p[2], p[3], p[4]))^2),
         lower   = c(
           p_max = max_y - 0.01*max_y,
           p_min = min_y - 0.01*max_y,
           r_max = -r_boundary,
           s     = 0),
         upper = c(
           p_max = max_y + 0.01*max_y,
           p_min = min_y + 0.01*max_y,
           r_max = r_boundary,
           s     = s_boundary),
         names = c(
           "p_max",
           "p_min",
           "r_max",
           "s"),
         popSize = 100,
         monitor = FALSE)@solution
}

Load some sample data

{r}
test_data <- readxl::read_excel(
  path = "/home/thomas/Téléchargements/Archive/21-07-02 GFP Gain 100/21-07-02-GFP gain 100.xlsx",
  sheet = "GFP 100 raw",
  na = "OVRFLW")

{r}
library(glue)
library(ggtext)
library(tidyverse)
test_data |> 
  select(Time, matches("^[A-Za-z]{1}\\d{1,2}")) |> 
  mutate(Time = (Time - first(Time)) |> as.numeric()/3600) |> 
  pivot_longer(cols = matches("^[A-Za-z]{1}\\d{1,2}"),
               names_to = "well")  |> 
  drop_na(value) ->  cleaned

cleaned|> 
  group_by(well) |> 
  nest(data = c(Time, value)) |> 
  mutate(fit = purrr::map(data, ~fit_richard2(.$value, .$Time))) |> 
  unnest_wider(fit, names_sep = "_") |> 
  rename(p_max = fit_1, p_min = fit_2,r_max = fit_3,s = fit_4) |>
  select(-data) -> parameters
  

parameters |> 
  mutate(annot = glue("rate<sub>max</sub> {sprintf('%0.2f', r_max)}")) -> annot

cleaned |>
  inner_join(parameters) |> 
  arrange(well) |> 
  rowwise() |> 
  mutate(fit = richard2(Time, p_max, p_min, r_max, s),
         max_r_t = ((log(abs(-(-p_min+(3*p_max)) / (p_min+p_max))) * (p_min-p_max))/(4*r_max))+s,
#         max_r_y = richard2(max_r_t, p_max, p_min, r_max, s),
#         test = max_r_y - ((p_max - p_min)/2 + p_min),
         annot = glue("rate<sub>max</sub> {sprintf('%0.2f', r_max)}")) |> 
  ggplot(aes(x = Time, y = value, group = well), color = "black") +
  geom_line() +
  geom_line(aes(x = Time, y = fit,  group = well), color = "red") +
  geom_vline(aes(xintercept = s), color = "blue") +
  geom_hline(aes(yintercept = p_max), color = "green") +
  geom_hline(aes(yintercept = p_min), color = "green") +
 # geom_richtext(aes(label = annot), x = 5, y = 4500, hjust = 1, vjust = 1) +
  facet_wrap(~well) +
  theme_classic()
  

cleaned |>
  inner_join(parameters) |> 
  arrange(well) |> 
  rowwise() |> 
  mutate(across(.col = c(value, p_max, p_min, r_max), .fns = ~./p_min),
         fit = richard2(Time, p_max, p_min, r_max, s)) |>
    mutate(condition = case_when(
    well %in% c("E1", "E2", "E3") ~ "A",
    well %in% c("E4", "E5", "E6") ~ "B",
    well %in% c("E7", "E8", "E9") ~ "C",
    well %in% c("F1", "F2", "F3") ~ "D",
    well %in% c("F4", "F5", "F6") ~ "E",
    well %in% c("F7", "F8", "F9") ~ "F",
    well %in% c("G1", "G2", "G3") ~ "G",
    well %in% c("G4", "G5", "G6") ~ "H",
    well %in% c("G7", "C8", "G9") ~ "I"
  )) |> 
  filter(condition %in% c("D", "E", "F")) |> 
  ggplot() +
  geom_line(aes(x = Time, y = value, group = well, color = condition)) +
  geom_line(aes(x = Time, y = fit,  group = well, color = condition),linetype = "dashed") +
#  geom_vline(aes(xintercept = s), color = "blue") +
#  geom_hline(aes(yintercept = p_max), color = "green") +
#  geom_hline(aes(yintercept = p_min), color = "green") +
 # geom_richtext(aes(label = annot), x = 5, y = 4500, hjust = 1, vjust = 1) +
  theme_classic()

cleaned |>
  inner_join(parameters) |> 
  arrange(well) |> 
  rowwise() |> 
  mutate(across(.col = c(value, p_max, p_min, r_max), .fns = ~./p_min),
         fit = richard2(Time, p_max, p_min, r_max, s)) |>
  mutate(condition = case_when(
    well %in% c("E1", "E2", "E3") ~ "A",
    well %in% c("E4", "E5", "E6") ~ "B",
    well %in% c("E7", "E8", "E9") ~ "C",
    well %in% c("F1", "F2", "F3") ~ "D",
    well %in% c("F4", "F5", "F6") ~ "E",
    well %in% c("F7", "F8", "F9") ~ "F",
    well %in% c("G1", "G2", "G3") ~ "G",
    well %in% c("G4", "G5", "G6") ~ "H",
    well %in% c("G7", "C8", "G9") ~ "I"
  )) |> 
  ggplot(aes(y = p_max, x = r_max, color = condition)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1)
  
y <- (p_max-p_min)/2
((log(((p_max-p_min) / ((p_max-p_min)/2-p_min)) - 1) * (p_min-p_max)) / (4*r_max)) + s

((log(-(-p_min+(3*p_max)) / (p_min+p_max)) * (p_min-p_max))/(4*r_max))+s

cleaned |>
  inner_join(parameters) |> 
  arrange(well) |> 
  rowwise() |> 
  mutate(across(.col = c(value, p_max, p_min, r_max), .fns = ~./p_min))  |> 
  mutate(condition = case_when(
      well %in% c("E1", "E2", "E3") ~ "A",
      well %in% c("E4", "E5", "E6") ~ "B",
      well %in% c("E7", "E8", "E9") ~ "C",
      well %in% c("F1", "F2", "F3") ~ "D",
      well %in% c("F4", "F5", "F6") ~ "E",
      well %in% c("F7", "F8", "F9") ~ "F",
      well %in% c("G1", "G2", "G3") ~ "G",
      well %in% c("G4", "G5", "G6") ~ "H",
      well %in% c("G7", "C8", "G9") ~ "I"
    )) |> 
  filter(condition %in% c("D", "E", "F")) |> 
  group_by(well) |> 
  summarise(across(.cols = where(is.numeric),.fns = mean)) |> 
  select(-c(well, Time, value, p_min)) |> 
  corrr::correlate() -> corr


corrr::rplot(corr)

$$
 t=\frac{\ln \left(\frac{A-Y}{Y-B}\right)\left(B-A\right)}{4r}+s \\
t=\frac{\ln \left(-\frac{-B+3A}{B+A}\right)\left(B-A\right)}{4r}+s
$$

{r}

x <- (98464.426-43614.260)
plyr::round_any(x, 10^(floor(log10(abs(x)))), f = ceiling)




