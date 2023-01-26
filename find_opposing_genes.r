library(tidyverse)

genome <- read_tsv("/home/thomas/BioinformaticsTools/Annotation_out/Annotated/mybb4.tsv",   comment = "#")


is_facing <- function(strand) {
    return(lag({{strand}}) == "+" & {{strand}} == "-")
}

min_distance <- 50
genome |>
    mutate(is_facing = if_else(is_facing(Strand), TRUE, FALSE),
            distance_between = Start - lag(Stop)) |>
    filter(is_facing, distance_between >= min_distance) |>
    arrange(desc(distance_between)) -> filtered

filtered |> write_tsv("~/Bureau/mybb4_opposing.tsv")
