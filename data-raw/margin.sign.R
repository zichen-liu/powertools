## code to prepare `delta.sign` dataset goes here
margin.sign <- data.frame(" " = c("Noninferiority", "Superiority by margin"),
                 "Higher is better" = c("-", "+"),
                 "Lower is better" = c("+", "-"),
                 check.names = F)

usethis::use_data(margin.sign, overwrite = TRUE)
