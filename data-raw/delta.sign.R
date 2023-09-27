## code to prepare `delta.sign` dataset goes here
delta.sign <- data.frame(" " = c("Noninferiority", "Superiority by margin"),
                 "Higher is better" = c("-", "+"),
                 "Lower is better" = c("+", "-"),
                 check.names = F)

usethis::use_data(delta.sign, overwrite = TRUE)
