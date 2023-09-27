## code to prepare `delta.sign` dataset goes here
delta <- data.frame(" " = c("Noninferiority", "Superiority by a margin"),
                 "Higher is better" = c("-", "+"),
                 "Lower is better" = c("+", "-"),
                 check.names = F)

usethis::use_data(delta, overwrite = TRUE)
