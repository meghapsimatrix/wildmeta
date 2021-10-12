if(auxiliary_dist == "Rademacher"){

  wts <- sample(c(-1, 1), size = length(num_cluster), replace = TRUE)

} else if(auxiliary_dist == "Mammen"){

 wts <- sample(c(- (sqrt(5) - 1)/2,  (sqrt(5) + 1)/2), size = length(num_cluster), replace = TRUE,
               prob = c((sqrt(5) + 1) /(2 * sqrt(5)), (sqrt(5) - 1) /(2 * sqrt(5)) ))

}

