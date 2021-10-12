if(auxiliary_dist == "Rademacher"){

  wts <- sample(c(-1, 1),
                size = length(num_cluster),
                replace = TRUE,
                prob = rep(1/2, 2))

} else if(auxiliary_dist == "Mammen"){

 wts <- sample(c(- (sqrt(5) - 1)/2,
                 (sqrt(5) + 1)/2),
               size = length(num_cluster),
               replace = TRUE,
               prob = c((sqrt(5) + 1) /(2 * sqrt(5)), (sqrt(5) - 1) /(2 * sqrt(5)) ))

} else if(auxiliary_dist == "Webb six"){

 wts <- sample(c(-sqrt(3/2),
                 -sqrt(2/2),
                 -sqrt(1/2),
                 sqrt(1/2),
                 sqrt(2/2),
                 sqrt(3/2)),
               size = length(num_cluster),
               replace = TRUE,
               prob = rep(1/6, 6))

} else if(auxiliary_dist = "uniform"){

  wts <- runif(n = length(num_cluster),
               min = -sqrt(3),
               max = sqrt(3))

} else if(auxiliary_dist = "standard normal"){

  wts <- rnorm(n = length(num_cluster))

}

