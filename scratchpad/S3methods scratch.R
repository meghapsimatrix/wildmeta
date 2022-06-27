my_method <- function(x, y) UseMethod("my_method")

my_method.integer <- function(x, y) {
  sum(x) + y
}

my_method.numeric <- function(x, y) {
  prod(x) * y
}

my_method.character <- function(x, y) {
  paste(paste(x, collapse = ""), y)
}

my_method(1:4, 40)
my_method(runif(4), 40)
my_method(LETTERS, 40)
my_method(NULL, 40)

lapply(101:103, my_method, x = 1:4)
lapply(101:103, my_method, x = runif(4))
lapply(101:103, my_method, x = LETTERS)

library(future.apply)
future_lapply(101:103, my_method, x = 1:4)
future_lapply(101:103, my_method, x = runif(4))
future_lapply(101:103, my_method, x = LETTERS)
