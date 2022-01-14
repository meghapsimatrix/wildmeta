library(clubSandwich)

full_dat <- readRDS("scratchpad/TD_DLD_full_data.rds")
r <- 0.4
Vmat <- impute_covariance_matrix(vi = full_dat$v, cluster = full_dat$sampleid, r = r)

# RQ1

RQ1_A <-
  rma.mv(yi = d, V = Vmat,
         random = ~ 1 | studyid / sampleid / esid,
         mods = ~ 0 + cat_ + DLD_age_yrs,
         data = full_dat, sparse = TRUE
  )

Cmat_1 <- constrain_equal("cat_", reg_ex = TRUE)
Wald_test(RQ1_A, constraints = Cmat_1, vcov = "CR2")
Wald_test_cwb(RQ1_A, constraints = Cmat_1, R = 99)

RQ1_B <-
  rma.mv(yi = d, V = Vmat,
         random = ~ 1 | studyid / sampleid / esid,
         mods = ~ 0 + cat_ + DLD_age_yrs + nonEnglish_task + bilingual + L2,
         data = full_dat, sparse = TRUE
  )
Wald_test(RQ1_B, constraints = Cmat_1, vcov = "CR2")
Wald_test_cwb(RQ1_B, constraints = Cmat_1, R = 99)


# RQ2

RQ2_A <-
  rma.mv(yi = d, V = Vmat,
         random = ~ 1 | studyid / sampleid / esid,
         mods = ~ 0 + con_ + DLD_age_yrs,
         data = full_dat, subset = cat_ == "macro" & !is.na(con_),
         sparse = TRUE
  )
Cmat_2 <- constrain_equal(2:4)
Wald_test(RQ2_A, constraints = Cmat_2, vcov = "CR2")
Wald_test_cwb(RQ2_A, constraints = Cmat_2, R = 199)

RQ2_B <-
  rma.mv(yi = d, V = Vmat,
         random = ~ 1 | studyid / sampleid / esid,
         mods = ~ 0 + con_ + DLD_age_yrs + nonEnglish_task + bilingual + L2,
         data = full_dat, subset = cat_ == "macro" & !is.na(con_),
         sparse = TRUE
  )
Wald_test(RQ2_B, constraints = Cmat_2, vcov = "CR2")
Wald_test_cwb(RQ2_B, constraints = Cmat_2, R = 99)

# RQ3

RQ3_A <-
  rma.mv(yi = d, V = Vmat,
         random = ~ 1 | studyid / sampleid / esid,
         mods = ~ 0 + con_ + DLD_age_yrs,
         data = full_dat, subset = cat_ == "micro" & !is.na(con_),
         sparse = TRUE
  )
Cmat_3 <- constrain_equal("con_", reg_ex = TRUE)
Wald_test(RQ3_A, constraints = Cmat_3, vcov = "CR2")
Wald_test_cwb(RQ3_A, constraints = Cmat_3, R = 99)

RQ3_B <-
  rma.mv(yi = d, V = Vmat,
         random = ~ 1 | studyid / sampleid / esid,
         mods = ~ 0 + con_ + DLD_age_yrs + nonEnglish_task + bilingual + L2,
         data = full_dat, subset = cat_ == "micro" & !is.na(con_),
         sparse = TRUE
  )
Wald_test(RQ3_B, constraints = Cmat_3, vcov = "CR2")
Wald_test_cwb(RQ3_B, constraints = Cmat_3, R = 99)

# RQ4

RQ4_A <-
  rma.mv(yi = d, V = Vmat,
         random = ~ 1 | studyid / sampleid / esid,
         mods = ~ 0 + cat_ + recruit_ + inclusion_ + DLD_age_yrs,
         data = full_dat, subset = !is.na(inclusion_),
         sparse = TRUE
  )
Cmat_4A <- constrain_zero("recruit_",reg_ex = TRUE)
Cmat_4B <- constrain_zero("inclusion_", reg_ex = TRUE)
Wald_test(RQ4_A,
          constraints = list(recruit = Cmat_4A, inclusion = Cmat_4B),
          vcov = "CR2")
Wald_test_cwb(RQ4_A, constraints = Cmat_4A, R = 99)
Wald_test_cwb(RQ4_A, constraints = Cmat_4B, R = 99)

RQ4_B <-
  rma.mv(yi = d, V = Vmat,
         random = ~ 1 | studyid / sampleid / esid,
         mods = ~ 0 + cat_ + recruit_ + inclusion_ + DLD_age_yrs + nonEnglish_task + bilingual + L2,
         data = full_dat, subset = !is.na(inclusion_),
         sparse = TRUE
  )
Wald_test(RQ4_B,
          constraints = list(recruit = Cmat_4A, inclusion = Cmat_4B),
          vcov = "CR2")
Wald_test_cwb(RQ4_B, constraints = Cmat_4A, R = 99)
Wald_test_cwb(RQ4_B, constraints = Cmat_4B, R = 99)


# RQ5

RQ5_A <-
  rma.mv(yi = d, V = Vmat,
         random = ~ 1 | studyid / sampleid / esid,
         mods = ~ 0 + cat_ + task_ + support_ + DLD_age_yrs,
         data = full_dat,
         sparse = TRUE
  )
Cmat_5A <- constrain_zero("task_",reg_ex = TRUE)
Cmat_5B <- constrain_zero("support_", reg_ex = TRUE)
Wald_test(RQ5_A,
          constraints = list(task = Cmat_5A, support = Cmat_5B),
          vcov = "CR2")
Wald_test_cwb(RQ5_A, constraints = Cmat_5A, R = 199)
Wald_test_cwb(RQ5_A, constraints = Cmat_5B, R = 199)

RQ5_B <-
  rma.mv(yi = d, V = Vmat,
         random = ~ 1 | studyid / sampleid / esid,
         mods = ~ 0 + cat_ + task_ + support_ +
           DLD_age_yrs + nonEnglish_task + bilingual + L2,
         data = full_dat,
         sparse = TRUE
  )
Wald_test(RQ5_B,
          constraints = list(task = Cmat_5A, support = Cmat_5B),
          vcov = "CR2")
Wald_test_cwb(RQ5_B, constraints = Cmat_5A, R = 199)
Wald_test_cwb(RQ5_B, constraints = Cmat_5B, R = 199)


# RQ6

RQ6_A <-
  rma.mv(yi = d, V = Vmat,
         random = ~ 1 | studyid / sampleid / esid,
         mods = ~ 0 + index_ + recruit_ + inclusion_ + DLD_age_yrs,
         data = full_dat, subset = include_index == "X",
         sparse = TRUE
  )
Cmat_6 <- constrain_equal("index_", reg_ex = TRUE)
Wald_test(RQ6_A, constraints = Cmat_6, vcov = "CR2")
Wald_test_cwb(RQ6_A, constraints = Cmat_6, R = 99)

RQ6_B_CHEp <-
  rma.mv(yi = d, V = Vmat,
         random = ~ 1 | studyid / sampleid / esid,
         mods = ~ 0 + index_ + recruit_ + inclusion_ + DLD_age_yrs + nonEnglish_task + bilingual + L2,
         data = full_dat, subset = include_index == "X",
         sparse = TRUE
  )
Wald_test(RQ6_B, constraints = Cmat_6, vcov = "CR2")
Wald_test_cwb(RQ6_B, constraints = Cmat_6, R = 99)
