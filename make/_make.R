id <- orderly::orderly_run("explore_huisman")
orderly::orderly_commit(id)

id <- orderly::orderly_run("explore_qPCR")
orderly::orderly_commit(id)

id <- orderly::orderly_run("fit_qPCR")
orderly::orderly_commit(id)

id <- orderly::orderly_run("explore_poisson-regression")
orderly::orderly_commit(id)

id <- orderly::orderly_run("docs_CoDA")
orderly::orderly_commit(id)

id <- orderly::orderly_run("explore_ppp-calc")
orderly::orderly_commit(id)

id <- orderly::orderly_run("sim_metagenomic-time-series")
orderly::orderly_commit(id)

id <- orderly::orderly_run("fit_poisson-regression")
orderly::orderly_commit(id)

id <- orderly::orderly_run("benchmark_poisson-regression")
orderly::orderly_commit(id)

id <- orderly::orderly_run("explore_compartmental")
orderly::orderly_commit(id)

id <- orderly::orderly_run("sim_airplane-time-series")
orderly::orderly_commit(id)
