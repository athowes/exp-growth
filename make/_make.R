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
