.onAttach <- function(libname, pkgname) {
  packageStartupMessage(pkgname, ' v', utils::packageVersion(pkgname),
                        ': OpenMP ', ifelse(chk_omp(), 'enabled', 'disabled'))
  register_autoplot_s3_methods()
}
if(getRversion() >= "2.15.1") {
  utils::globalVariables(c(".data","level"))
}
