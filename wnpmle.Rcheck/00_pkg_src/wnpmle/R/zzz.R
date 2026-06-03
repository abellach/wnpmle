.onLoad <- function(libname, pkgname) {
  tmb_dir <- system.file("tmb", package = pkgname)

  for (model in c("fn_BC_tmb", "fn_log_tmb")) {
    dll_file <- file.path(tmb_dir, paste0(model, .Platform$dynlib.ext))

    # only compile if DLL is missing (should be pre-compiled in package)
    if (!file.exists(dll_file)) {
      cpp_file <- file.path(tmb_dir, paste0(model, ".cpp"))
      TMB::compile(cpp_file, silent = TRUE)
    }

    dyn.load(TMB::dynlib(file.path(tmb_dir, model)))
  }
}

.onUnload <- function(libpath) {
  tmb_dir <- system.file("tmb", package = "wnpmle")
  for (model in c("fn_BC_tmb", "fn_log_tmb")) {
    dll <- TMB::dynlib(file.path(tmb_dir, model))
    if (is.loaded(model)) dyn.unload(dll)
  }
}
