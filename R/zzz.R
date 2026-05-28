.onLoad <- function(libname, pkgname) {
  tmb_dir <- system.file("tmb", package = pkgname)

  for (model in c("fn_BC_tmb", "fn_log_tmb")) {
    cpp_file <- file.path(tmb_dir, paste0(model, ".cpp"))
    dll_file  <- file.path(tmb_dir, paste0(model, .Platform$dynlib.ext))

    # compile only if DLL doesn't exist or is older than the source
    if (!file.exists(dll_file) ||
        file.mtime(dll_file) < file.mtime(cpp_file)) {
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
