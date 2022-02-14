
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Attaching expTSNE version ",
                        packageDescription("expTSNE")$Version, ".")
  #When adding new options here, also add them to the "names" setMethod below
  options("ET_FORCE_CACHE_OVERWRITE" = FALSE)
  options("ET_CACHE_VERSION" = "v1")
  options("ET_CACHE_PATH" = "~/.cache")
  ET_OPTIONS <<- new("ET_OPTIONS")
}

setClass("ET_OPTIONS", representation = list(
  is_valid = "logical"
  ))

setMethod("names", "ET_OPTIONS",
          function(x)
          {
            c(
              "ET_FORCE_CACHE_OVERWRITE",
              "ET_CACHE_VERSION",
              "ET_CACHE_PATH",
              "mc.cores"
            )
          })


setMethod("$", "ET_OPTIONS",
          function(x, name)
          {
            getOption(name)
          })

setReplaceMethod("$", "ET_OPTIONS",
                 function(x, name, value)
                 {
                   warn_msg = "This assignment is not supported.  No effect."
                   value = list(value)
                   names(value) = name
                   do.call("options", value)
                   x
                 })



setMethod("show", "ET_OPTIONS",
          function(object)
          {
            message("Use the $ accessor (i.e. ET_OPTIONS$ET_CACHE_PATH) to get/set ET relevant options.")
            message("Use names(ET_OPTIONS) to view all options.")
          })
