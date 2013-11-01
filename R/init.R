.setCAMERAOptions <- function(pkgname, camera.opt=NA) {

  if (! any(is.na(camera.opt))) {
    if (class(camera.opt) != "BioCPkg")
      stop("obviously invalid package options !")

    BioC <- getOption("BioC")
    BioC$CAMERA <- camera.opt
    options("BioC"=BioC)
    return()
  }

  ## add xcms specific options
  ## (not unlike what is done in 'affy')
  if (is.null(getOption("BioC"))) {
    BioC <- list()
    class(BioC) <- "BioCOptions"
    options("BioC"=BioC)
  }

  ## all calcPC methods
  start <- nchar("calcPC.")
  all.camera <- ls(asNamespace(pkgname))
  calcPC.methods <-  substr(all.camera[grep("calcPC\\..*", all.camera)], start+1, 100)

  ## default for the methods
  calcPC.method <- "hcs"

  ## all combineCalc methods
  start <- nchar("combineCalc.")
  all.camera <- ls(asNamespace(pkgname))
  combineCalc.methods <-  substr(all.camera[grep("combineCalc\\..*", all.camera)], start+1, 100)

  ## default for the methods
  combineCalc.method <- "sum"
 
  camera.opt <- list(calcPC.methods=calcPC.methods, calcPC.method=calcPC.method,combineCalc.methods=combineCalc.methods,combineCalc.method = combineCalc.method)
                   
  class(camera.opt) <- "BioCPkg"

  BioC <- getOption("BioC")
  BioC$CAMERA <- camera.opt
  options("BioC"=BioC)
}
