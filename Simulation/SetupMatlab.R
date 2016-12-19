# This is the function for prepare matlab function
PrepareMatlab = function(dir, matlab = NULL) {
    oldwd <- getwd()
    setwd(dir)
    require("R.matlab")
    if (is.null(matlab) || !isOpen(matlab)) {
        Matlab$startServer()
        matlab <- Matlab()
        open(matlab)
    }
    nameList = list.files(path = "matlab")
    fileList <- paste("matlab",nameList , sep = "/")
    n = length(fileList)
    pb <- txtProgressBar(1, n, style = 3)
    for (i in 1:n) {
        str = readChar(fileList[i], file.info(fileList[i])$size)
        setFunction(matlab, str)
        setTxtProgressBar(pb, i)
    }
    setwd(oldwd)
    matlab
}
