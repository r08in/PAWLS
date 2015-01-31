#This is the function for prepare matlab function
#
PrepareMatlab=function(dir,matlab=NULL)
{
  require("R.matlab")
  if(is.null(matlab)||!isOpen(matlab))
  {
    Matlab$startServer()
    matlab <- Matlab()
    open(matlab)
  }
  nameList=list.files(path=dir)
  fileList=paste(dir,nameList,sep="\\")
  n=length(fileList)
  for(i in 1:n)
  {
    str=readChar(fileList[i], file.info(fileList[i])$size)
    setFunction(matlab,str)
  }
  matlab
}