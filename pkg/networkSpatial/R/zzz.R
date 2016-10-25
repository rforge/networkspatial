######################################################################
#
# zzz.R
#
# Written by Carter T. Butts <buttsc@uci.edu>.
#
# Last Modified 05/22/09
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/networkSpatial package
#
# .First.lib is run when the package is loaded with library(networkSpatial)
#
######################################################################


#.First.lib <- function(lib, pkg){
#   library.dynam("networkSpatial", pkg, lib)
#   require(network, quietly=TRUE)
#    if(R.version$major=="1"){
#     ehelp <- help(package="networkSpatial")$info[[2]][[2]]
#     cat(paste("\n'",ehelp[4],"'\n",
#               "Version ",ehelp[2],
#               " created on ",ehelp[3],".\n", sep=""))
#    }else{
#     ehelp <- help(package="networkSpatial")$info[[1]]
#     cat(paste("\n",substring(ehelp[4],first=16),"\n",
#               "Version ",substring(ehelp[2],first=16),
#               " created on ",
#                substring(ehelp[3],first=16),".\n", sep=""))
#    }
#    cat(paste("copyright (c) 2009, Carter T. Butts, University of California-Irvine\n",sep=""))
#    cat('For citation information, type citation("networkSpatial").\n')
#    cat('Type help(package="networkSpatial") to get started.\n')
#}

#.onLoad <- function(libname, pkgname){
#  library.dynam("networkSpatial", package=pkgname, lib.loc=libname)
#}

.onAttach <- function(libname, pkgname){
  	temp<-packageDescription("networkSpatial")
  	msg<-paste(temp$Package,": ",temp$Title,"\n",
      "Version ",temp$Version,
      " created on ",
      temp$Date,".\n", sep="")
  msg<-paste(msg,"copyright (c) 2013, Carter T. Butts, University of California-Irvine\n",sep="")
  msg<-paste(msg,"                    Zack W. Almquist, University of Minnesota\n",sep="")
  msg<-paste(msg,'For citation information, type citation("networkSpatial").\n')
  msg<-paste(msg,'Type help(package="networkSpatial") to get started.\n')
  packageStartupMessage(msg)
}

