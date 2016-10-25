/*
######################################################################
#
# networkapi.c
#
# Written by Carter T. Butts <buttsc@uci.edu>
# Last Modified 5/07/16
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Provides support for the R/network package API.
#
# This file is known to be compatible with network package version 1.3.  
# It should be compatible with subsequent versions, but updates may be 
# necessary in rare cases.
#
# This file contains the registration routine needed to use the 
# C-level network package API.
#
######################################################################
*/

/*INCLUSIONS----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include "networkapi.h"


/*INTERNAL FUNCTIONS--------------------------------------------------------*/

void netRegisterFunctions(void)
/*Register functions for the network package API.  This function must be called before using any API routines, since these routines will not otherwise be defined within the local namespace.*/
{
  /*Register access routines*/
  netGetEdgeAttrib_ptr = (SEXP (*)(SEXP, int, const char*)) R_GetCCallable("network", "getEdgeAttribute");
  netGetEdgeIDs_ptr = (SEXP (*)(SEXP, int, int, const char*, int)) R_GetCCallable("network", "getEdgeIDs");
  netGetEdges_ptr = (SEXP (*)(SEXP, int, int, const char*, int)) R_GetCCallable("network", "getEdges");
  netGetNeighborhood_ptr = (SEXP (*)(SEXP, int, const char*, int)) R_GetCCallable("network", "getNeighborhood");
  netGetNetAttrib_ptr = (SEXP (*)(SEXP, const char*)) R_GetCCallable("network", "getNetworkAttribute");
  netHasLoops_ptr = (int (*)(SEXP)) R_GetCCallable("network", "hasLoops");
  netIsAdj_ptr = (int (*)(SEXP, int, int, int)) R_GetCCallable("network", "isAdjacent");
  netIsDir_ptr = (int (*)(SEXP)) R_GetCCallable("network", "isDirected");
  netIsHyper_ptr = (int (*)(SEXP)) R_GetCCallable("network", "isHyper");
  netIsLoop_ptr = (int (*)(SEXP, SEXP)) R_GetCCallable("network", "isLoop");
  netIsMulti_ptr = (int (*)(SEXP)) R_GetCCallable("network", "isMultiplex");
  netIsNetwork_ptr = (int (*)(SEXP)) R_GetCCallable("network", "isNetwork");
  netNetEdgecount_ptr = (int (*)(SEXP, int)) R_GetCCallable("network", "networkEdgecount");
  netNetSize_ptr = (int (*)(SEXP)) R_GetCCallable("network", "networkSize");

  //netGetEdgeAttrib_ptr = (SEXP (*)(SEXP, int, const char*)) R_FindSymbol("getEdgeAttribute", "network", NULL);
  //netGetEdgeIDs_ptr = (SEXP (*)(SEXP, int, int, const char*, int))R_FindSymbol("getEdgeIDs", "network", NULL);
  //netGetEdges_ptr = (SEXP (*)(SEXP, int, int, const char*, int))R_FindSymbol("getEdges", "network", NULL);
  //netGetNeighborhood_ptr = (SEXP (*)(SEXP, int, const char*, int))R_FindSymbol("getNeighborhood", "network", NULL);
  //netGetNetAttrib_ptr = (SEXP (*)(SEXP, const char*))R_FindSymbol("getNetworkAttribute", "network", NULL);
  //netHasLoops_ptr = (int (*)(SEXP))R_FindSymbol("hasLoops", "network", NULL);
  //netIsAdj_ptr = (int (*)(SEXP, int, int, int))R_FindSymbol("isAdjacent", "network", NULL);
  //netIsDir_ptr = (int (*)(SEXP))R_FindSymbol("isDirected", "network", NULL);
  //netIsHyper_ptr = (int (*)(SEXP))R_FindSymbol("isHyper", "network", NULL);
  //netIsLoop_ptr = (int (*)(SEXP, SEXP))R_FindSymbol("isLoop", "network", NULL);
  //netIsMulti_ptr = (int (*)(SEXP))R_FindSymbol("isMultiplex", "network", NULL);
  //netIsNetwork_ptr = (int (*)(SEXP))R_FindSymbol("isNetwork", "network", NULL);
  //netNetEdgecount_ptr = (int (*)(SEXP, int))R_FindSymbol("networkEdgecount", "network", NULL);
  //netNetSize_ptr = (int (*)(SEXP))R_FindSymbol("networkSize", "network", NULL);

  /*Register modification routines*/
  netAddEdge_ptr = (SEXP (*)(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP)) R_GetCCallable("network", "addEdge_R");
  netAddEdges_ptr = (SEXP (*)(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP)) R_GetCCallable("network", "addEdges_R");
  netDelEdgeAttrib_ptr = (SEXP (*)(SEXP, int, const char*)) R_GetCCallable("network", "deleteEdgeAttribute");
  netDelVertexAttrib_ptr = (SEXP (*)(SEXP, int, const char*)) R_GetCCallable("network", "deleteVertexAttribute");
  netDelNetAttrib_ptr = (SEXP (*)(SEXP, const char*)) R_GetCCallable("network", "deleteNetworkAttribute");
  netSetNetAttrib_ptr = (SEXP (*)(SEXP, const char*, SEXP)) R_GetCCallable("network", "setNetworkAttribute");
  netSetVertexAttrib_ptr = (SEXP (*)(SEXP, const char*, SEXP, int)) R_GetCCallable("network", "setVertexAttribute");

  //netAddEdge_ptr = (SEXP (*)(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP))R_FindSymbol("addEdge_R", "network", NULL);
  //netAddEdges_ptr = (SEXP (*)(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP))R_FindSymbol("addEdges_R", "network", NULL);
  //netDelEdgeAttrib_ptr = (SEXP (*)(SEXP, int, const char*))R_FindSymbol("deleteEdgeAttribute", "network", NULL);
  //netDelNetAttrib_ptr = (SEXP (*)(SEXP, const char*))R_FindSymbol("deleteNetworkAttribute", "network", NULL);
  //netDelVertexAttrib_ptr = (SEXP (*)(SEXP, int, const char*))R_FindSymbol("deleteVertexAttribute", "network", NULL);
  //netSetNetAttrib_ptr = (SEXP (*)(SEXP, const char*, SEXP))R_FindSymbol("setNetworkAttribute", "network", NULL);
  //netSetVertexAttrib_ptr = (SEXP (*)(SEXP, const char*, SEXP, int))R_FindSymbol("setVertexAttribute", "network", NULL);
}

