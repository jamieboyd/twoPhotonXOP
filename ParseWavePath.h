/* Last Modified 2014/09/23 by Jamie Boyd 
   added conditional typedef for UInt16 */

#ifndef PARSEWAVEPATH_H
#define PARSEWAVEPATH_H

#include "XOPStandardHeaders.h"

// typedefs for paramater string handling
typedef char DFPATH [MAXCMDLEN + 1];
typedef char WVNAME [MAX_OBJ_NAME + 1];

/* Prototypes */
void ParseWavePath (Handle fullPath, DFPATH dataFolderName, WVNAME waveName);
waveHndl* ParseWaveListPaths (Handle pathsList, UInt16 *nWavesPtr);

//#define IPARSEWAVE_INF0         //to print out wavenames and paths for checking results of ParseWavePath
 #undef IPARSEWAVE_INF0

#endif
