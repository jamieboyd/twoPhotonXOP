#include "twoPhoton.h"

/* ------------------------------Minimum/Maximum/Avg/Median Intensity Projections----------------------
 Code for making Minimum/Maximum/Avg/Median Intensity Projections along X,Y, and Z axes
 Last Modified 2025/06/25 by Jamie Boyd
 
 ProjectAllFrames and ProjectSpecFrames make a projection image of a 3D wave, either into a 2D wave,
 or by collapsing the input 3D wave into 2D. The differences between ProjectAllFrames and ProjectSpecFrames:
  1 ProjectAllFrames makes a project image of ALL the layers in a 3D wave, while ProjectSpecFrames does a projection over a
  specified range of frames, with a starting and ending element
  2 ProjectAllFrames has a path to an output wave as parameter, and will create the needed wave if it doesn't exist
  3 ProjectAllFrames output is always to a 2D wave, flattening the output wave to 2D if it has 3D (if input is same as output, e.g) ProjectSpecFrames can output a projection to a given X-Y plane in a 3D wave
  4 ProjectAllFrames has a flag for overwriting output wave if it already exists. ProjectSpecFrames does not

 ProjectXslice, ProjectYSlice, and ProjectZSLice are used to get a single slice from a 3D wave and place it in an
 existing 2D wave of the right dimensions
 ------------------------------------------------------------------------------------------------------- */
 
/* The following templates get a single slice or makes a minimum or maximum intensity projection for one of
 the 8 types of wave data in the X,Y,or Z dimension, for either of the Project functions
 srcWaveStart: pointer to the start of the data in the 3D input wave
 desWaveStart: pointer to the start of the data in the 2D or 3D outPut wave
 projMode: 1 for maximum intensity projection, 0 for minimum intensity projection
 example:
 X projection - looking for max/min with a range of columns at a single YZ location
 result is wave with dim 0 = ySize, dim 1 =zSize
 using Igor terminology of x being columns and y being rows and z being layers
 Last Modified 2014/01/27 by Jamie Boyd */

template <typename T> void doProjectX(T *srcWaveStart, T *destWaveStart, UInt8 projMode, CountInt xSize, CountInt ySize, CountInt zSize, CountInt startP, CountInt endP, UInt8 ti, UInt8 tN){
    // number of YZ locations to look at
    CountInt yzLocs = ySize * zSize;
    // Number of ZY locations to do for this thread. divide YZ locations among threads. truncated to an integer
    CountInt tPoints = yzLocs/tN;
    // which point to start this thread on depends on thread number
    CountInt inStartPos = startP + (ti * tPoints * xSize);
    CountInt outStartPos = ti * tPoints;
    // the last thread gets any left-over YZ locations
    if (ti == tN - 1) tPoints += (yzLocs % tN);
    // number of columns to do within each YZ location.
    CountInt ColumnsToDo = endP - startP + 1;
    // number to skip to get to start of next column at this ZY location
    // number of points to get to next YZloc
    CountInt toNextYZ = xSize - ColumnsToDo;
    // pointers for loops
    T *srcWave = srcWaveStart + inStartPos;
    T *destWave = destWaveStart + outStartPos;
    T *destWaveEnd;
    T* lastColumn; // last column at this YZ location
    if (ColumnsToDo == 1){ // getting a single slice
        // adjust toNextYZ for lack of inner loop
        toNextYZ += 1;
        // loop through YZ locations
        for (destWaveEnd = destWave + tPoints; destWave < destWaveEnd; destWave++, srcWave += toNextYZ){
            // set destination to chosen column at this source ZY location
            *destWave = *srcWave;
        }
    }else{ // doing a projection
        switch(projMode){
            case 0: // Minimum projection
                // loop through YZ locations looking for minimum
                for (destWaveEnd = destWave + tPoints; destWave < destWaveEnd; destWave++, srcWave += toNextYZ){
                    // set destination to first column at this source ZY location
                    *destWave = *srcWave;
                    // check remining columns in source to find max
                    for (lastColumn= srcWave + ColumnsToDo, srcWave++ ; srcWave < lastColumn; srcWave++){
                        if (*srcWave < *destWave)
                            *destWave = *srcWave;
                    }
                }
                break;
            case 1: // Maximum projection
                // loop through YZ locations looking for maximum
                for (destWaveEnd = destWave + tPoints; destWave < destWaveEnd; destWave++, srcWave += toNextYZ){
                    // set destination to first column at this source ZY location
                    *destWave = *srcWave;
                    // check remining columns in source to find min
                    for (lastColumn= srcWave + ColumnsToDo, srcWave++ ; srcWave < lastColumn; srcWave++){
                        if (*srcWave > *destWave)
                            *destWave = *srcWave;
                    }
                }
                break;
            case 2: // Average projection
                double temp;
                // loop through YZ locations calculating average, with double precision float
                for (destWaveEnd = destWave + tPoints; destWave < destWaveEnd; destWave++, srcWave += toNextYZ){
                    for (lastColumn= srcWave + ColumnsToDo, temp=0; srcWave < lastColumn; srcWave++){
                        temp += *srcWave;
                    }
                    *destWave = temp/ColumnsToDo;
                }
                break;
            case 3: // Median projection
                T *bufferStart = (T*) WMNewPtr (ColumnsToDo * sizeof(T));
                T *bufferPos;
                // loop through YZ locations filling a buffer for calculating median
                for (destWaveEnd = destWave + tPoints; destWave < destWaveEnd; destWave++, srcWave += toNextYZ){
                    for (lastColumn= srcWave + ColumnsToDo, bufferPos =bufferStart; srcWave < lastColumn; srcWave++, bufferPos++){
                        *bufferPos = *srcWave;
                    }
                    *destWave = medianT (ColumnsToDo, bufferStart);
                }
                WMDisposePtr ((Ptr)bufferStart);
                break;
        }
    }
}

/* Y projection - looking for max/min within a range of rows at each XZ location
 result is wave with dim 0 = xSize, dim 1 =zSize
 Last Modified 2014/01/27 by Jamie Boyd */
template <typename T> void doProjectY(T *srcWaveStart, T *destWaveStart, UInt8 projMode, CountInt xSize, CountInt ySize, CountInt zSize, CountInt startP, CountInt endP, UInt8 ti, UInt8 tN){
    // number of layers to do in each thread = layers/threads.  truncated to an integer
    // each layer needs to be in a single thread because offset to next XZ location is different at end of a layer
    CountInt tLayers = zSize/tN;
    // which point to start this thread on depends on thread number * layers Per Thread . ti is 0 based
    CountInt inStartPos = (startP * xSize) + (ti * tLayers * xSize * ySize);
    CountInt outStartPos = ti * tLayers * xSize;
    // the last thread gets any left-over layers
    if (ti == tN - 1) tLayers += (zSize % tN);
    // number of rows to look at in each XZ location
    CountInt rowsToDo = endP - startP + 1;
    CountInt toNextRow = xSize; // number to add to get to next row at this YZ location
    CountInt toNextColumn = 1 - (rowsToDo * xSize); // number to add to get to next YZ location in this column
    CountInt toNextLayer = (xSize * ySize) - xSize; // number to add to get to first YZ location in the next layer
    T *srcWave = srcWaveStart + inStartPos;
    T *destWave = destWaveStart + outStartPos;
    T* lastRow; // last row at this XZ location
    T* lastColumn; //last row of last column in this layer
    T* lastLayer; // last row of last column of last layer
    if (rowsToDo == 1){ // getting a single slice
        // adjust toNextColumn for lack of innermost loop
        toNextColumn += toNextRow;
        // loop through layers (z)
        for (lastLayer = srcWave + (tLayers * xSize * ySize); srcWave < lastLayer; srcWave += toNextLayer){
            // loop through columns (x) in this layer
            for (lastColumn = srcWave + xSize; srcWave < lastColumn; srcWave += toNextColumn, destWave++){
                // set destination to chosen column in source row
                *destWave = *srcWave;
            }
        }
    }else{ // doing a projection
        switch (projMode) {
            case 0:  // minimum projection
                for (lastLayer = srcWave + (tLayers * xSize * ySize); srcWave < lastLayer; srcWave += toNextLayer){
                    // loop through columns (x) in this layer
                    for (lastColumn = srcWave + xSize; srcWave < lastColumn; srcWave += toNextColumn, destWave++){
                        // set destination to first row in source column
                        *destWave = *srcWave;
                        // loop through remaining rows (y) in this column (x)
                        for (lastRow = srcWave + (rowsToDo * xSize), srcWave += toNextRow;srcWave < lastRow; srcWave += toNextRow){
                            if (*srcWave < *destWave) *destWave = *srcWave;
                        }
                    }
                }
                break;
            case 1:  // max projection
                for (lastLayer = srcWave + (tLayers * xSize * ySize); srcWave < lastLayer; srcWave += toNextLayer){
                    // loop through columns (x) in this layer
                    for (lastColumn = srcWave + xSize; srcWave < lastColumn; srcWave += toNextColumn, destWave++){
                        // set destination to first row in source column
                        *destWave = *srcWave;
                        // loop through remaining rows (y) in this column (x)
                        for (lastRow = srcWave + (rowsToDo * xSize), srcWave += toNextRow;srcWave < lastRow; srcWave += toNextRow){
                            if (*srcWave > *destWave) *destWave = *srcWave;
                        }
                    }
                }
                break;
            case 2: // Avg projection
                double temp;
                for (lastLayer = srcWave + (tLayers * xSize * ySize); srcWave < lastLayer; srcWave += toNextLayer){
                    // loop through columns (x) in this layer
                    for (lastColumn = srcWave + xSize; srcWave < lastColumn; srcWave += toNextColumn, destWave++){
                        // loop through rows (y) in this column (x)
                        for (lastRow = srcWave + (rowsToDo * xSize), temp=0;srcWave < lastRow; srcWave += toNextRow){
                            temp += *srcWave;
                        }
                        *destWave = temp/rowsToDo;
                    }
                }
                break;
            case 3: // median proj
                T *bufferStart = (T*)WMNewPtr(rowsToDo * sizeof(T));
                T *bufferPos;
                for (lastLayer = srcWave + (tLayers * xSize * ySize); srcWave < lastLayer; srcWave += toNextLayer){
                    // loop through columns (x) in this layer
                    for (lastColumn = srcWave + xSize; srcWave < lastColumn; srcWave += toNextColumn, destWave++){
                        // loop through rows (y) in this column (x)
                        for (lastRow = srcWave + (rowsToDo * xSize), bufferPos=bufferStart;srcWave < lastRow; srcWave += toNextRow, bufferPos++){
                            *bufferPos = *srcWave;
                        }
                        *destWave = medianT (rowsToDo, bufferStart);
                    }
                }
                WMDisposePtr ((Ptr)bufferStart);
                break;
        }
        
    }
}

/* Z projection - looking for max/min in range of layers at same XY location
 result is wave with dim 0 = xSize, dim 1 =ySize
 Last Modified 2013/01/27 by Jamie Boyd  */
template <typename T> void doProjectZ(T *srcWaveStart, T *destWaveStart, UInt8 projMode, CountInt xSize, CountInt ySize, CountInt startP, CountInt endP, UInt8 ti, UInt8 tN){
    // total number of XY locations to process
    CountInt xyLocs = (xSize * ySize);
    // divide XY locations among threads. XYlocs per thread = total number of XYlocs/number of threads, truncated to an integer
    CountInt tPoints = xyLocs/tN;
    // which point to start this thread on depends on thread number * XY locations Per Thread . ti is 0 based
    CountInt inStartPos = (startP * xSize * ySize) + (ti * tPoints);
    CountInt outStartPos = ti * tPoints;
    // the last thread gets any left-over XY locations
    if (ti == tN - 1) tPoints += (xyLocs % tN);
    // number of layers to do within each XYloc, and number to skip to get to start of next layer at this XYlocs
    CountInt layersToDo = endP - startP + 1;
    CountInt toNextLayer = xSize * ySize;
    // number of points to get to next XYloc
    CountInt toNextXY = 1 - (layersToDo * xSize * ySize);
    // pointers for loops
    T *srcWave = srcWaveStart + inStartPos;
    T *destWave = destWaveStart + outStartPos;
    T *destWaveEnd;
    T* lastLayer; // last layer at this XY location
    if (layersToDo == 1){ // getting a single slice
        //adjust toNextXY for lack of innermost loop
        toNextXY += toNextLayer;
        for (destWaveEnd = destWave + tPoints; destWave < destWaveEnd; destWave++, srcWave += toNextXY){
            *destWave = *srcWave; // set destination to chosen layer in source
        }
    }else{ // a max or min projection
        switch (projMode) {
            case 0: // minimum projection
                for (destWaveEnd = destWave + tPoints; destWave < destWaveEnd; destWave++, srcWave += toNextXY){
                    *destWave = *srcWave; // set destination to first layer in source
                    // check remining layers in source to find min
                    for (lastLayer= srcWave + (layersToDo * toNextLayer), srcWave += toNextLayer ; srcWave < lastLayer; srcWave += toNextLayer){
                        if (*srcWave < *destWave)
                            *destWave = *srcWave;
                    }
                }
                break;
            case 1: // maximum projection
                for (destWaveEnd = destWave + tPoints; destWave < destWaveEnd; destWave++, srcWave += toNextXY){
                    *destWave = *srcWave; // set destination to first layer in source
                    // check remining layers in source to find max
                    for (lastLayer= srcWave + (layersToDo * toNextLayer), srcWave += toNextLayer ; srcWave < lastLayer; srcWave += toNextLayer){
                        if (*srcWave > *destWave)
                            *destWave = *srcWave;
                    }
                }
                break;
            case 2: // average projection
                double temp;
                for (destWaveEnd = destWave + tPoints; destWave < destWaveEnd; destWave++, srcWave += toNextXY){
                    for (lastLayer= srcWave + (layersToDo * toNextLayer), temp =0; srcWave < lastLayer; srcWave += toNextLayer){
                        temp += *srcWave;
                    }
                    *destWave = temp/layersToDo;
                }
                break;
            case 3: // median projection
                T *bufferStart = (T*) WMNewPtr (layersToDo * sizeof(T));
                T *bufferPos;
                for (destWaveEnd = destWave + tPoints; destWave < destWaveEnd; destWave++, srcWave += toNextXY){
                    for (lastLayer= srcWave + (layersToDo * toNextLayer), bufferPos = bufferStart; srcWave < lastLayer; srcWave += toNextLayer, bufferPos++){
                        *bufferPos = *srcWave;
                    }
                    *destWave = medianT (layersToDo, bufferStart);
                }
                WMDisposePtr ((Ptr)bufferStart);
                break;
        }
    }
}


/*  Structure to pass data to each ProjectThread
 Last Modified 2013/07/16 by Jamie Boyd */
typedef struct ProjectThreadParams{
    int inPutWaveType; // WM codes for wave types
    char* inPutDataStartPtr; // pointer to start of data in input wave
    char* outPutDataStartPtr; // pointer to start of data in output wave (may be same as outPut wave)
    UInt8 flatDim; // 0 for X projection, 1 for Y projection, 2 for Z projection
    UInt8 projMode; // 1 for a maximum intensity projection, 0 for a minimum intensity projection
    CountInt xSize; // size of X dimension (number of rows)
    CountInt ySize; // size of Y dimension (number of columns)
    CountInt zSize; // size of Z dimension (number of layers)
    CountInt startP; // starting point for projection
    CountInt endP; // end point for projection
    UInt8 ti; // number of this thread, starting from 0
    UInt8 tN; // total number of threads
} ProjectThreadParams, *ProjectThreadParamsPtr;


/* Each thread to do a maximum or minimum intensity projection through X,Y, or Z starts with this function
 Last Modified 2013/07/16 by Jamie Boyd */
void* ProjectThread (void* threadarg){
    struct ProjectThreadParams* p;
    p = (struct ProjectThreadParams*) threadarg;
    switch (p->flatDim){
        case 0:  //Project X
            switch(p->inPutWaveType){
                case NT_I8:
                    doProjectX((char*)p->inPutDataStartPtr, (char*)p->outPutDataStartPtr, p->projMode, p->xSize, p->ySize, p->zSize, p->startP, p->endP, p->ti, p->tN);
                    break;
                case (NT_I8 | NT_UNSIGNED):
                    doProjectX((unsigned char*)p->inPutDataStartPtr , (unsigned char*)p->outPutDataStartPtr,  p->projMode, p->xSize, p->ySize, p->zSize, p->startP, p->endP, p->ti, p->tN);
                    break;
                case NT_I16:
                    doProjectX((short*)p->inPutDataStartPtr, (short*)p->outPutDataStartPtr, p->projMode, p->xSize, p->ySize, p->zSize, p->startP, p->endP, p->ti, p->tN);
                    break;
                case (NT_I16 | NT_UNSIGNED):
                    doProjectX((unsigned short*)p->inPutDataStartPtr , (unsigned short*)p->outPutDataStartPtr,  p->projMode, p->xSize, p->ySize, p->zSize, p->startP, p->endP, p->ti, p->tN);
                    break;
                case NT_I32:
                    doProjectX((long*)p->inPutDataStartPtr, (long*)p->outPutDataStartPtr, p->projMode, p->xSize, p->ySize, p->zSize, p->startP, p->endP, p->ti, p->tN);
                    break;
                case (NT_I32| NT_UNSIGNED):
                    doProjectX((unsigned long*)p->inPutDataStartPtr, (unsigned long*)p->outPutDataStartPtr, p->projMode, p->xSize, p->ySize, p->zSize, p->startP, p->endP, p->ti, p->tN);
                    break;
                case NT_FP32:
                    doProjectX((float*)p->inPutDataStartPtr, (float*)p->outPutDataStartPtr, p->projMode, p->xSize, p->ySize, p->zSize, p->startP, p->endP, p->ti, p->tN);
                    break;
                case NT_FP64:
                    doProjectX((double*)p->inPutDataStartPtr, (double*)p->outPutDataStartPtr,  p->projMode, p->xSize, p->ySize, p->zSize, p->startP, p->endP, p->ti, p->tN);
                    break;
            }
            break;
        case 1: //Project Y
            switch(p->inPutWaveType){
                case NT_I8:
                    doProjectY((char*)p->inPutDataStartPtr, (char*)p->outPutDataStartPtr, p->projMode, p->xSize, p->ySize, p->zSize, p->startP, p->endP, p->ti, p->tN);
                    break;
                case (NT_I8 | NT_UNSIGNED):
                    doProjectY((unsigned char*)p->inPutDataStartPtr , (unsigned char*)p->outPutDataStartPtr,  p->projMode, p->xSize, p->ySize, p->zSize, p->startP, p->endP, p->ti, p->tN);
                    break;
                case NT_I16:
                    doProjectY((short*)p->inPutDataStartPtr, (short*)p->outPutDataStartPtr, p->projMode, p->xSize, p->ySize, p->zSize, p->startP, p->endP, p->ti, p->tN);
                    break;
                case (NT_I16 | NT_UNSIGNED):
                    doProjectY((unsigned short*)p->inPutDataStartPtr , (unsigned short*)p->outPutDataStartPtr,  p->projMode, p->xSize, p->ySize, p->zSize, p->startP, p->endP, p->ti, p->tN);
                    break;
                case NT_I32:
                    doProjectY((long*)p->inPutDataStartPtr, (long*)p->outPutDataStartPtr, p->projMode, p->xSize, p->ySize, p->zSize, p->startP, p->endP, p->ti, p->tN);
                    break;
                case (NT_I32| NT_UNSIGNED):
                    doProjectY((unsigned long*)p->inPutDataStartPtr, (unsigned long*)p->outPutDataStartPtr, p->projMode, p->xSize, p->ySize, p->zSize, p->startP, p->endP, p->ti, p->tN);
                    break;
                case NT_FP32:
                    doProjectY((float*)p->inPutDataStartPtr, (float*)p->outPutDataStartPtr, p->projMode, p->xSize, p->ySize, p->zSize, p->startP, p->endP, p->ti, p->tN);
                    break;
                case NT_FP64:
                    doProjectY((double*)p->inPutDataStartPtr, (double*)p->outPutDataStartPtr,  p->projMode, p->xSize, p->ySize, p->zSize, p->startP, p->endP, p->ti, p->tN);
                    break;
            }
            break;
        case 2: // Project Z
            switch(p->inPutWaveType){
                case NT_I8:
                    doProjectZ((char*)p->inPutDataStartPtr, (char*)p->outPutDataStartPtr, p->projMode, p->xSize, p->ySize, p->startP, p->endP, p->ti, p->tN);
                    break;
                case (NT_I8 | NT_UNSIGNED):
                    doProjectZ((unsigned char*)p->inPutDataStartPtr , (unsigned char*)p->outPutDataStartPtr,  p->projMode, p->xSize, p->ySize, p->startP, p->endP, p->ti, p->tN);
                    break;
                case NT_I16:
                    doProjectZ((short*)p->inPutDataStartPtr, (short*)p->outPutDataStartPtr, p->projMode, p->xSize, p->ySize, p->startP, p->endP, p->ti, p->tN);
                    break;
                case (NT_I16 | NT_UNSIGNED):
                    doProjectZ((unsigned short*)p->inPutDataStartPtr , (unsigned short*)p->outPutDataStartPtr,  p->projMode, p->xSize, p->ySize, p->startP, p->endP, p->ti, p->tN);
                    break;
                case NT_I32:
                    doProjectZ((long*)p->inPutDataStartPtr, (long*)p->outPutDataStartPtr, p->projMode, p->xSize, p->ySize, p->startP, p->endP, p->ti, p->tN);
                    break;
                case (NT_I32| NT_UNSIGNED):
                    doProjectZ((unsigned long*)p->inPutDataStartPtr, (unsigned long*)p->outPutDataStartPtr, p->projMode, p->xSize, p->ySize, p->startP, p->endP, p->ti, p->tN);
                    break;
                case NT_FP32:
                    doProjectZ((float*)p->inPutDataStartPtr, (float*)p->outPutDataStartPtr, p->projMode, p->xSize, p->ySize, p->startP, p->endP, p->ti, p->tN);
                    break;
                case NT_FP64:
                    doProjectZ((double*)p->inPutDataStartPtr, (double*)p->outPutDataStartPtr,  p->projMode, p->xSize, p->ySize, p->startP, p->endP, p->ti, p->tN);
                    break;
            }
            break;
    }
    return 0;
}


/* ProjectSpecFrames XOP entry function
 Makes a projection image of a specified range of columns or rows or layers in the input wave and puts
 the result in a specified layer of the output wave.
 ProjectSpecFramesParams
 double projMode            0 = mimimum intensity projection, 1 =maximum intensity projection, 2 = avg, 3 = median
 double flatDimension       the dimension that we are projecting along
 double outPutLayer         the layer in the output wave that receives the projection
 waveHndl outPutWaveH       A handle to output wave
 double inPutEndLayer       end of range of layers to project
 double inPutStartLayer     start of range of layers to project
 waveHndl inPutWaveH        handle to the input wave
 Last Modified 2025/06/25 by Jamie Boyd */
extern "C" int ProjectSpecFrames (ProjectSpecFramesParamsPtr p){
    int result =0;                                  // The error returned from various Wavemetrics functions
    waveHndl inPutWaveH = NIL, outPutWaveH = NIL;   // Handles to the input and output waves
    int inPutWaveType, outPutWaveType;              // Wavetypes numeric codes for things like 32 bit floating point, 16 bit int, etc
    int inPutDimensions,outPutDimensions;           // number of numDimensions in input and output waves
    CountInt inPutDimensionSizes[MAX_DIMENSIONS+1]; // an array used to hold the width, height, layers, and chunk sizes of input wave
    CountInt outPutDimensionSizes[MAX_DIMENSIONS+1]; //an array used to hold the width, height, layers, and chunk sizes of output wave
    CountInt inPutWaveOffset, outPutWaveOffset;      //offset in bytes from begnning of handle to a wave to the actual data
    char* srcWaveStart, *destWaveStart;              // Pointers to start of data in the inut and output waves.
    CountInt startP = p->inPutStartLayer;            // starting Z (or X or Y) positon for projection
    CountInt endP = p->inPutEndLayer;                // ending Z (or X or Y) positon for projection
    CountInt outPutP = p->outPutLayer;               // Z (or X or Y) layer in output wave that get result of projection
    UInt8 flatDimension = p->flatDimension;          // dimension along which flattening occurs. 0 = X, 1 =Y, 2 = Z
    UInt8 iThread;                                   // number of each thread, starting from 0
    UInt8 nThreads;                                  // total number of threads
    ProjectThreadParamsPtr paramArrayPtr = nullptr;  // array of params for threading
    pthread_t* threadsPtr;                           // array of pthreads
    try {
        // check projection mode
        if ((p->projMode < 0) || (p->projMode > 3)) throw result = BADDSTYPE;
        // Get handle to input and output waves make sure they exist.
        inPutWaveH = p->inPutWaveH;
        if (inPutWaveH == NIL) throw result= NON_EXISTENT_WAVE;
        outPutWaveH = p->outPutWaveH;
        if (outPutWaveH == NIL) throw result = NON_EXISTENT_WAVE;
        // Get waves data type
        inPutWaveType = WaveType(inPutWaveH);
        outPutWaveType =  WaveType(outPutWaveH);
        if (inPutWaveType != outPutWaveType) throw result = NOTSAMEWAVETYPE;
        // Check that we don't have a text wave
        if ((inPutWaveType==TEXT_WAVE_TYPE) ||  (outPutWaveType==TEXT_WAVE_TYPE)) throw result = NOTEXTWAVES;
        // Get number of used dimensions in waves.
        if (MDGetWaveDimensions(inPutWaveH, &inPutDimensions, inPutDimensionSizes)) throw result = WAVEERROR_NOS;
        if (MDGetWaveDimensions(outPutWaveH, &outPutDimensions, outPutDimensionSizes)) throw result = WAVEERROR_NOS;
        // Check that input wave is 3D and output wave is 2D or 3D
        if (inPutDimensions != 3) throw result = INPUTNEEDS_3D_WAVE;
        if ((outPutDimensions != 2) && (outPutDimensions != 3)) throw result = OUTPUTNEEDS_2D3D_WAVE;
        // Check that the dimensionality of the input and output waves match up for the dimension we are trying to collapse
        switch (flatDimension){
            case 0:    // X projection  output is dim [0] = y-size, dim [1] = z-Size
                if ((outPutDimensionSizes[0] != inPutDimensionSizes [1]) || (outPutDimensionSizes[1] != inPutDimensionSizes [2])) throw result = NOTSAMEDIMSIZE;
                break;
            case 1:    // Y projection  output is dim [0] = x-size, dim [1] = z-Size
                if ((outPutDimensionSizes[0] != inPutDimensionSizes [0]) || (outPutDimensionSizes[1] != inPutDimensionSizes [2]))
                    throw result = NOTSAMEDIMSIZE;
                break;
            case 2:    //Z projection  output is dim [0] = x-size, dim [1] = y-Size
                if ((outPutDimensionSizes[0] != inPutDimensionSizes [0]) || (outPutDimensionSizes[1] != inPutDimensionSizes [1]))
                    throw result = NOTSAMEDIMSIZE;
                break;
            default:
                throw result = BADDIMENSION;
                break;
        }
        // check that range of input dimension is valid
        if ((startP < 0) || (startP >= inPutDimensionSizes [flatDimension])) throw result = INPUT_RANGE;
        if ((endP < 0) || (endP >= inPutDimensionSizes [flatDimension])) throw result = INPUT_RANGE;
        if (startP > endP) throw result = INPUT_RANGE;
        //check that output layer is within the size of the output wave
        if (outPutDimensions == 2){  //output is 2D wave so frame has to be 0
            if (outPutP != 0) throw result = INVALIDOUTPUTFRAME;
        }else{
            if ((outPutP < 0) || (outPutP >= outPutDimensionSizes [2])) throw result = INVALIDOUTPUTFRAME;
        }
        // Get the offsets to the data in the input
        if (MDAccessNumericWaveData(inPutWaveH, kMDWaveAccessMode0, &inPutWaveOffset)) throw result = WAVEERROR_NOS;
        srcWaveStart = (char*)(*inPutWaveH) + inPutWaveOffset;
        if (MDAccessNumericWaveData(outPutWaveH, kMDWaveAccessMode0, &outPutWaveOffset)) throw result = WAVEERROR_NOS;
        destWaveStart = (char*)(*outPutWaveH) + outPutWaveOffset;
        // make an array of parameter structures
        nThreads = gNumProcessors;
        paramArrayPtr = (ProjectThreadParamsPtr)WMNewPtr (nThreads * sizeof(ProjectThreadParams));
        if (paramArrayPtr == nullptr) throw result = MEMFAIL;
        // make an array of pthread_t
        threadsPtr =(pthread_t*)WMNewPtr(nThreads * sizeof(pthread_t));
        if (threadsPtr == nullptr) throw result = MEMFAIL;
    }catch (int result){
        if (threadsPtr != nullptr) WMDisposePtr ((Ptr)threadsPtr);
        if (paramArrayPtr != nullptr) WMDisposePtr ((Ptr)paramArrayPtr);
        p -> result = (double)(result - FIRST_XOP_ERR);
#ifdef NO_IGOR_ERR
        return (0);
#else
        return (result);
#endif
    }
    // fill paramater array
    switch(outPutWaveType){
        case NT_I8:
        case (NT_I8 | NT_UNSIGNED):
            destWaveStart += sizeof(char) * (outPutP * outPutDimensionSizes[0] * outPutDimensionSizes [1]);
            break;
        case NT_I16:
        case (NT_I16 | NT_UNSIGNED):
            destWaveStart += sizeof(short) * (outPutP * outPutDimensionSizes[0] * outPutDimensionSizes [1]);
            break;
        case NT_I32:
        case (NT_I32 | NT_UNSIGNED):
            destWaveStart += sizeof(long) * (outPutP * outPutDimensionSizes[0] * outPutDimensionSizes [1]);
            break;
        case NT_FP32:
            destWaveStart += sizeof(float) * (outPutP * outPutDimensionSizes[0] * outPutDimensionSizes [1]);
            break;
        case NT_FP64:
            destWaveStart += sizeof(double) * (outPutP * outPutDimensionSizes[0] * outPutDimensionSizes [1]);
            break;
    }
    for (iThread = 0; iThread < nThreads; iThread++){
        paramArrayPtr[iThread].inPutWaveType = inPutWaveType;
        paramArrayPtr[iThread].inPutDataStartPtr = srcWaveStart;
        paramArrayPtr[iThread].outPutDataStartPtr = destWaveStart;
        paramArrayPtr[iThread].flatDim = flatDimension;
        paramArrayPtr[iThread].projMode = p->projMode;
        paramArrayPtr[iThread].xSize = inPutDimensionSizes [0];
        paramArrayPtr[iThread].ySize = inPutDimensionSizes [1];
        paramArrayPtr[iThread].zSize =inPutDimensionSizes [2];
        paramArrayPtr[iThread].startP = startP;
        paramArrayPtr[iThread].endP= endP;
        paramArrayPtr[iThread].ti=iThread; // number of this thread, starting from 0
        paramArrayPtr[iThread].tN =nThreads; // total number of threads
    }
    // create the threads
    for (iThread = 0; iThread < nThreads; iThread++){
        pthread_create (&threadsPtr[iThread], NULL, ProjectThread, (void *) &paramArrayPtr[iThread]);
    }
    //join threads
    for (iThread = 0; iThread < nThreads; iThread++){
        pthread_join (threadsPtr[iThread], NULL);
    }
    WMDisposePtr ((Ptr)threadsPtr);     // free memory for pThreads Array
    WMDisposePtr ((Ptr)paramArrayPtr);   // Free paramaterArray memory
    WaveHandleModified(outPutWaveH);     // Inform Igor that we have changed the output wave.
    p->result= (0);
    return (0);
}

/* ProjectAllFrames XOP entry function
 Makes a projection Image of all columns or rows or layers in the input wave and makes a 2D wave to receive result
 ProjectAllFramesParams
 waveHndl inPutWaveH         handle to the input wave
 Handle outPutPath           A handle to a string containing path to output wave we want to make
 double flatDimension        Which dimension we want to collapse on, 0 for x, 1 for y, 2 for z
 double overwrite            0 to give errors when wave already exists. non-zero to cheerfully overwrite existing wave.
 double projMode             0 = mimimum intensity projection, 1 =maximum intensity projection, 2 = avg, 3 = median
 Last Modified 2025/06/25 by Jamie Boyd */
extern "C" int ProjectAllFrames (ProjectAllFramesParamsPtr p){
    int result = 0;                                     // for error codes
    waveHndl inPutWaveH, outPutWaveH;                   // Handles to the input and output waves
    int waveType;                                       //  Wavetypes numeric codes for things like 32 bit floating point, 16 bit int, etc
    int numDimensions;                                  //number of Dimensions in input and output waves
    CountInt inPutDimensionSizes[MAX_DIMENSIONS+1];     // an array used to hold wave width, height, layers, and chunk sizes of input wave
    CountInt outPutDimensionSizes[MAX_DIMENSIONS+1];    // an array used to hold wave width, height, layers, and chunk sizes of output wave
    UInt16 outPutPathLen;                               // Length of the path to the target folder (output path - wave name)
    DataFolderHandle outPutDFHandle=nullptr;            // Handle to the datafolder where we will put the output wave
    DataFolderHandle inPutDFHandle=nullptr;             // Handle to datafolder for input wave, used to test if overwriting a wave
    DFPATH inPutPath, outPutPath;                       // C string to hold data folder path of output wave
    WVNAME inPutWaveName, outPutWaveName;               // C string to hold name of output wave
    CountInt inPutWaveOffset, outPutWaveOffset;         //offset in bytes from begnning of handle to a wave to the actual data
    char* srcWaveStart, *destWaveStart;                 // Pointers to start of data in the inut and output waves.
    UInt8 flatDimension = p->flatDimension;             // dimension along which flattening occurs. 0 = X, 1 =Y, 2 = Z
    UInt8 flatten = 0;                                  // set if overwriting input wave, i.e., we flatten it
    CountInt endP;                                      // start point is always = 0, end point is dimension size - 1
    UInt8 iThread;                                      // number of each thread, starting from 0
    UInt8 nThreads;                                     // total number of threads
    ProjectThreadParamsPtr paramArrayPtr = nullptr;     // array of params for threading
    pthread_t* threadsPtr = nullptr;                    // array of pthreads
    int overWrite= p->overwrite;
    try{
        // check projection mode
        if ((p->projMode < 0) || (p->projMode > 3)) throw result = BADDSTYPE;
        // Get handle to input wave make sure it exists.
        inPutWaveH = p->inPutWaveH;
        if (inPutWaveH == NIL) throw result = NON_EXISTENT_WAVE;
        // Get wave data type
        waveType = WaveType(inPutWaveH);
        // Check that we don't have a text wave
        if (waveType==TEXT_WAVE_TYPE) throw result = NOTEXTWAVES;
        // Get number of used numDimensions in input wave.
        if (MDGetWaveDimensions(inPutWaveH, &numDimensions, inPutDimensionSizes)) throw result = WAVEERROR_NOS;
        // Check that input wave is 3D
        if (numDimensions != 3) throw result = INPUTNEEDS_3D_WAVE;
        // Set the sizes of the output wave as appropriate for the dimension we want to squish
        switch (flatDimension){
            case 0:    // X projection  output is dim [0] = y-size, dim [1] = z-Size
                outPutDimensionSizes [0] = inPutDimensionSizes [1];
                outPutDimensionSizes [1] = inPutDimensionSizes [2];
                endP = inPutDimensionSizes [0] - 1;
                break;
            case 1:// Y projection  output is dim [0] = x-size, dim [1] = z-Size
                outPutDimensionSizes [0] = inPutDimensionSizes [0];
                outPutDimensionSizes [1] = inPutDimensionSizes [2];
                endP = inPutDimensionSizes [1] - 1;
                break;
            case 2:    //z dimension
                outPutDimensionSizes [0] = inPutDimensionSizes [0];
                outPutDimensionSizes [1] = inPutDimensionSizes [1];
                endP = inPutDimensionSizes [2] - 1;
                break;
            default:
                throw result = BADDIMENSION;
                break;
        }
        outPutDimensionSizes [2] = 0;
        outPutDimensionSizes [3] = 0;
        // If outPutPath is empty string, we are overwriting existing wave
        outPutPathLen = WMGetHandleSize (p->outPutPath);
        if (outPutPathLen == 0){
            if (overWrite == NO_OVERWITE) throw result = OVERWRITEALERT;
            outPutWaveH = inPutWaveH;
            flatten = 1;
        }else{ // Parse outPut path
            ParseWavePath (p->outPutPath, outPutPath, outPutWaveName);
            // Clean up wave name: no liberal names
            CleanupName (0, outPutWaveName, MAX_OBJ_NAME);
            //check that data folder is valid and get a handle to the datafolder
            if (GetNamedDataFolder (NULL, outPutPath, &outPutDFHandle))throw result = WAVEERROR_NOS;
            //Test for overwriting
            WaveName (inPutWaveH, inPutWaveName);
            GetWavesDataFolder (inPutWaveH, &inPutDFHandle);
            GetDataFolderNameOrPath (inPutDFHandle, 1, inPutPath);
            if ((CmpStr (inPutPath,outPutPath) ==0) && (CmpStr (inPutWaveName,outPutWaveName) ==0)){    // Then we would overrite input wave
                if (overWrite == NO_OVERWITE) throw result = OVERWRITEALERT;
                flatten = 1;
                outPutWaveH = inPutWaveH;
            }else{
                // make the output wave
                if (MDMakeWave (&outPutWaveH, outPutWaveName, outPutDFHandle, outPutDimensionSizes,waveType, overWrite)) throw result = WAVEERROR_NOS;
            }
        }
        // Get the offsets to the data in the input and output waves
        if (MDAccessNumericWaveData(inPutWaveH, kMDWaveAccessMode0, &inPutWaveOffset)) throw result = WAVEERROR_NOS;
        srcWaveStart = (char*)(*inPutWaveH) + inPutWaveOffset;
        if (MDAccessNumericWaveData(outPutWaveH, kMDWaveAccessMode0, &outPutWaveOffset)) throw result = WAVEERROR_NOS;
        destWaveStart = (char*)(*outPutWaveH) + outPutWaveOffset;
        // make an array of parameter structures
        nThreads = gNumProcessors;
        paramArrayPtr = (ProjectThreadParamsPtr)WMNewPtr (nThreads * sizeof(ProjectThreadParams));
        if (paramArrayPtr == nullptr) throw result = MEMFAIL;
        // make an array of pthread_t
        threadsPtr =(pthread_t*)WMNewPtr(nThreads * sizeof(pthread_t));
        if (threadsPtr == nullptr) throw result = MEMFAIL;
    }catch (int (result)){
        if (threadsPtr != nullptr) WMDisposePtr ((Ptr)threadsPtr);
        if (paramArrayPtr != nullptr) WMDisposePtr ((Ptr)paramArrayPtr);
        p -> result = (double)(result - FIRST_XOP_ERR);
#ifdef NO_IGOR_ERR
        return (0);
#else
        return (result);
#endif
    }
    // fill param array struct
    for (iThread = 0; iThread < nThreads; iThread++){
        paramArrayPtr[iThread].inPutWaveType = waveType;
        paramArrayPtr[iThread].inPutDataStartPtr = srcWaveStart;
        paramArrayPtr[iThread].outPutDataStartPtr = destWaveStart;
        paramArrayPtr[iThread].flatDim = flatDimension;
        paramArrayPtr[iThread].projMode = p->projMode;
        paramArrayPtr[iThread].xSize = inPutDimensionSizes [0];
        paramArrayPtr[iThread].ySize = inPutDimensionSizes [1];
        paramArrayPtr[iThread].zSize =inPutDimensionSizes [2];
        paramArrayPtr[iThread].startP = 0;
        paramArrayPtr[iThread].endP= endP;
        paramArrayPtr[iThread].ti=iThread; // number of this thread, starting from 0
        paramArrayPtr[iThread].tN =nThreads; // total number of threads
    }
    // create the threads
    for (iThread = 0; iThread < nThreads; iThread++){
        pthread_create (&threadsPtr[iThread], NULL, ProjectThread, (void *) &paramArrayPtr[iThread]);
    }
    //join the threads
    for (iThread = 0; iThread < nThreads; iThread++){
        pthread_join (threadsPtr[iThread], NULL);
    }
    WMDisposePtr ((Ptr)threadsPtr);     // free memory for pThreads Array
    WMDisposePtr ((Ptr)paramArrayPtr);   // Free paramaterArray memory
    if (flatten){
        MDChangeWave (outPutWaveH, -1, outPutDimensionSizes);
    }
    WaveHandleModified(outPutWaveH);     // Inform Igor that we have changed the output wave.
    if (p->outPutPath)
        WMDisposeHandle(p->outPutPath);
    p->result = (0);
    return (0);
}

/* ProjectXSlice XOP entry function
 Gets an X slice from a 3D wave and puts it a pre-existing 2D wave of the right dimensions
 typedef struct ProjectSliceParams
 waveHndl inPutWaveH    handle to the input wave
 waveHndl outPutWaveH   handle to the output wave
 double slice           X, Y or Z slice to get
 Last Modified: 2025/06/25 by Jamie Boyd */
extern "C" int ProjectXSlice (ProjectSliceParamsPtr p) {
    int result =0;                                  // The error returned from various Wavemetrics functions
    waveHndl inPutWaveH = NIL, outPutWaveH = NIL;   // Handles to the input and output waves
    int inPutWaveType, outPutWaveType;              // Wavetypes numeric codes for things like 32 bit floating point, 16 bit int, etc
    int inPutDimensions,outPutDimensions;           // number of numDimensions in input and output waves
    CountInt inPutDimensionSizes[MAX_DIMENSIONS+1]; // an array used to hold the width, height, layers, and chunk sizes of input wave
    CountInt outPutDimensionSizes[MAX_DIMENSIONS+1]; //an array used to hold the width, height, layers, and chunk sizes of output wave
    CountInt inPutWaveOffset, outPutWaveOffset;      //offset in bytes from begnning of handle to a wave to the actual data
    char* srcWaveStart, *destWaveStart;              // Pointers to start of data in the inut and output waves.
    CountInt slice =  p->slice;                      // input layer
    UInt8 iThread;                                   // number of each thread, starting from 0
    UInt8 nThreads;                                  // total number of threads
    ProjectThreadParamsPtr paramArrayPtr = nullptr;  // array of params for threading
    pthread_t* threadsPtr = nullptr;                 // array of pthreads
    try {
        // Get handle to input and output waves make sure they exist.
        inPutWaveH = p->inPutWaveH;
        if (inPutWaveH == NIL) throw result= NON_EXISTENT_WAVE;
        outPutWaveH = p->outPutWaveH;
        if (outPutWaveH == NIL) throw result = NON_EXISTENT_WAVE;
        // Get waves data type
        inPutWaveType = WaveType(inPutWaveH);
        outPutWaveType =  WaveType(outPutWaveH);
        if (inPutWaveType != outPutWaveType) throw result = NOTSAMEWAVETYPE;
        // Check that we don't have a text wave
        if ((inPutWaveType==TEXT_WAVE_TYPE) || (outPutWaveType==TEXT_WAVE_TYPE)) throw result = NOTEXTWAVES;
        // Get number of used dimensions in waves.
        if (MDGetWaveDimensions(inPutWaveH, &inPutDimensions, inPutDimensionSizes)) throw result = WAVEERROR_NOS;
        if (MDGetWaveDimensions(outPutWaveH, &outPutDimensions, outPutDimensionSizes)) throw result = WAVEERROR_NOS;
        // Check that input wave is 3D and output wave is 2D or 3D
        if (inPutDimensions != 3) throw result = INPUTNEEDS_3D_WAVE;
        if ((outPutDimensions != 2) && (outPutDimensions != 3)) throw result = OUTPUTNEEDS_2D3D_WAVE;
        // Check that the dimensionality of the input and output waves match up for the dimension we are trying to collapse
        if ((outPutDimensionSizes[0] != inPutDimensionSizes [1]) || (outPutDimensionSizes[1] != inPutDimensionSizes [2])) throw result = NOTSAMEDIMSIZE;
        // check that slice is valid
        if ((slice < 0) || (slice >= inPutDimensionSizes [0])) throw result = INPUT_RANGE;
        // Get the offsets to the data in the input
        if (MDAccessNumericWaveData(inPutWaveH, kMDWaveAccessMode0, &inPutWaveOffset)) throw result = WAVEERROR_NOS;
        srcWaveStart = (char*)(*inPutWaveH) + inPutWaveOffset;
        if (MDAccessNumericWaveData(outPutWaveH, kMDWaveAccessMode0, &outPutWaveOffset)) throw result = WAVEERROR_NOS;
        destWaveStart = (char*)(*outPutWaveH) + outPutWaveOffset;
        // make an array of parameter structures
        nThreads = gNumProcessors;
        paramArrayPtr = (ProjectThreadParamsPtr)WMNewPtr (nThreads * sizeof(ProjectThreadParams));
        if (paramArrayPtr == nullptr) throw result = MEMFAIL;
        // make an array of pthread_t
        threadsPtr =(pthread_t*)WMNewPtr(nThreads * sizeof(pthread_t));
        if (threadsPtr == nullptr) throw result = MEMFAIL;
    }catch (int result){
        if (threadsPtr != nullptr) WMDisposePtr ((Ptr)threadsPtr);
        if (paramArrayPtr != nullptr) WMDisposePtr ((Ptr)paramArrayPtr);
        p -> result = (double)(result - FIRST_XOP_ERR);
        return result;
    }
    // fill param array struct
    for (iThread = 0; iThread < nThreads; iThread++){
        paramArrayPtr[iThread].inPutWaveType = inPutWaveType;
        paramArrayPtr[iThread].inPutDataStartPtr = srcWaveStart;
        paramArrayPtr[iThread].outPutDataStartPtr = destWaveStart;
        paramArrayPtr[iThread].flatDim = 0;
        paramArrayPtr[iThread].projMode = 0;
        paramArrayPtr[iThread].xSize = inPutDimensionSizes [0];
        paramArrayPtr[iThread].ySize = inPutDimensionSizes [1];
        paramArrayPtr[iThread].zSize =inPutDimensionSizes [2];
        paramArrayPtr[iThread].startP = slice;
        paramArrayPtr[iThread].endP= slice;
        paramArrayPtr[iThread].ti=iThread; // number of this thread, starting from 0
        paramArrayPtr[iThread].tN =nThreads; // total number of threads
    }
    // create the threads
    for (iThread = 0; iThread < nThreads; iThread++){
        pthread_create (&threadsPtr[iThread], NULL, ProjectThread, (void *) &paramArrayPtr[iThread]);
    }
    //join the threads
    for (iThread = 0; iThread < nThreads; iThread++){
        pthread_join (threadsPtr[iThread], NULL);
    }
    WMDisposePtr ((Ptr)threadsPtr);     // free memory for pThreads Array
    WMDisposePtr ((Ptr)paramArrayPtr);   // Free paramaterArray memory
    WaveHandleModified(outPutWaveH);     // Inform Igor that we have changed the output wave.
    p->result= (0);
    return (0);
}


/* ProjectYSlice XOP entry function
 Gets a Y slice from a 3D wave and puts it a pre-existing 2D wave of the right dimensions
 Last Modified: 2025/06/25 by Jamie Boyd */
extern "C" int ProjectYSlice (ProjectSliceParamsPtr p){
    int result =0;                                  // The error returned from various Wavemetrics functions
    waveHndl inPutWaveH = NIL, outPutWaveH = NIL;   // Handles to the input and output waves
    int inPutWaveType, outPutWaveType;              // Wavetypes numeric codes for things like 32 bit floating point, 16 bit int, etc
    int inPutDimensions,outPutDimensions;           // number of numDimensions in input and output waves
    CountInt inPutDimensionSizes[MAX_DIMENSIONS+1]; // an array used to hold the width, height, layers, and chunk sizes of input wave
    CountInt outPutDimensionSizes[MAX_DIMENSIONS+1]; //an array used to hold the width, height, layers, and chunk sizes of output wave
    CountInt inPutWaveOffset, outPutWaveOffset;      //offset in bytes from begnning of handle to a wave to the actual data
    char* srcWaveStart, *destWaveStart;              // Pointers to start of data in the inut and output waves.
    CountInt slice =  p->slice;                      // input layer
    UInt8 iThread;                                   // number of each thread, starting from 0
    UInt8 nThreads;                                  // total number of threads
    ProjectThreadParamsPtr paramArrayPtr = nullptr;  // array of params for threading
    pthread_t* threadsPtr = nullptr;                 // array of pthreads
    try {
        // Get handle to input and output waves make sure they exist.
        inPutWaveH = p->inPutWaveH;
        if (inPutWaveH == NIL) throw result= NON_EXISTENT_WAVE;
        outPutWaveH = p->outPutWaveH;
        if (outPutWaveH == NIL) throw result = NON_EXISTENT_WAVE;
        // Get waves data type
        inPutWaveType = WaveType(inPutWaveH);
        outPutWaveType =  WaveType(outPutWaveH);
        if (inPutWaveType != outPutWaveType) throw result = NOTSAMEWAVETYPE;
        // Check that we don't have a text wave
        if ((inPutWaveType==TEXT_WAVE_TYPE) || (outPutWaveType==TEXT_WAVE_TYPE)) throw result = NOTEXTWAVES;
        // Get number of used dimensions in waves.
        if (MDGetWaveDimensions(inPutWaveH, &inPutDimensions, inPutDimensionSizes)) throw result = WAVEERROR_NOS;
        if (MDGetWaveDimensions(outPutWaveH, &outPutDimensions, outPutDimensionSizes)) throw result = WAVEERROR_NOS;
        // Check that input wave is 3D and output wave is 2D or 3D
        if (inPutDimensions != 3) throw result = INPUTNEEDS_3D_WAVE;
        if ((outPutDimensions != 2) && (outPutDimensions != 3)) throw result = OUTPUTNEEDS_2D3D_WAVE;
        // Check that the dimensionality of the input and output waves match up for the dimension we are trying to collapse
        if ((outPutDimensionSizes[0] != inPutDimensionSizes [0]) || (outPutDimensionSizes[1] != inPutDimensionSizes [2])) throw result = NOTSAMEDIMSIZE;
        // check that slice is valid in range of Y
        if ((slice < 0) || (slice >= inPutDimensionSizes [1])) throw result = INPUT_RANGE;
        // Get the offsets to the data in the input
        if (MDAccessNumericWaveData(inPutWaveH, kMDWaveAccessMode0, &inPutWaveOffset)) throw result = WAVEERROR_NOS;
        srcWaveStart = (char*)(*inPutWaveH) + inPutWaveOffset;
        if (MDAccessNumericWaveData(outPutWaveH, kMDWaveAccessMode0, &outPutWaveOffset)) throw result = WAVEERROR_NOS;
        destWaveStart = (char*)(*outPutWaveH) + outPutWaveOffset;
        // make an array of parameter structures
        nThreads = gNumProcessors;
        paramArrayPtr = (ProjectThreadParamsPtr)WMNewPtr (nThreads * sizeof(ProjectThreadParams));
        if (paramArrayPtr == nullptr) throw result = MEMFAIL;
        // make an array of pthread_t
        threadsPtr =(pthread_t*)WMNewPtr(nThreads * sizeof(pthread_t));
        if (threadsPtr == nullptr) throw result = MEMFAIL;
    }catch (int result){
        if (threadsPtr != nullptr) WMDisposePtr ((Ptr)threadsPtr);
        if (paramArrayPtr != nullptr) WMDisposePtr ((Ptr)paramArrayPtr);
        p -> result = (double)(result - FIRST_XOP_ERR);
        return result;
    }
    // fill param array struct
    for (iThread = 0; iThread < nThreads; iThread++){
        paramArrayPtr[iThread].inPutWaveType = inPutWaveType;
        paramArrayPtr[iThread].inPutDataStartPtr = srcWaveStart;
        paramArrayPtr[iThread].outPutDataStartPtr = destWaveStart;
        paramArrayPtr[iThread].flatDim = 1;
        paramArrayPtr[iThread].projMode = 0;
        paramArrayPtr[iThread].xSize = inPutDimensionSizes [0];
        paramArrayPtr[iThread].ySize = inPutDimensionSizes [1];
        paramArrayPtr[iThread].zSize =inPutDimensionSizes [2];
        paramArrayPtr[iThread].startP = slice;
        paramArrayPtr[iThread].endP= slice;
        paramArrayPtr[iThread].ti=iThread; // number of this thread, starting from 0
        paramArrayPtr[iThread].tN =nThreads; // total number of threads
    }
    // create the threads
    for (iThread = 0; iThread < nThreads; iThread++){
        pthread_create (&threadsPtr[iThread], NULL, ProjectThread, (void *) &paramArrayPtr[iThread]);
    }
    //join the threads
    for (iThread = 0; iThread < nThreads; iThread++){
        pthread_join (threadsPtr[iThread], NULL);
    }
    WMDisposePtr ((Ptr)threadsPtr);     // free memory for pThreads Array
    WMDisposePtr ((Ptr)paramArrayPtr);   // Free paramaterArray memory
    WaveHandleModified(outPutWaveH);     // Inform Igor that we have changed the output wave.
    p->result= (0);
    return (0);
}

/* ProjectZSlice XOP entry function
 Gets a Z slice from a 3D wave and puts it a pre-existing 2D wave of the right dimensions
 Last Modified: 2025/06/25 by Jamie Boyd */
extern "C" int ProjectZSlice (ProjectSliceParamsPtr p){
    int result =0;                                  // The error returned from various Wavemetrics functions
    waveHndl inPutWaveH = NIL, outPutWaveH = NIL;   // Handles to the input and output waves
    int inPutWaveType, outPutWaveType;              // Wavetypes numeric codes for things like 32 bit floating point, 16 bit int, etc
    int inPutDimensions,outPutDimensions;           // number of numDimensions in input and output waves
    CountInt inPutDimensionSizes[MAX_DIMENSIONS+1]; // an array used to hold the width, height, layers, and chunk sizes of input wave
    CountInt outPutDimensionSizes[MAX_DIMENSIONS+1]; //an array used to hold the width, height, layers, and chunk sizes of output wave
    CountInt inPutWaveOffset, outPutWaveOffset;      //offset in bytes from begnning of handle to a wave to the actual data
    char* srcWaveStart, *destWaveStart;              // Pointers to start of data in the inut and output waves.
    CountInt slice =  p->slice;                      // input layer
    UInt8 iThread;                                   // number of each thread, starting from 0
    UInt8 nThreads;                                  // total number of threads
    ProjectThreadParamsPtr paramArrayPtr = nullptr;  // array of params for threading
    pthread_t* threadsPtr = nullptr;                 // array of pthreads
    try {
        // Get handle to input and output waves make sure they exist.
        inPutWaveH = p->inPutWaveH;
        if (inPutWaveH == NIL) throw result= NON_EXISTENT_WAVE;
        outPutWaveH = p->outPutWaveH;
        if (outPutWaveH == NIL) throw result = NON_EXISTENT_WAVE;
        // Get waves data type
        inPutWaveType = WaveType(inPutWaveH);
        outPutWaveType =  WaveType(outPutWaveH);
        if (inPutWaveType != outPutWaveType) throw result = NOTSAMEWAVETYPE;
        // Check that we don't have a text wave
        if ((inPutWaveType==TEXT_WAVE_TYPE) || (outPutWaveType==TEXT_WAVE_TYPE)) throw result = NOTEXTWAVES;
        // Get number of used dimensions in waves.
        if (MDGetWaveDimensions(inPutWaveH, &inPutDimensions, inPutDimensionSizes)) throw result = WAVEERROR_NOS;
        if (MDGetWaveDimensions(outPutWaveH, &outPutDimensions, outPutDimensionSizes)) throw result = WAVEERROR_NOS;
        // Check that input wave is 3D and output wave is 2D or 3D
        if (inPutDimensions != 3) throw result = INPUTNEEDS_3D_WAVE;
        if ((outPutDimensions != 2) && (outPutDimensions != 3)) throw result = OUTPUTNEEDS_2D3D_WAVE;
        // Check that the dimensionality of the input and output waves match up for the dimension we are trying to collapse
        if ((outPutDimensionSizes[0] != inPutDimensionSizes [0]) || (outPutDimensionSizes[1] != inPutDimensionSizes [1])) throw result = NOTSAMEDIMSIZE;
        // check that slice is valid in range of Z
        if ((slice < 0) || (slice >= inPutDimensionSizes [2])) throw result = INPUT_RANGE;
        // Get the offsets to the data in the input
        if (MDAccessNumericWaveData(inPutWaveH, kMDWaveAccessMode0, &inPutWaveOffset)) throw result = WAVEERROR_NOS;
        srcWaveStart = (char*)(*inPutWaveH) + inPutWaveOffset;
        if (MDAccessNumericWaveData(outPutWaveH, kMDWaveAccessMode0, &outPutWaveOffset)) throw result = WAVEERROR_NOS;
        destWaveStart = (char*)(*outPutWaveH) + outPutWaveOffset;
        // make an array of parameter structures
        nThreads = gNumProcessors;
        paramArrayPtr = (ProjectThreadParamsPtr)WMNewPtr (nThreads * sizeof(ProjectThreadParams));
        if (paramArrayPtr == nullptr) throw result = MEMFAIL;
        // make an array of pthread_t
        threadsPtr =(pthread_t*)WMNewPtr(nThreads * sizeof(pthread_t));
        if (threadsPtr == nullptr) throw result = MEMFAIL;
    }catch (int result){
        if (threadsPtr != nullptr) WMDisposePtr ((Ptr)threadsPtr);
        if (paramArrayPtr != nullptr) WMDisposePtr ((Ptr)paramArrayPtr);
        p -> result = (double)(result - FIRST_XOP_ERR);
        return result;
    }
    // fill param array struct
    for (iThread = 0; iThread < nThreads; iThread++){
        paramArrayPtr[iThread].inPutWaveType = inPutWaveType;
        paramArrayPtr[iThread].inPutDataStartPtr = srcWaveStart;
        paramArrayPtr[iThread].outPutDataStartPtr = destWaveStart;
        paramArrayPtr[iThread].flatDim = 2;     // for Z
        paramArrayPtr[iThread].projMode = 0;    // doesn't matter for a single frame
        paramArrayPtr[iThread].xSize = inPutDimensionSizes [0];
        paramArrayPtr[iThread].ySize = inPutDimensionSizes [1];
        paramArrayPtr[iThread].zSize =inPutDimensionSizes [2];
        paramArrayPtr[iThread].startP = slice;
        paramArrayPtr[iThread].endP= slice;
        paramArrayPtr[iThread].ti=iThread; // number of this thread, starting from 0
        paramArrayPtr[iThread].tN =nThreads; // total number of threads
    }
    // create the threads
    for (iThread = 0; iThread < nThreads; iThread++){
        pthread_create (&threadsPtr[iThread], NULL, ProjectThread, (void *) &paramArrayPtr[iThread]);
    }
    //join the threads
    for (iThread = 0; iThread < nThreads; iThread++){
        pthread_join (threadsPtr[iThread], NULL);
    }
    WMDisposePtr ((Ptr)threadsPtr);     // free memory for pThreads Array
    WMDisposePtr ((Ptr)paramArrayPtr);   // Free paramaterArray memory
    WaveHandleModified(outPutWaveH);     // Inform Igor that we have changed the output wave.
    p->result= (0);
    return (0);
}
