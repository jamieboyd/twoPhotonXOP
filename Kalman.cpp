#include "twoPhoton.h"

/* ---------------------------------------------Kalman.cpp ----------------------------------------------------------
 Code for Kalman averaging of frames in a 3d wave, 2d wave, or a list of waves
 Last Modified 2026/09/30 by Jamie Boyd - getting native windows threading working for KalmanAllFrames, KalmanSpecFrames, KalmanListFrames
 Modified 2026/09/29 by Jamie Boyd - got native windows threading working, cleaned up a bit
 Modified 2026/08/11 by Jamie Boyd - changed long to SInt32 and unsigned long to SInt32 and added support for 64 bitInteger waves
 Modified 2026/08/02 by Jamie Boyd - start at using native windows threading instead of pThreadsWin32 - not working yet, hence the "dirty" includes
 ---------------------------------------------------------------------------------------------------------------*/

/* ************************************************** KalmanT **************************************************************
 The following template is used to handle any one of the 8 types of wave data, for any of the three Kalman averaging functions
 that work on 3D waves (KalmanAllFrames, KalmanSpecFrames, KalmanWaveToFrame)
 Last modified 2026/09/29 by Jamie Boyd */
template <typename T> int KalmanT(T *srcWaveStart, T *destWaveStart, CountInt pixPerThread, CountInt pixToNextFrame, CountInt numLayers, float multiplier) {
    T* destWaveEnd = destWaveStart + pixPerThread;    // End of the output layer we are putting the average into
    T* srcWave, *destWave;    // pointers used to iterate through source wave and destination wave
    /* If multiplier is < 1, use standard averaging across layers with a floating point temporary value */
    if (multiplier < 1){
        double tempVal;
        CountInt pixPerLayer = pixPerThread + pixToNextFrame;
        T* srcWaveEnd = srcWaveStart + (numLayers * pixPerLayer);
        CountInt toNextSrcPix = (numLayers * pixPerLayer) - 1;
        for (srcWave = srcWaveStart, destWave = destWaveStart ; destWave < destWaveEnd ; srcWave -= toNextSrcPix, destWave ++){
            for (tempVal =0; srcWave < srcWaveEnd; srcWave += pixPerLayer)
            tempVal += *srcWave;
            *destWave = tempVal/numLayers;
        }
    }else{
        CountInt layer; // used to iterate through layers
        // Special stuff for first frame
        if (srcWaveStart == destWaveStart){    //this happens when collapsing a wave into the first frame
            if (multiplier > 1){
                //Multiply the first layer by Multiplier
                for (destWave = destWaveStart; destWave < destWaveEnd; destWave++){
                    *destWave *= multiplier;
                }
            }
            // position src wave pointer at start of 2nd frame
            srcWave = srcWaveStart + pixPerThread + pixToNextFrame;
        }else{ //srcwave and dsetwave are different
            if (multiplier > 1){ // To increase precision in averaging in special instances where, for example, 16 bit waves contain less than 16 bits of data
                // Set output layer = first layer of input wave * Multiplier
                for (destWave = destWaveStart, srcWave = srcWaveStart; destWave < destWaveEnd; destWave++, srcWave++){
                    *destWave = *srcWave * multiplier;
                }
            }else{ // first frame when No Multiplier, and Src and dest are different
                for (destWave = destWaveStart, srcWave = srcWaveStart; destWave < destWaveEnd; destWave++, srcWave++){
                    *destWave = *srcWave;
                }
            }
            // advance src wave pointer to start of next frame
            srcWave += pixToNextFrame;
        }
        //For each remaining layer in input wave, iterate through, averaging the input value into the output layer
        if (multiplier > 1){
            for(layer=1; layer < numLayers; layer++, srcWave += pixToNextFrame) {
                for (destWave = destWaveStart; destWave < destWaveEnd; destWave++, srcWave++){
                    *destWave = ((*destWave * layer) + *srcWave * multiplier)/(layer + 1);
                }
            }
            // Divide output layer by Multiplier
            for (destWave = destWaveStart;destWave < destWaveEnd; destWave++){
                *destWave /= multiplier;
            }
        }else{ // no multiplier
            for(layer=1; layer < numLayers; layer++, srcWave += pixToNextFrame) {
                for (destWave = destWaveStart; destWave < destWaveEnd; destWave++, srcWave++){
                    *destWave = ((*destWave * layer) + *srcWave)/(layer + 1);
                }
            }
        }
    }
    return 0;
}

/* *************************************************** KalmanThreadParams *************************************************************
 Structure to pass data to each KalmanThread for functions that work on 3D waves (KalmanAllFrames, KalmanSpecFrames, KalmanWaveToFrame)
 Last modified 2026/09/29 by Jamie Boyd */
typedef struct KalmanThreadParams{
    int inPutWaveType;          // Wavemetrics code for data type of wave
    char* inPutDataStartPtr;    // byte pointer to start of data in input wave - cast to pointer for input wave type
    CountInt startLayer;        // number of layer to start with Kalman Averaging
    CountInt numLayers;         // number of layers to Kalman average, starting from from start layer
    char* outPutDataStartPtr;   // byte pointer to start of data in output wave - cast to pointer for input wave type
    CountInt outPutLayer;       // number of layer in output wave where result of averging is put
    float multiplier;           // multiply data by this value when averaging to maintain precision (if data is small rel. to type)
    CountInt layerSize;         // number of pixels in a single layer (numCols * numRows)
    CountInt threadOffset;      // offset in pixels from start of start layer for this thread
    CountInt threadPix;         // number of pixels for this thread to do
 } KalmanThreadParams, *KalmanThreadParamsPtr;


/* ********************************* KalmanThread *******************************************************************************
 Each thread to do Kalman averaging for (KalmanAllFrames, KalmanSpecFrames, KalmanWaveToFrame) starts with this function
 Last Modified 2013/07/16 by Jamie Boyd */
#ifdef __GNUC__
void* KalmanThread (void* threadarg){
#endif
#ifdef _WINDOWS_
DWORD WINAPI KalmanThread(LPVOID threadarg) {
#endif
    KalmanThreadParams *p = (KalmanThreadParams*) threadarg;
    CountInt startLayer = p->startLayer;
    CountInt numLayers = p->numLayers;
    CountInt outPutLayer = p->outPutLayer;
    float multiplier = p->multiplier;
    CountInt layerSize = p->layerSize;
    CountInt threadOffset = p->threadOffset;
    CountInt threadPix = p->threadPix;
    CountInt pixToNextFrame = layerSize - threadPix;
    switch (p->inPutWaveType) {
    case NT_I8:
        KalmanT ((char*)p->inPutDataStartPtr + (startLayer * layerSize) + threadOffset, (char*)p->outPutDataStartPtr + (outPutLayer * layerSize) + threadOffset, threadPix, pixToNextFrame, numLayers, multiplier);
        break;
    case (NT_I8 | NT_UNSIGNED):
       KalmanT ((unsigned char*)p->inPutDataStartPtr + (startLayer * layerSize) + threadOffset, (unsigned char*)p->outPutDataStartPtr + (outPutLayer * layerSize) + threadOffset, threadPix, pixToNextFrame, numLayers, multiplier);
        break;
    case NT_I16:
        KalmanT ((short*)p->inPutDataStartPtr + (startLayer * layerSize) + threadOffset, (short*) p->outPutDataStartPtr + (outPutLayer * layerSize) + threadOffset, threadPix, pixToNextFrame, numLayers, multiplier);
        break;
    case (NT_I16 | NT_UNSIGNED):
        KalmanT ((unsigned short*)p->inPutDataStartPtr + (startLayer * layerSize) + threadOffset, (unsigned short*) p->outPutDataStartPtr + (outPutLayer * layerSize) + threadOffset, threadPix, pixToNextFrame, numLayers, multiplier);
        break;
    case NT_I32:
        KalmanT ((SInt32*)p->inPutDataStartPtr + (startLayer * layerSize) + threadOffset, (SInt32*) p->outPutDataStartPtr + (outPutLayer * layerSize) + threadOffset, threadPix, pixToNextFrame, numLayers, multiplier);
        break;
    case (NT_I32| NT_UNSIGNED):
        KalmanT ((UInt32*)p->inPutDataStartPtr + (startLayer * layerSize)+ threadOffset, (UInt32*) p->outPutDataStartPtr + (outPutLayer * layerSize) + threadOffset, threadPix, pixToNextFrame, numLayers, multiplier);
        break;
    case NT_I64:
        KalmanT((SInt64*)p->inPutDataStartPtr + (startLayer * layerSize) + threadOffset, (SInt64*)p->outPutDataStartPtr + (outPutLayer * layerSize) + threadOffset, threadPix, pixToNextFrame, numLayers, multiplier);
        break;
    case (NT_I64 | NT_UNSIGNED):
       KalmanT((UInt64*)p->inPutDataStartPtr + (startLayer * layerSize) + threadOffset, (UInt64*)p->outPutDataStartPtr + (outPutLayer * layerSize) + threadOffset, threadPix, pixToNextFrame, numLayers, multiplier);
        break;
    case NT_FP32:
        KalmanT ((float*)p->inPutDataStartPtr + (startLayer * layerSize) + threadOffset, (float*) p->outPutDataStartPtr + (outPutLayer * layerSize) + threadOffset, threadPix, pixToNextFrame, numLayers, multiplier);
        break;
    case NT_FP64:
        KalmanT ((double*)p->inPutDataStartPtr + (startLayer * layerSize) + threadOffset, (double*) p->outPutDataStartPtr + (outPutLayer * layerSize) + threadOffset, threadPix, pixToNextFrame, numLayers, multiplier);
        break;
    }
#ifdef __GNUC__
    return nullptr;
#endif
#ifdef _WINDOWS_
    return 0;
#endif
}

/* ******************************************** KalmanAllFrames ********************************************************************
 KalmanAllFrames XOP entry function
 Does Kalman averaging across all layers in a 3D wave and places results in a new 2D wave that it makes from provided name
 KalmanAllFramesParams
 waveHndl inPutWaveH    handle to a 3D input wave
 Handle outPutPath      A handle to a string containing path to output wave we want to make
 double multiplier      Multiplier for,e.g., 16 bit waves containing less than 16 bits of data
 double overWrite       0 to give errors when wave already exists. non-zero to overwrite existing wave without warning.
 result                 0 for success, error code for failure
Last Modified 2026/09/30 by Jamie Boyd  */
extern "C" int KalmanAllFrames(KalmanAllFramesParamsPtr p) {
    int result = 0;                 // The error returned from various Wavemetrics functions
    waveHndl inPutWaveH = nullptr;  // Handle to the input wave
    waveHndl outPutWaveH = nullptr; // Handle to the output wave that we will make
    int inPutWaveType;              // Wavemetrics numeric code for data type of wave, input and output must be the same
    int inPutDimensions;            // The number of dimensions used in the input wave, must be 3
    CountInt inPutOffset;           // offset in bytes from begnning of handle to the actual data in input wave
    CountInt outPutOffset;          // offset in bytes from begnning of handle to the actual data in
    CountInt inputDimSizes[MAX_DIMENSIONS+1];
                                    // An array used to hold the sizes of each dimension of the input wave (rows, columns, layers)
    UInt16 outPutPathLen;           //Length of the path to the target folder (output path:wave name)
    DataFolderHandle outPutDFHandle;// Handle to the datafolder where we will put the output wave
    DFPATH outPutPath;              // string to hold data folder path of output wave
    WVNAME outPutWaveName;          // C-style string to hold name of output wave
    char *inPutDataStartPtr;        // pointer to where data starts, past header, di units, etc. in input Wave
    char *outPutDataStartPtr;       // pointer to where data starts, past header, di units, etc. in output Wave
    float multiplier = (float)(p->multiplier);
                                    // multiplier useful for integer waves containing less than their full range of data
    UInt8 overWriteOK = (UInt8)p->overWrite;
                                    // pass 0 to not overwrite output wave if it already exists, 1 to overwrite old waves
    UInt8 isOverWritingInput;       // will be set if output wave is overwriting the input wave (we will flatten the wave to 2D)
    CountInt numLayers;             // number of layers in the input wave
    CountInt layerSize;             // number of points in a single image frame (rows x columns)
   // variables for threading
    UInt8 nThreads;                 // number of threads used for processing
    UInt8 iThread;                  // variable used to iterate through threads
    CountInt threadPix;             // number of pixels to be done by a thread
    KalmanThreadParamsPtr paramArrayPtr = nullptr;
                                    // Pointer to an array of KalmanThreadParams structures
#ifdef __GNUC__
    pthread_t* threadsPtr = nullptr;   // pointer to an array of pThread_t structures we will make
#endif
#ifdef  _WINDOWS_
    HANDLE *threadsPtr = nullptr;      // pointer to an array of Windows handles, which will point to thread structures we make
    DWORD *threadIDsPtr = nullptr;     // pointer to an array of DWORD (32-bit unsigned integer) that Windows uses for thread IDs
#endif
    try {
        // Get handle to input wave.
        inPutWaveH = p ->inPutWaveH;
        if(inPutWaveH == nullptr) throw result = NON_EXISTENT_WAVE;
        // get wave data type and check that we don't have a text wave
        inPutWaveType = WaveType(inPutWaveH);
        if (inPutWaveType==TEXT_WAVE_TYPE) throw result = NOTEXTWAVES;
        //Get number of used dimensions in input wave.
        if (MDGetWaveDimensions(inPutWaveH, &inPutDimensions, inputDimSizes))throw result = WAVEERROR_NOS;
        // Check that input wave is 3D
        if (inPutDimensions != 3) throw result = INPUTNEEDS_3D_WAVE;
        // Save number of layers as we will be resizing dimensions array to re-use it for making 2D wave
        numLayers = inputDimSizes [LAYERS];
        // If outPutPath is empty string, we are overwriting existing wave, which may or may  not be OK
        outPutPathLen = WMGetHandleSize (p->outPutPath);
        if (outPutPathLen == 0){
            if (!(overWriteOK)) throw result = OVERWRITEALERT;
            outPutWaveH = inPutWaveH;
            isOverWritingInput = 1;
        }else{ // we have a passed in string we need to parse
            isOverWritingInput = 0;
            // Parse passed in outPut path into dataFolder path and waveName
            ParseWavePath (p->outPutPath, outPutPath, outPutWaveName);
            // Get a handle to the output datafolder
            if (GetNamedDataFolder (NULL, outPutPath, &outPutDFHandle))throw result = WAVEERROR_NOS;
            // get a handle to the output wave, if it exists
            outPutWaveH = FetchWaveFromDataFolder(outPutDFHandle, outPutWaveName);
            // if outPutwave does not exist, make it a 2D wave, same frame size and wavetype as input
            if (outPutWaveH != nullptr){
                if (outPutWaveH == inPutWaveH){
                    isOverWritingInput = 1;
                }
                if (!(overWriteOK)) throw result = OVERWRITEALERT; //check if overwriting existing wave is OK
            }
            // if not overwritingInput, make the output wave, overwriting it if it exists
            if (!(isOverWritingInput)){
                inputDimSizes [LAYERS] = 0;     // now inputDimSizes descibes a 2D wave with same columns and layers as input wave
                //No liberal wave names for output wave, maybe we don't really need to enforce non-liberal names?
                CleanupName (0, outPutWaveName, MAX_OBJ_NAME);
                if (MDMakeWave (&outPutWaveH, outPutWaveName, outPutDFHandle, inputDimSizes, inPutWaveType, overWriteOK)) throw result = WAVEERROR_NOS;
            }
        } // we now have handles to input wave and output wave
        //Get data offsets for the 2 waves (1 wave, if overwriting)
        if (MDAccessNumericWaveData(inPutWaveH, kMDWaveAccessMode0, &inPutOffset) != 0) throw result = WAVEERROR_NOS;
        inPutDataStartPtr = (char*)(*inPutWaveH) + inPutOffset;
        if (isOverWritingInput){
            outPutOffset = inPutOffset;
            outPutDataStartPtr = inPutDataStartPtr;
        }else{
            if (MDAccessNumericWaveData(outPutWaveH, kMDWaveAccessMode0, &outPutOffset)) throw result = WAVEERROR_NOS;
            outPutDataStartPtr =  (char*)(*outPutWaveH) + outPutOffset;
        }
        // make threads for each processor core, but make sure each thread gets at least one pixel
        layerSize = inputDimSizes [COLUMNS] * inputDimSizes [ROWS];
        if (layerSize < gNumProcessors){
            nThreads = layerSize;
        }else{
            nThreads = gNumProcessors;
        }
        threadPix = layerSize/nThreads;     // number of pixels done by a thread
        // make array of parameter structures pointed to by paramArrayPtr
        paramArrayPtr = (KalmanThreadParamsPtr)WMNewPtr (nThreads * sizeof(KalmanThreadParams));
        if (paramArrayPtr == nullptr) throw result = MEMFAIL;
#ifdef __GNUC__
        // set threadsPtr to an array of pthread_t for MacOS
        threadsPtr =(pthread_t*)WMNewPtr(nThreads * sizeof(pthread_t));
        if (threadsPtr == nullptr) throw result = MEMFAIL;
#endif
#ifdef _WINDOWS_
        // set threadsPtr to an array of Handles for Windows
        threadsPtr = (HANDLE*)WMNewPtr(nThreads * sizeof(HANDLE));
        // also set threadIDsPtr to an array of DWORDs for thread IDs
        threadIDsPtr = (DWORD*)WMNewPtr(nThreads * sizeof(DWORD));
        if ((threadsPtr == nullptr) || (threadIDsPtr == nullptr)) throw result = MEMFAIL;
#endif
        // catch error before starting threads - threads don't return errors
    }catch (int result){
        if (paramArrayPtr != nullptr) WMDisposePtr ((Ptr)paramArrayPtr);
        if (threadsPtr != nullptr) WMDisposePtr ((Ptr)threadsPtr);
#ifdef _WINDOWS_
        if (threadIDsPtr != nullptr) WMDisposePtr((Ptr)threadIDsPtr);
#endif
        WMDisposeHandle(p->outPutPath);  // dispose passed in string paramater
        p -> result = (double)(result - FIRST_XOP_ERR);
#ifdef NO_IGOR_ERR
        return (0);
#else
        return (result);
#endif
    } // Now we have both wavehandles and all structures we need created
    // fill the array of KalmanThreadParams
    for (iThread = 0; iThread < nThreads; iThread++){
        paramArrayPtr[iThread].inPutWaveType = inPutWaveType;
        paramArrayPtr[iThread].inPutDataStartPtr = inPutDataStartPtr;
        paramArrayPtr[iThread].startLayer = 0;
        paramArrayPtr[iThread].numLayers = numLayers;
        paramArrayPtr[iThread].outPutDataStartPtr = outPutDataStartPtr;
        paramArrayPtr[iThread].outPutLayer = 0;
        paramArrayPtr[iThread].multiplier = multiplier;
        paramArrayPtr[iThread].layerSize = layerSize;
        paramArrayPtr[iThread].threadOffset = threadPix * iThread;
        paramArrayPtr[iThread].threadPix = threadPix;
     }
    // last thread gets any left-over pixels
    paramArrayPtr[nThreads -1].threadPix += layerSize % nThreads;
#ifdef __GNUC__
    // create the threads
    for (iThread = 0; iThread < nThreads; iThread++){
        pthread_create(&threadsPtr[iThread], NULL, KalmanThread, (void*)&paramArrayPtr[iThread]);
    }
    // wait for threads to finish
    for (iThread = 0; iThread < nThreads; iThread++) {
        pthread_join(threadsPtr[iThread], NULL);
    }
#endif
#ifdef _WINDOWS_
    // create the threads with Default security attributes, Default stack size, KalmanThread function, a pointer to a parameter structure, Default creation flags, and a pointer to a threadId
    for (iThread = 0; iThread < nThreads; iThread++){
        threadsPtr[iThread] = CreateThread( NULL, 0, (LPTHREAD_START_ROUTINE)KalmanThread, (LPVOID)&paramArrayPtr[iThread], 0, &threadIDsPtr[iThread]);
    }
    // wait for the threads to finish
    WaitForMultipleObjects(nThreads, threadsPtr, TRUE, INFINITE);
    // Close thread handles
    for (iThread = 0; iThread < nThreads; iThread++) {
        CloseHandle(threadsPtr[iThread]);
    }
    // free threadIDs, only used for windows
    WMDisposePtr((Ptr)threadIDsPtr);    // freee threadIDs
#endif
    // free thread structures, same for both Windows and MacOS
     WMDisposePtr((Ptr)threadsPtr);     // free memory for threads
     WMDisposePtr((Ptr)paramArrayPtr);  // Free paramaterArray memory
    // collapse 3D wave into 2D, if overwriting input wave
    if (isOverWritingInput){    //then collapsing a 3D wave to 2 D
        inputDimSizes [0] = -1;
        inputDimSizes [1] = -1;
        inputDimSizes [2] = 0;
        inputDimSizes [3] = 0;
        MDChangeWave (outPutWaveH, -1, inputDimSizes);
    }
    // dispose passed in string parameter
    WMDisposeHandle(p->outPutPath);
    // Inform Igor that we have changed output wave
    WaveHandleModified(outPutWaveH);
    p -> result = (0);
    return (0);
}

/* ******************************************* KalmanSpecFrames *********************************************************************
 KalmanSpecFrames XOP entry function
 Averages a specified range of layers of the input wave into a specified layer of the output wave
 KalmanSpecFramesParams
 inPutWaveH          handle to input wave
 startLayer          start of layers to average for input wave
 endLayer            end of layers to average
 outPutWaveH         handle to output wave
 outPutLayer         layer of output wave to receive results of averaging
 multiplier          Multiplier, as for 16 bit waves containing less than 16 bits of data
 result              0 or error code
 Last Modified 2026/09/29 by Jamie Boyd */
extern "C" int KalmanSpecFrames(KalmanSpecFramesParamsPtr p) {
    int result = 0;                     // The error returned from various Wavemetrics functions
    waveHndl inPutWaveH = nullptr;          // Handle to the input wave
    waveHndl outPutWaveH = nullptr;         // Handle to the output wave
    int inPutWaveType;                  // Wavemetrics numeric code for data type of wave, for input wave
    int outPutWaveType;                 // Wavemetrics numeric code for data type of wave, for output wave
    int inPutDimensions;                // number of dimensions in input wave
    int outPutDimensions;               // number of dimensions in output wave
    CountInt inputDimSizes[MAX_DIMENSIONS+1];
                                        // An array used to hold the sizes of each dimension of the input wave
    CountInt outPutDimSizes[MAX_DIMENSIONS+1];
                                        // An array used to hold the sizes of each dimension of the output wave
    CountInt inPutOffset;               // offset in bytes from begnning of handle to the actual data in inputWave
    CountInt outPutOffset;              // offset in bytes from begnning of handle to the actual data in outputWave
    CountInt startLayer, endLayer, layersToDo;  //vaiables for iterating through the data.
    CountInt outPutLayer;               //The layer of the output wave that gets the result
    char *inPutDataStartPtr;            // pointer to start of data in inputwave
    char *outPutDataStartPtr;           // pointer to start of data in outputwave
    float multiplier = (float)(p->multiplier);    // multiplier for integer waves containing less than their full range of data
    UInt8 iThread, nThreads;
    CountInt layerSize;
    CountInt threadPix;
    KalmanThreadParamsPtr paramArrayPtr = nullptr;
#ifdef __GNUC__
    pthread_t* threadsPtr = nullptr;
#endif
#ifdef _WINDOWS_
    HANDLE* threadsPtr = nullptr;
    DWORD *threadIDsPtr = nullptr;
#endif
    try {
        // Get handles to input wave and kernel. Make sure both waves exist.
        inPutWaveH = p ->inPutWaveH;
        outPutWaveH =  p->outPutWaveH;
        if ((inPutWaveH == nullptr) || (outPutWaveH == nullptr)) throw result = NON_EXISTENT_WAVE;
        // get wave data types and check that datatypes are the same and that neither is a text wave
        inPutWaveType = WaveType(inPutWaveH);
        outPutWaveType = WaveType(outPutWaveH);
        if (inPutWaveType != outPutWaveType) throw result = NOTSAMEWAVETYPE;
        if ((inPutWaveType==TEXT_WAVE_TYPE) || (outPutWaveType==TEXT_WAVE_TYPE))throw result = NOTEXTWAVES;
        // Get number of used dimensions in waves.
        if (MDGetWaveDimensions(inPutWaveH, &inPutDimensions, inputDimSizes))throw result=WAVEERROR_NOS;
        if (MDGetWaveDimensions(outPutWaveH, &outPutDimensions, outPutDimSizes))throw result=WAVEERROR_NOS;
        // Check that input wave is 3D and output wave is 2D or 3D
        if (inPutDimensions != 3) throw result = INPUTNEEDS_3D_WAVE;
        if (!((outPutDimensions == 2) || (outPutDimensions == 3))) throw result = OUTPUTNEEDS_2D3D_WAVE;
        // Check that X and Y dimensions of the 2 waves are the same size.
        if (!((inputDimSizes[ROWS] == outPutDimSizes [ROWS]) && (inputDimSizes[COLUMNS] == outPutDimSizes [COLUMNS]))) throw result = NOTSAMEDIMSIZE;
        // outPut layer must be 0 to use 2 D wave as output wave
        outPutLayer = (CountInt)p ->outPutLayer;
        if (((outPutLayer != 0) && (outPutLayer > outPutDimSizes [2] - 1)) || (outPutLayer < 0)) throw result = INVALIDOUTPUTFRAME;
        startLayer = p -> startLayer;
        endLayer = p->endLayer;
        // swap startLayer and endLayer, if they are reversed
        if (endLayer < startLayer) {
            CountInt temp;
            SWAP(startLayer, endLayer);
        }
        // clip start Layer to 0, if negative
        if (startLayer < 0) startLayer = 0;
         // Clip endlayer to the last layer of the input wave
        if (endLayer > inputDimSizes [2] -1) endLayer = inputDimSizes [2] -1;
        // Calculate number of layers to do
        layersToDo = endLayer - startLayer + 1;
        // As X and Y are same size for input and output waves, we need only look at the input wave to get points per layer
        layerSize = inputDimSizes[ROWS] * inputDimSizes[COLUMNS];
        //Get data offsets for the 2 waves
        if (MDAccessNumericWaveData(inPutWaveH, kMDWaveAccessMode0, &inPutOffset)) throw result = WAVEERROR_NOS;
        inPutDataStartPtr = (char*)(*inPutWaveH) + inPutOffset;
        if (MDAccessNumericWaveData(outPutWaveH, kMDWaveAccessMode0, &outPutOffset))throw result = WAVEERROR_NOS;
        outPutDataStartPtr =  (char*)(*outPutWaveH) + outPutOffset;
        // multiprocessor initialization
        layerSize = inputDimSizes [COLUMNS] * inputDimSizes [ROWS];
        if (layerSize < gNumProcessors){
            nThreads = layerSize;
        }else{
            nThreads = gNumProcessors;
        }
        threadPix = layerSize/nThreads;     // number of pixels done by a thread, note integer truncation
        // make an array of parameter structures
        paramArrayPtr = (KalmanThreadParamsPtr)WMNewPtr (nThreads * sizeof(KalmanThreadParams));
        if (paramArrayPtr == nullptr) throw result = MEMFAIL;
#ifdef __GNUC__
        // make array for pthread structures
        threadsPtr =(pthread_t*)WMNewPtr(nThreads * sizeof(pthread_t));
        if (threadsPtr == nullptr) throw result = MEMFAIL;
#endif
#ifdef _WINDOWS_
        // make an array of handles for threads
        threadsPtr = (HANDLE*)WMNewPtr(nThreads * sizeof(HANDLE));
        // also make array of DWORD for threadIDs
        threadIDsPtr = (DWORD*)WMNewPtr(nThreads * sizeof(DWORD));
        if ((threadsPtr == nullptr) || (threadIDsPtr == nullptr)) throw result = MEMFAIL;
#endif
    }catch (int result){
        if (paramArrayPtr != nullptr) WMDisposePtr ((Ptr)paramArrayPtr);
        if (threadsPtr != nullptr) WMDisposePtr ((Ptr)threadsPtr);
#ifdef _WINDOWS_
        if (threadIDsPtr != nullptr) WMDisposePtr((Ptr)threadIDsPtr);
#endif
        p -> result = (double)(result - FIRST_XOP_ERR);
#ifdef NO_IGOR_ERR
    return (0);
#else
    return (result);
#endif
    }
    // Now we have all wavehandles and all the memory we need allocated
    // fill the paramater array for each thread
    for (iThread = 0; iThread < nThreads; iThread++){
        paramArrayPtr[iThread].inPutWaveType = inPutWaveType;
        paramArrayPtr[iThread].inPutDataStartPtr = inPutDataStartPtr;
        paramArrayPtr[iThread].startLayer = startLayer;
        paramArrayPtr[iThread].numLayers= endLayer - startLayer + 1;
        paramArrayPtr[iThread].outPutDataStartPtr = outPutDataStartPtr;
        paramArrayPtr[iThread].outPutLayer = outPutLayer;
        paramArrayPtr[iThread].multiplier = multiplier;
        paramArrayPtr[iThread].layerSize = layerSize;
        paramArrayPtr[iThread].threadOffset = threadPix * iThread;
        paramArrayPtr[iThread].threadPix = threadPix;
    }
    // last thread gets any left-over pixels
    paramArrayPtr[nThreads -1].threadPix += layerSize % nThreads;
#ifdef __GNUC__
    // create the threads
    for (iThread = 0; iThread < nThreads; iThread++){
        pthread_create (&threadsPtr[iThread], NULL, KalmanThread, (void *) &paramArrayPtr[iThread]);
    }
    // Wait till all the threads are finished
    for (iThread = 0; iThread < nThreads; iThread++){
        pthread_join (threadsPtr[iThread], NULL);
    }
#endif
#ifdef _WINDOWS_
    for (iThread = 0; iThread < nThreads; iThread++){
        // create the threads with Default security attributes, Default Stack size, KalmanThread function, pointer to a parameter array, Default creation flags, and a pointer to a threadId
        threadsPtr[iThread] = CreateThread(NULL, 0, (LPTHREAD_START_ROUTINE)KalmanThread, (LPVOID)&paramArrayPtr[iThread], 0, &threadIDsPtr[iThread]);
    }
    // wait for the threads to finish
    WaitForMultipleObjects(nThreads, threadsPtr, TRUE, INFINITE);
    // Close thread handles
    for (iThread = 0; iThread < nThreads; iThread++) {
        CloseHandle(threadsPtr[iThread]);
    }
    // free thread ID array, only used on Windows
    WMDisposePtr((Ptr)threadIDsPtr);    // freee thread IDs
#endif
    // free thread structures
    WMDisposePtr((Ptr)threadsPtr);     // free memory for thread pointers Array
    WMDisposePtr((Ptr)paramArrayPtr);  // Free paramaterArray memory
    // Inform Igor that we have changed the output wave.
    WaveHandleModified(outPutWaveH);
    p -> result = (0);
    return (0);
}

/* ************************************************ KalmanWaveToFrame ****************************************************************************
 KalmanWaveToFrame XOP entry function
 Collapses a 3D input wave into a single 2D frame. You can get the same result with KalmanAllFrames by
 using "" as outPut String
 KalmanWaveToFrameParams
 inPutWaveH     handle to input wave
 multiplier     Multiplier, as for  for 16 bit waves containing less than 16 bits of data
 result         0 or error code
 Last Modified: 2026/09/29 by Jamie Boyd */
int KalmanWaveToFrame (KalmanWaveToFrameParamsPtr p) {
    int result=0;                   // The error returned from various Wavemetrics functions
    waveHndl inPutWaveH = nullptr;      // handle to the input wave
    int inPutWaveType;              //  Wavetypes numeric codes for things like 32 bit floating point, 16 bit int, etc
    int inPutDimensions;            // number of dimensions in input wave
    CountInt inputDimSizes[MAX_DIMENSIONS+1];
                                    // an array used to hold the width, height, layers, and chunk sizes
    CountInt numLayers;
    CountInt inPutOffset;           //offset in bytes from begnning of handle to a wave to the actual data - size of headers, etc.
    char *inPutDataStartPtr;        // char pointer to start of data in input wave, after the header
    float multiplier = (float)(p->multiplier); // multiplier for integer waves containing less than their full range of data
    UInt8 iThread, nThreads;
    CountInt layerSize;
    CountInt threadPix;
    KalmanThreadParamsPtr paramArrayPtr = nullptr;
    //**** declare pthreads or HANDLE pointer
    #ifdef __GNUC__
    pthread_t *threadsPtr = nullptr;
    #endif
    #ifdef  _WINDOWS_
    HANDLE* threadsPtr = nullptr;       // WIndows Handles to point to thread structures
    DWORD* threadIDsPtr = nullptr;     // pointer to an array of DWORD (32-bit unsigned integer) that Windows uses for thread IDs
    #endif
    try {
        // Get handle to input wave. Make sure input wave exists.
        inPutWaveH = p->inPutWaveH;
        if(inPutWaveH == nullptr) throw result = NON_EXISTENT_WAVE;
        // get wave data type and check that we don't have a text wave
        inPutWaveType = WaveType(inPutWaveH);
        if (inPutWaveType == TEXT_WAVE_TYPE) throw result = NOTEXTWAVES;
        //Get number of used dimensions in input wave.
        if (MDGetWaveDimensions(inPutWaveH, &inPutDimensions, inputDimSizes))throw result = WAVEERROR_NOS;
        // Check that input wave is 3D
        if (inPutDimensions != 3) throw result = INPUTNEEDS_3D_WAVE;
        //Get data offset for the wave and  make pointer to start of data
        if (MDAccessNumericWaveData(inPutWaveH, kMDWaveAccessMode0, &inPutOffset)) throw result = WAVEERROR_NOS;
        inPutDataStartPtr = (char*)(*inPutWaveH) + inPutOffset;
        // Save z Size as we will be resizing dimensions array
        numLayers = inputDimSizes [LAYERS];
        layerSize = inputDimSizes [COLUMNS] * inputDimSizes [ROWS];
        if (layerSize < gNumProcessors){
            nThreads = layerSize;
        }else{
            nThreads = gNumProcessors;
        }
        threadPix = layerSize/nThreads;     // number of pixels done by a thread, note integer truncation
        // make an array of parameter structures
        paramArrayPtr = (KalmanThreadParamsPtr)WMNewPtr (nThreads * sizeof(KalmanThreadParams));
        if (paramArrayPtr == nullptr) throw result = MEMFAIL;
        #ifdef __GNUC__
        // make array of threads
        threadsPtr =(pthread_t*)WMNewPtr(nThreads * sizeof(pthread_t));
        #endif
        #ifdef _WINDOWS_
        threadsPtr = (HANDLE*)WMNewPtr(nThreads * sizeof(HANDLE));
        threadIDsPtr = (DWORD*)WMNewPtr(nThreads * sizeof(DWORD));
        #endif
        if (threadsPtr == nullptr) throw result = MEMFAIL;
    }catch (int result){
        if (paramArrayPtr != nullptr)  WMDisposePtr ((Ptr)paramArrayPtr);
        if (threadsPtr != nullptr) WMDisposePtr ((Ptr)threadsPtr);
        p -> result = (double)(result - FIRST_XOP_ERR);
#ifdef NO_IGOR_ERR
    return (0);
#else
    return (result);
#endif
    }
    // fill pramaters array
    for (iThread = 0; iThread < nThreads; iThread++){
        paramArrayPtr[iThread].inPutWaveType = inPutWaveType;
        paramArrayPtr[iThread].inPutDataStartPtr = inPutDataStartPtr;
        paramArrayPtr[iThread].startLayer = 0;
        paramArrayPtr[iThread].numLayers = numLayers;
        paramArrayPtr[iThread].outPutDataStartPtr = inPutDataStartPtr;
        paramArrayPtr[iThread].outPutLayer = 0;
        paramArrayPtr[iThread].multiplier = multiplier;
        paramArrayPtr[iThread].layerSize = layerSize;
        paramArrayPtr[iThread].threadOffset = threadPix * iThread;
        paramArrayPtr[iThread].threadPix = threadPix;
    }
    // last thread gets any left-over pixels
    paramArrayPtr[nThreads -1].threadPix += layerSize % nThreads;
#ifdef __GNUC__
    for (iThread = 0; iThread < nThreads; iThread++){
        pthread_create (&threadsPtr[iThread], NULL, KalmanThread, (void *) &paramArrayPtr[iThread]);
    }
    // Wait till all the threads are finished
    for (iThread = 0; iThread < nThreads; iThread++){
        pthread_join (threadsPtr[iThread], NULL);
    }
#endif
#ifdef _WINDOWS_
    for (iThread = 0; iThread < nThreads; iThread++){
        // create the threads with Default security attributes, Default Stack size, KalmanThread function, pointer to a parameter array, Default creation flags, and a pointer to a threadId
        threadsPtr[iThread] = CreateThread( NULL, 0, (LPTHREAD_START_ROUTINE)KalmanThread, (LPVOID)&paramArrayPtr[iThread], 0, &threadIDsPtr[iThread]);
    }
    // wait for the threads to finish
    WaitForMultipleObjects(nThreads, threadsPtr, TRUE, INFINITE);
    // Close thread handles
    for (iThread = 0; iThread < nThreads; iThread++) {
        CloseHandle(threadsPtr[iThread]);
    }
    WMDisposePtr((Ptr)threadIDsPtr);    // freee threadIDs, only used for Windows
#endif
    // free thread structures
    WMDisposePtr((Ptr)threadsPtr);     // free memory for threads
    WMDisposePtr((Ptr)paramArrayPtr);  // Free paramaterArray memory
    // Redimension wave
    inputDimSizes [0] = -1;
    inputDimSizes [1] = -1;
    inputDimSizes [2] = 0;
    inputDimSizes [3] = 0;
    MDChangeWave (inPutWaveH, -1, inputDimSizes);
    WaveHandleModified(inPutWaveH);     // Inform Igor that we have changed the input wave.
    p -> result = (0);
    return (0);
}

/* ********************************************** KalmanNextT ****************************************************************
Template for doing sequential Kalman averaging, src wave is the new wave,
 dest wave is the old wave already averaged iKal times
 Last Modified 2026/09/29 by Jamie Boyd */
template <typename T> void KalmanNextT (T* srcWaveStart, T* destWaveStart, CountInt nPoints, UInt16 iKal) {
    T* srcWavePtr;
    T* destWavePtr;
    T* destWaveEnd = destWaveStart + nPoints;
    if (iKal == 0){
        for (srcWavePtr = srcWaveStart, destWavePtr = destWaveStart; destWavePtr < destWaveEnd; srcWavePtr++, destWavePtr++)
            *destWavePtr = *srcWavePtr;
    }else{
        for (srcWavePtr = srcWaveStart, destWavePtr = destWaveStart; destWavePtr < destWaveEnd; srcWavePtr++, destWavePtr++)
            *destWavePtr = (*destWavePtr  * iKal +  *srcWavePtr)/(iKal + 1);
    }
}

/* ***************************************************** KalmanNextThreadParams ***************************************************************************
Structure to pass data to each KalmanNext Thread
 Last Modified: 2014/01/28 by Jamie Boyd */
typedef struct KalmanNextThreadParams{
    int inPutWaveType;
    char* inPutDataStartPtr;
    char* outPutDataStartPtr;
    UInt16 iKal;
    CountInt threadOffset;              // offset in pixels from start of start layer for this thread
    CountInt threadPix;                 // number of pixels for this thread to do
} KalmanNextThreadParams, *KalmanNextThreadParamsPtr;


/* ************************************************ KalmanNextThread *********************************************************************
 Each thread to do sequential Kalmaning starts with this function
 Last Modified 2014/01/29 by Jamie Boyd */
#ifdef __GNUC__
void* KalmanNextThread (void* threadarg){
#endif
#ifdef _WINDOWS_
DWORD WINAPI KalmanNextThread(LPVOID threadarg) {
#endif
    struct KalmanNextThreadParams *p;
    p = (struct KalmanNextThreadParams*) threadarg;
    switch (p->inPutWaveType) {
        case NT_I8:
            KalmanNextT ((char*)p->inPutDataStartPtr + p->threadOffset, (char*)p->outPutDataStartPtr + p->threadOffset, p->threadPix, p->iKal);
            break;
        case (NT_I8 | NT_UNSIGNED):
            KalmanNextT ((unsigned char*)p->inPutDataStartPtr + p->threadOffset, (unsigned char*)p->outPutDataStartPtr + p->threadOffset, p->threadPix, p->iKal);
            break;
        case NT_I16:
            KalmanNextT ((short*)p->inPutDataStartPtr + p->threadOffset, (short*)p->outPutDataStartPtr + p->threadOffset, p->threadPix, p->iKal);
            break;
        case (NT_I16 | NT_UNSIGNED):
            KalmanNextT (((unsigned short*)p->inPutDataStartPtr) + p->threadOffset, ((unsigned short*)p->outPutDataStartPtr) + p->threadOffset, p->threadPix, p->iKal);
            break;
        case NT_I32:
            KalmanNextT ((SInt32*)p->inPutDataStartPtr + p->threadOffset, (SInt32*)p->outPutDataStartPtr + p->threadOffset, p->threadPix, p->iKal);
            break;
        case (NT_I32| NT_UNSIGNED):
            KalmanNextT ((UInt32*)p->inPutDataStartPtr + p->threadOffset, (UInt32*)p->outPutDataStartPtr + p->threadOffset, p->threadPix, p->iKal);
            break;
        case NT_I64:
            KalmanNextT((SInt64*)p->inPutDataStartPtr + p->threadOffset, (SInt64*)p->outPutDataStartPtr + p->threadOffset, p->threadPix, p->iKal);
            break;
        case (NT_I64 | NT_UNSIGNED):
            KalmanNextT((UInt64*)p->inPutDataStartPtr + p->threadOffset, (UInt64*)p->outPutDataStartPtr + p->threadOffset, p->threadPix, p->iKal);
            break;
        case NT_FP32:
            KalmanNextT ((float*)p->inPutDataStartPtr + p->threadOffset, (float*)p->outPutDataStartPtr + p->threadOffset, p->threadPix, p->iKal);
            break;
        case NT_FP64:
            KalmanNextT ((double*)p->inPutDataStartPtr + p->threadOffset, (double*)p->outPutDataStartPtr + p->threadOffset, p->threadPix, p->iKal);
            break;
    }
#ifdef __GNUC__
    return nullptr;
#endif
#ifdef _WINDOWS_
    return 0;
#endif
}

/* ********************************************* KalmanNext ***************************************************************************
 KalmanNext XOP entry function
 Averages a wave into an already averaged wave. Both waves must have same data type and same dimensions.
 KalmanNextParams:
 inPutWaveH     handle to input wave
 outPutWaveH    handle to output wave
 iKal           index of wave are we adding
 result         0 or error code
 Last Modified 2026/09/30 by Jamie Boyd */
extern "C" int KalmanNext (KalmanNextParamsPtr p) {
    int result = 0;                                 // The error returned from various Wavemetrics functions
    waveHndl outPutWaveH = nullptr;                 // handle to output wave
    waveHndl inPutWaveH = nullptr;                  // handle to input wave
    int inPutWaveType, outPutWaveType;              // Wavetypes numeric codes for things like 32 bit floating point, 16 bit int, etc
    int inPutDimensions, outPutDimensions;          // number of numDimensions in input and output waves
    CountInt inputDimSizes[MAX_DIMENSIONS+1];       // an array used to hold the width, height, layers, and chunk sizes
    CountInt outPutDimSizes[MAX_DIMENSIONS+1];
    char* outPutDataStartPtr;                       // Pointer to start of output wave
    char* inPutDataStartPtr;                        // Pointer to start of input wave
    CountInt inPutOffset, outPutOffset;             //offset in bytes from begnning of handle to a wave to the actual data - size of headers, units, etc.
    UInt16 iKal = (UInt16)p->iKal;
    UInt8 iThread, nThreads;
    CountInt layerSize;
    CountInt threadPix;
    KalmanNextThreadParamsPtr paramArrayPtr = nullptr;
#ifdef __GNUC__
    pthread_t* threadsPtr = nullptr;
#endif
#ifdef _WINDOWS_
    HANDLE* threadsPtr = nullptr;
    DWORD *threadIDsPtr = nullptr;     // pointer to an array of DWORD (32-bit unsigned integer) that Windows uses for thread IDs
#endif
    try {
        // Get handle to input wave. Make sure input wave exists.
        inPutWaveH = p->inPutWaveH;
        if(inPutWaveH == nullptr)throw result = NON_EXISTENT_WAVE;
        // get wave data type and check that we don't have a text wave
        inPutWaveType = WaveType(inPutWaveH);
        if (inPutWaveType==TEXT_WAVE_TYPE) throw result = NOTEXTWAVES;
        //Get number of used dimensions in input wave.
        if (MDGetWaveDimensions(inPutWaveH, &inPutDimensions, inputDimSizes))throw result = WAVEERROR_NOS;
        if (inPutDimensions != 2) throw result = INPUTNEEDS_2D_WAVE;
        // Get handle to outPut wave. Make sure outPut wave exists.
        outPutWaveH = p->outPutWaveH;
        if(outPutWaveH == nullptr) throw result = NON_EXISTENT_WAVE;
        // get wave data type and check that we don't have a text wave
        outPutWaveType = WaveType(outPutWaveH);
        if (outPutWaveType==TEXT_WAVE_TYPE) throw result = NOTEXTWAVES;
        //Get number of used dimensions in outPut wave.
        if (MDGetWaveDimensions(outPutWaveH, &outPutDimensions, outPutDimSizes))throw result = WAVEERROR_NOS;
        // check that waves are the same types and dimension sizes
        if (inPutWaveType != outPutWaveType) throw result = NOTSAMEDIMSIZE;
        if (inPutDimensions != outPutDimensions) throw result = NOTSAMEDIMSIZE;
        // check sizes of each dimension, and calculate total number of points as well
        if (!((inputDimSizes [ROWS] = outPutDimSizes [ROWS]) && (inputDimSizes [COLUMNS] = outPutDimSizes [COLUMNS]))) throw result = NOTSAMEDIMSIZE;
        //Get data offset for the waves
        if (MDAccessNumericWaveData(inPutWaveH, kMDWaveAccessMode0, &inPutOffset)) throw result = WAVEERROR_NOS;
        inPutDataStartPtr = (char*)(*inPutWaveH) + inPutOffset;
        if (MDAccessNumericWaveData(outPutWaveH, kMDWaveAccessMode0, &outPutOffset)) throw result = WAVEERROR_NOS;
        outPutDataStartPtr = (char*)(*outPutWaveH) + outPutOffset;
        // threads
        layerSize = inputDimSizes [COLUMNS] * inputDimSizes [ROWS];
        if (layerSize < gNumProcessors){
            nThreads = layerSize;
        }else{
            nThreads = gNumProcessors;
        }
        threadPix = layerSize/nThreads;     // number of pixels done by a thread, note integer truncation
        // make an array of parameter structures
        paramArrayPtr = (KalmanNextThreadParamsPtr)WMNewPtr(nThreads * sizeof(KalmanNextThreadParams));
        if (paramArrayPtr == nullptr) throw result = MEMFAIL;
        // **** make array of threads
#ifdef __GNUC__
        threadsPtr =(pthread_t*)WMNewPtr(nThreads * sizeof(pthread_t));
        if (threadsPtr == nullptr) throw result = MEMFAIL;
#endif
#ifdef _WINDOWS_
        threadsPtr = (HANDLE*)WMNewPtr(nThreads * sizeof(HANDLE));
        threadIDsPtr = (DWORD*)WMNewPtr(nThreads * sizeof(DWORD));
        if ((threadsPtr == nullptr) || (threadIDsPtr == nullptr)) throw result = MEMFAIL;
#endif
       
    }catch (int result){
        if (paramArrayPtr != nullptr)  WMDisposePtr ((Ptr)paramArrayPtr);
        if (threadsPtr != nullptr) WMDisposePtr ((Ptr)threadsPtr);
#ifdef _WINDOWS_
        if (threadIDsPtr != nullptr) WMDisposePtr ((Ptr)threadIDsPtr);
#endif
        p -> result = (double)(result - FIRST_XOP_ERR);
#ifdef NO_IGOR_ERR
    return (0);
#else
    return (result);
#endif
    }
    // fill pramaters array
    for (iThread = 0; iThread < nThreads; iThread++){
        paramArrayPtr[iThread].inPutWaveType = inPutWaveType;
        paramArrayPtr[iThread].inPutDataStartPtr = inPutDataStartPtr;
        paramArrayPtr[iThread].outPutDataStartPtr = outPutDataStartPtr;
        paramArrayPtr[iThread].iKal = iKal;
        paramArrayPtr[iThread].threadOffset = threadPix * iThread;
        paramArrayPtr[iThread].threadPix = threadPix;
    }
    // last thread gets any left-over pixels
    paramArrayPtr[nThreads -1].threadPix += layerSize % nThreads;
#ifdef __GNUC__
    for (iThread = 0; iThread < nThreads; iThread++){
        pthread_create (&threadsPtr[iThread], NULL, KalmanNextThread, (void *) &paramArrayPtr[iThread]);
    }
    // Wait till all the threads are finished
    for (iThread = 0; iThread < nThreads; iThread++){
        pthread_join (threadsPtr[iThread], NULL);
    }
#endif
#ifdef _WINDOWS_
    for (iThread = 0; iThread < nThreads; iThread++){
        // create the threads with Default security attributes, Default Stack size, KalmanThread function, pointer to a parameter array, Default creation flags, and a pointer to a DWORD to hold threadId
        threadsPtr[iThread] = CreateThread( NULL, 0, (LPTHREAD_START_ROUTINE)KalmanNextThread, (LPVOID)&paramArrayPtr[iThread], 0, &threadIDsPtr[iThread]);
    }
    // wait for the threads to finish
    WaitForMultipleObjects(nThreads, threadsPtr, TRUE, INFINITE);
    // Close thread handles
    for (iThread = 0; iThread < nThreads; iThread++) {
        CloseHandle(threadsPtr[iThread]);
    }
    // free thread structures
    WMDisposePtr((Ptr)threadIDsPtr);    // freee threadIDs, only used on Windows
#endif
    WMDisposePtr((Ptr)threadsPtr);     // free memory for thread pointers Array
    WMDisposePtr((Ptr)paramArrayPtr);  // Free paramaterArray memory
    // Inform Igor that we have changed the input wave.
    WaveHandleModified(outPutWaveH);
    p -> result = (0);
    return (0);
}


/* ********************************************** KalmanListT **********************************************
 Template for handling all data types for KalmanList function
Last Modified 2013/07/16 by Jamie Boyd  */
template <typename T> int KalmanListT (T** srcWaveStarts, T* destWaveStart, UInt16 nWaves, CountInt startPos, CountInt endPos, float multiplier) {
    UInt16 iWave;
    CountInt iPos;
    T** srcWave;
    T** srcWaveEnd = srcWaveStarts + nWaves;
    /* If multiplier is < 1, use standard averaging across waves with a floating point temporary value */
    if (multiplier < 1){
        double tempVal;
        for (iPos = startPos; iPos < endPos; iPos++){
            for (tempVal = 0, srcWave = srcWaveStarts; srcWave < srcWaveEnd; srcWave++){
                tempVal += *(*srcWave + iPos);
            }
            *(destWaveStart + iPos) = (tempVal/nWaves);
        }
    }else{ //Kalman
        // Do Special stuff for first wave
        if (*srcWaveStarts == destWaveStart){    //this happens when collapsing a wave into the first frame
            if (multiplier > 1){
                //Multiply the first wave by Multiplier
                for (iPos = startPos; iPos < endPos; iPos++){
                    *(destWaveStart + iPos) *= multiplier;
                }
            }
        }else{ //srcwave and destwave are different
            if (multiplier > 1){
                // Set output wave = first input wave * Multiplier
                for (iPos=startPos ; iPos < endPos; iPos++){
                    *(destWaveStart + iPos) = *(*srcWaveStarts + iPos) * multiplier;
                }
            }else{ // first wave when No Multiplier, and Src and dest are different
                for (iPos = startPos ; iPos < endPos ; iPos++){
                    *(destWaveStart + iPos) = *(*srcWaveStarts + iPos);
                }
            }
        }
        //For each remaining wave in input list, iterate through, averaging the input value into the output wave
        if (multiplier > 1){
            for(iWave=1, srcWave = srcWaveStarts + 1; iWave < nWaves ; iWave++, srcWave++) {
                for (iPos = startPos ; iPos < endPos; iPos++){
                    *(destWaveStart + iPos) = ((*(destWaveStart + iPos) * iWave) + *(*srcWave + iPos) * multiplier)/(iWave + 1);
                }
            }
            for (iPos =startPos; iPos < endPos; iPos ++){
                *(destWaveStart + iPos) /= multiplier;
            }
        }else{ // no multiplier
            for(iWave=1, srcWave = srcWaveStarts + 1; iWave < nWaves; iWave++, srcWave++) {
                for (iPos = startPos; iPos < endPos ; iPos++){
                    *(destWaveStart + iPos) = ((*(destWaveStart + iPos) * iWave) +  *(*srcWave + iPos))/(iWave + 1);
                }
            }
        }
    }
    return 0;
}


/* Structure to pass data to each KalmanListThread
Last Modified 2013/07/16 by Jamie Boyd */
typedef struct KalmanListThreadParams{
    int inPutWaveType;
    Ptr* inPutDataStartsPtr;
    char* outPutDataStartPtr;
    UInt16 nWaves;
    float multiplier;
    CountInt threadOffset;
    CountInt threadPix;
   // UInt8 ti; // number of this thread, starting from 0
    //UInt8 tN; // total number of threads
} KalmanListThreadParams, *KalmanListThreadParamsPtr;


/* Each thread to average a list of waves starts with this function
Last Modified 2013/07/16 by Jamie Boyd */
#ifdef __GNUC__
    void* KalmanListThread (void* threadarg){
#endif
#ifdef _WINDOWS_
DWORD WINAPI KalmanListThread(LPVOID threadarg) {
#endif
    struct KalmanListThreadParams* p;
    p = (struct KalmanListThreadParams*) threadarg;
   // CountInt nPnts= p->nPnts;
   // float multiplier = p->multiplier;
    //UInt8 ti = p->ti;
    //UInt8 tN = p->tN;
    //CountInt pntsPerThread = nPnts/tN;
    //CountInt startPos = ti * pntsPerThread; // which point to start this thread on depends on thread number * points per thread. ti is 0 based
   // if (ti == (tN - 1)) pntsPerThread += (nPnts % tN); // last thread gets any extra points
    switch (p->inPutWaveType) {
    case NT_I8:
        KalmanListT ((char**)p->inPutDataStartsPtr, (char*)p->outPutDataStartPtr, p->nWaves, p->threadOffset, p->threadOffset + p->threadPix, p->multiplier);
        break;
    case (NT_I8 | NT_UNSIGNED):
        KalmanListT ((unsigned char**)p->inPutDataStartsPtr, (unsigned char*)p->outPutDataStartPtr, p->nWaves, p->threadOffset, p->threadOffset + p->threadPix, p->multiplier);
        break;
    case NT_I16:
        KalmanListT ((short**)p->inPutDataStartsPtr, (short*)p->outPutDataStartPtr, p->nWaves, p->threadOffset, p->threadOffset + p->threadPix, p->multiplier);
    case (NT_I16 | NT_UNSIGNED):
        KalmanListT ((unsigned short**)p->inPutDataStartsPtr, (unsigned short*)p->outPutDataStartPtr, p->nWaves, p->threadOffset, p->threadOffset + p->threadPix, p->multiplier);
        break;
    case NT_I32:
        KalmanListT ((SInt32**)p->inPutDataStartsPtr, (SInt32*)p->outPutDataStartPtr, p->nWaves, p->threadOffset, p->threadOffset + p->threadPix, p->multiplier);
        break;
    case (NT_I32| NT_UNSIGNED):
        KalmanListT ((UInt32**)p->inPutDataStartsPtr, (UInt32*)p->outPutDataStartPtr, p->nWaves, p->threadOffset, p->threadOffset + p->threadPix, p->multiplier);
        break;
    case NT_I64:
        KalmanListT((SInt64**)p->inPutDataStartsPtr, (SInt64*)p->outPutDataStartPtr, p->nWaves, p->threadOffset, p->threadOffset + p->threadPix, p->multiplier);
        break;
    case (NT_I64 | NT_UNSIGNED):
        KalmanListT((UInt64**)p->inPutDataStartsPtr, (UInt64*)p->outPutDataStartPtr, p->nWaves, p->threadOffset, p->threadOffset + p->threadPix, p->multiplier);
        break;
    case NT_FP32:
        KalmanListT ((float**)p->inPutDataStartsPtr, (float*)p->outPutDataStartPtr, p->nWaves, p->threadOffset, p->threadOffset + p->threadPix, p->multiplier);
        break;
    case NT_FP64:
        KalmanListT ((double**)p->inPutDataStartsPtr, (double*)p->outPutDataStartPtr, p->nWaves, p->threadOffset, p->threadOffset + p->threadPix, p->multiplier);
        break;
    }
#ifdef __GNUC__
    return nullptr;
#endif
#ifdef _WINDOWS_
    return 0;
#endif
}

/* *************************************************************** KalmanList ***************************************************************************************
 KalmanList XOP entry function
 Averages a semicolon-separated list of 2D waves. Each wave must have same data type and same dimensions. This is not used by the
 twoPhoton acquisition code so we print some more information in the error cases with XOPNotice
 KalmanListParams
 inPutList          semicolon separated list of input waves, with paths
 outPutPath         path and wavename of output wave to make
 multiplier         Multiplier for 16 bit waves containing less than 16 bits of data
 overwrite          0 to give errors when output wave already exists. non-zero to overwrite existing wave.
 result             0 for success, else error code
 Last Modified 2014/02/13 by Jamie Boyd */
extern "C" int KalmanList (KalmanListParamsPtr p) {
    int result = 0;                             // The error returned from various Wavemetrics functions
    waveHndl outPutWaveH = nullptr;             // handle to output wave
    waveHndl* handleList = nullptr;             // pointer to an array of handles for input waves
    DFPATH inPutPath, outPutPath ;              // string to hold data folder path of input wave
    WVNAME inPutWaveName, outPutWaveName;       // C strings to hold names of input and output waves
    DataFolderHandle inPutDFHandle;             // Handle to datafolders of input and output waves
    DataFolderHandle outPutDFHandle;            // Handles to datafolders of input and output waves
    int inPutWaveType;                          //  Wavetypes numeric codes for things like 32 bit floating point, 16 bit int, etc
    int inPutDimensions;                        // number of numDimensions in input and output waves
    CountInt inputDimSizes[MAX_DIMENSIONS+1];   // an array used to hold the width, height, layers, and chunk sizes
    char* outPutDataStartPtr;                   // Pointer to start of data in output wave
    Ptr* inPutDataStartsPtr = nullptr;          // Pointer to an array of pointers for starts of input data
    UInt8 overWriteOK = p->overwrite;           // if it is O.K. to overwrite an existing output wave
    UInt8 isOverWriting = 0;                    // 0 if using a separate output wave, 1 for overwriting output wave
    float multiplier = p->multiplier;
    UInt16 numWaves;                            //number of input waves in the input list
    CountInt waveOffset;                        //offset in bytes from begnning of handle to a wave to the actual data - size of headers, units, etc.
    CountInt threadPix;
    CountInt layerSize;
    UInt8 iThread, nThreads;
    KalmanListThreadParamsPtr paramArrayPtr = nullptr;
    //**** declare pthreads or HANDLE pointer
    #ifdef __GNUC__
    pthread_t* threadsPtr = nullptr;
    #endif
    #ifdef  _WINDOWS_
    HANDLE* threadsPtr = nullptr;
    DWORD *threadIDsPtr = nullptr;     // pointer to an array of DWORD (32-bit unsigned integer) that Windows uses for thread IDs
    #endif
    try {
        // Check that input string exists
        if (WMGetHandleSize (p->inPutList) == 0) throw result = NON_EXISTENT_WAVE;
        // If outPutPath is empty string, we are overwriting first wave in list with results
        if (WMGetHandleSize (p->outPutPath) == 0){
            if (!(overWriteOK)) throw result = OVERWRITEALERT;
            isOverWriting = 1;
        }else{ // Parse outPut path
            ParseWavePath (p->outPutPath, outPutPath, outPutWaveName);
            //check that data folder is valid and get a handle to the datafolder
            if (GetNamedDataFolder (NULL, outPutPath, &outPutDFHandle))throw result = WAVEERROR_NOS;
        }
        // parse input list into an array of waveHandles
        handleList = ParseWaveListPaths (p->inPutList, &numWaves);
        // check that dimension sizes and wave types (no text waves) are the same for each wave in array
#ifdef IPARSEWAVE_INF0
        char XOPbuffer [256]; // string used for XOPAlert if we find a bad wave
#endif
        WVNAME tInPutWaveName; // temp wave name  for each wave in array
        DFPATH tInPutPath; // temp datafolder path for each wave in array
        DataFolderHandle tInPutDFHandle;
        int tInPutWaveType; // temp value for wave type of each wave in the array
        int tInPutDimensions;    // temp number of numDimensions for each wave in array
        CountInt tinputDimSizes[MAX_DIMENSIONS+1];    // temp width, height, layers, and chunk sizes for each wave in array
        // get info for first wave in list
        if (handleList [0] == NULL){
#ifdef IPARSEWAVE_INF0
            sprintf(XOPbuffer, "The specification for wave %d in the input list was bad.\r", 0);
            XOPNotice (XOPbuffer);
#endif
            throw result = BADWAVEINLIST;
        }
        inPutWaveType = WaveType(handleList[0]);
        if (inPutWaveType==TEXT_WAVE_TYPE){
#ifdef IPARSEWAVE_INF0
            sprintf(XOPbuffer, "Wave %d in the input list was a text wave.\r", 0);
            XOPNotice (XOPbuffer);
#endif
            throw result = NOTEXTWAVES;
        }
        // Get wave dimensions and calculate number of points
        if (MDGetWaveDimensions(handleList[0], &inPutDimensions, inputDimSizes))throw result = WAVEERROR_NOS;
        if (inPutDimensions != 2) throw result = INPUTNEEDS_2D_WAVE;
        // check to see if output wave is the same as the 1st input wave
        WaveName (handleList[0], inPutWaveName);
        GetWavesDataFolder (handleList[0], &inPutDFHandle);
        GetDataFolderNameOrPath (inPutDFHandle, 1, inPutPath);
        if (isOverWriting == 0){
            if ((!(CmpStr (inPutPath,outPutPath))) && (!(CmpStr (inPutWaveName,outPutWaveName)))) throw result = OVERWRITEALERT;
        }
        // check values for other waves in array against values for first wave
        for (int iw = 1; iw < numWaves; iw++){
            int id;
            // check that handle is good
            if (handleList [iw] == nullptr){
#ifdef IPARSEWAVE_INF0
                sprintf(XOPbuffer, "The specification for wave %d in the input list was bad.\r", iw);
                XOPNotice (XOPbuffer);
#endif
                throw result = BADWAVEINLIST;
            }
            // check input type
            tInPutWaveType = WaveType(handleList[iw]);
            if (tInPutWaveType==TEXT_WAVE_TYPE){
#ifdef IPARSEWAVE_INF0
                sprintf(XOPbuffer, "Wave %d in the input list was a text wave.\r", iw);
                XOPNotice (XOPbuffer);
#endif
                throw result = NOTEXTWAVES;
            }
            if (tInPutWaveType != inPutWaveType) throw result = NOTSAMEWAVETYPE;
            // check number of dimensions
            if (MDGetWaveDimensions(handleList[iw], &tInPutDimensions, tinputDimSizes))throw result = WAVEERROR_NOS;
            if (tInPutDimensions != inPutDimensions){
#ifdef IPARSEWAVE_INF0
                sprintf(XOPbuffer, "The number of dimensions of wave %d in the input list did not match the number of dimenisons of the first wave in the list.\r", iw);
                XOPNotice (XOPbuffer);
#endif
                throw result = NOTSAMEDIMSIZE;
            }
            // check sizes of each dimension
            for (id=0; id < MAX_DIMENSIONS; id +=1){
                if (tinputDimSizes [id] != inputDimSizes [id]){
#ifdef IPARSEWAVE_INF0
                    sprintf(XOPbuffer, "The %d dimension size of wave %d in the input list did not match the corresponding dimensions size of the first wave in the list.\r", id, iw);
                    XOPNotice (XOPbuffer);
#endif
                    throw result = NOTSAMEDIMSIZE;
                }
            }
            // Check wavename for overwriting output wave
            WaveName (handleList[iw], tInPutWaveName);
            GetWavesDataFolder (handleList[iw], &tInPutDFHandle);
            GetDataFolderNameOrPath (tInPutDFHandle, 1, tInPutPath);
            if ((!(CmpStr (tInPutPath, outPutPath))) && (!(CmpStr (tInPutWaveName, outPutWaveName)))){
#ifdef IPARSEWAVE_INF0
                sprintf(XOPbuffer, "The output wave specified would overwrite the %d wave in the input list.\r", iw);
                XOPNotice (XOPbuffer);
#endif
                throw result = OVERWRITEALERT;
            }
        }
        // make the output wave, unless overwriting first wave in list
        //No liberal wave names for output wave
        if (isOverWriting == 0){
            CleanupName (0, outPutWaveName, MAX_OBJ_NAME);
            if ( MDMakeWave (&outPutWaveH, outPutWaveName, outPutDFHandle, inputDimSizes, inPutWaveType, overWriteOK)) throw result = WAVEERROR_NOS;
        }
        // get offsets to data for input waves
        inPutDataStartsPtr = (Ptr*) WMNewPtr (numWaves * sizeof (Ptr));
        for (int iw = 0; iw < numWaves; iw++){
            if (MDAccessNumericWaveData(handleList[iw], kMDWaveAccessMode0, &waveOffset)) throw result = WAVEERROR_NOS;
            *(inPutDataStartsPtr + iw) = (char*)(*handleList[iw]) + waveOffset;
        }
        // get offset for outPut wave
        if (isOverWriting) {
            outPutDataStartPtr = *inPutDataStartsPtr;
        }else{
            if (MDAccessNumericWaveData(outPutWaveH, kMDWaveAccessMode0, &waveOffset)) throw result = WAVEERROR_NOS;
            outPutDataStartPtr =  (char*)(*outPutWaveH) + waveOffset;
        }
        // multiprocessor init
        // threads
        layerSize = inputDimSizes[COLUMNS] * inputDimSizes[ROWS];
        if (layerSize < gNumProcessors) {
            nThreads = layerSize;
        }
        else {
            nThreads = gNumProcessors;
        }
        threadPix = layerSize / nThreads;     // number of pixels done by a thread, note integer truncation
        paramArrayPtr = (KalmanListThreadParamsPtr)WMNewPtr(nThreads * sizeof(KalmanListThreadParams));
        if (paramArrayPtr == nullptr) throw result = MEMFAIL;
       
        // **** make array of threads
#ifdef __GNUC__
        threadsPtr =(pthread_t*)WMNewPtr(nThreads * sizeof(pthread_t));
        if (threadsPtr == nullptr) throw result = MEMFAIL;
#endif
#ifdef _WINDOWS_
        threadsPtr = (HANDLE*)WMNewPtr(nThreads * sizeof(HANDLE));
        threadIDsPtr = (DWORD*)WMNewPtr(nThreads * sizeof(DWORD));
        if ((threadsPtr == nullptr) || (threadIDsPtr == nullptr)) throw result = MEMFAIL;
#endif
       
    }catch (int result){
        if (threadsPtr != nullptr) WMDisposePtr((Ptr)threadsPtr);
        if (paramArrayPtr != nullptr) WMDisposePtr ((Ptr)paramArrayPtr);
#ifdef _WINDOWS_
        if (threadIDsPtr != nullptr) WMDisposePtr ((Ptr)threadIDsPtr);
#endif
        if (inPutDataStartsPtr != nullptr) WMDisposePtr ((Ptr)inPutDataStartsPtr);
        WMDisposeHandle(p->inPutList);
        WMDisposeHandle(p->outPutPath);
        // set result
        p -> result = (result - FIRST_XOP_ERR);    // XFUNC error code
        return (result);
    }
    // fill threadArray
    for (iThread = 0; iThread < nThreads; iThread++){
        paramArrayPtr[iThread].inPutWaveType = inPutWaveType;
        paramArrayPtr[iThread].inPutDataStartsPtr = inPutDataStartsPtr;
        paramArrayPtr[iThread].outPutDataStartPtr = outPutDataStartPtr;
        paramArrayPtr[iThread].nWaves=numWaves;
        paramArrayPtr[iThread].multiplier = multiplier;
        paramArrayPtr[iThread].threadOffset = iThread * threadPix;
        paramArrayPtr[iThread].threadPix = threadPix;
     }
    paramArrayPtr[nThreads - 1].threadPix += layerSize % nThreads;
#ifdef __GNUC__
    for (iThread = 0; iThread < nThreads; iThread++){
        pthread_create (&threadsPtr[iThread], NULL, KalmanListThread, (void *) &paramArrayPtr[iThread]);
    }
    // Wait till all the threads are finished
    for (iThread = 0; iThread < nThreads; iThread++){
        pthread_join (threadsPtr[iThread], NULL);
    }
#endif
#ifdef _WINDOWS_

    for (iThread = 0; iThread < nThreads; iThread++) {
        // create the threads with Default security attributes, Default Stack size, KalmanThread function, pointer to a parameter array, Default creation flags, and a pointer to a DWORD to hold threadId
        threadsPtr[iThread] = CreateThread(NULL, 0, (LPTHREAD_START_ROUTINE)KalmanListThread, (LPVOID) & paramArrayPtr[iThread], 0, &threadIDsPtr[iThread]);
    }
    // wait for the threads to finish
    WaitForMultipleObjects(nThreads, threadsPtr, TRUE, INFINITE);
    // Close thread handles
    for (iThread = 0; iThread < nThreads; iThread++) {
        CloseHandle(threadsPtr[iThread]);
    }
    WMDisposePtr((Ptr) threadIDsPtr);
#endif
    WMDisposePtr ((Ptr)threadsPtr);         // free memory for threads Array
    WMDisposePtr ((Ptr)paramArrayPtr);      // Free paramaterArray memory
  
    WMDisposePtr ((Ptr)inPutDataStartsPtr); // free pointers to data starts
    WMDisposeHandle(p->inPutList);          // free inPutList input string
    WMDisposeHandle(p->outPutPath);         // free outPutPath input string
    // Inform Igor that we have changed the wave.
    WaveHandleModified(outPutWaveH);
    // set result
    p -> result = (0);
    return (0);
}

