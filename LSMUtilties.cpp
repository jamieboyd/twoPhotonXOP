#include "twoPhoton.h"
#include "math.h"

/* ------------------------------LSM Utilities --------------------------------------------------
utility functions specialized for Laser Scanning Microscope data acquisition
Last Modified 2025/06/23 by Jamie Boyd
 -------------------------------------------------------------------------------------------------*/
 
/* --------------------------- GetSetNumProcessors--------------------------------------------
sanity check function for getting number of processors should return same as ThreadProcessorCount
Last Modified 2025/06/24 by Jamie Boyd */
extern "C" int GetSetNumProcessors(GetSetNumProcessorsParamsPtr p){
    gNumProcessors = num_processors();
    p->result = (double)gNumProcessors;
    return (0);
}

/* -------------------------------- SwapEven -------------------------------------------------------
 horizontally swaps every other line in an image or series of images.
 Used after doing back and forth scanning, where every other line is scanned from the opposite direction.
 ----------------------------------------------------------------------------------------------------- */

/* Template for SwapEven function
Last Modified 2013/07/15 by Jamie Boyd */
template <typename T> void SwapEvenT (T *dataStartPtr, CountInt numLines, CountInt lineLen) {
    // make Pointer to end of the data
    T* dataEndPtr = dataStartPtr + (numLines * lineLen);
    // make pointers for start and end of each line
    T* lineStartPtr;
    T* lineEndPtr;
    //temporary value for swapping
    T temp;
    // calculate distance to start and end of next line
    UInt16 toNextStart = (2 * lineLen) - (lineLen/2);    //due to rounding of integer math, it doesn't matter if LineLen is even or odd
    UInt16 toNextEnd = (2 * lineLen) + (lineLen/2);
    //iterate through every other line, flipping it around
    for(lineStartPtr = dataStartPtr + lineLen, lineEndPtr = lineStartPtr + lineLen -1; lineEndPtr < dataEndPtr; lineStartPtr += toNextStart,lineEndPtr += toNextEnd) {
        while (lineStartPtr < lineEndPtr){
            SWAP (*lineStartPtr, *lineEndPtr);
            lineStartPtr++;
            lineEndPtr--;
        }
    }
}

/* Structure to pass data to each SwapEvenThread
Last Modified 2013/07/15 by Jamie Boyd */
typedef struct SwapEvenThreadParams{
    int inPutWaveType;
    char* dataStartPtr;
    CountInt numLines;
    CountInt lineLen;
    UInt8 ti; // number of this thread, starting from 0
    UInt8 tN; // total number of threads (255 "should be enough for anyone")
}SwapEvenThreadParams, *SwapEvenThreadParamsPtr;


/* Each thread to swap a range of rows starts with this function
Last Modified 2025/06/13 by Jamie Boyd */
void* SwapEvenThread (void* threadarg){
    struct SwapEvenThreadParams* p;
    p = (struct SwapEvenThreadParams*) threadarg;
    CountInt tLines = (p-> numLines/p->tN);
    if (tLines % 2) tLines--; // make sure number of lines per thread is even
    CountInt startPos = p->ti * tLines * p->lineLen;
    // the last thread gets any left-over lines
    if (p->ti == p->tN - 1) tLines +=  (p->numLines % p->tN);
    // call the right template function for the wave types
    switch (p->inPutWaveType) {
    case NT_I8:
        SwapEvenT ((char*) p->dataStartPtr + startPos, tLines, p->lineLen);
        break;
    case (NT_I8 | NT_UNSIGNED):
        SwapEvenT ((unsigned char*)p->dataStartPtr + startPos, tLines, p->lineLen);
        break;
    case NT_I16:
        SwapEvenT ((short*)p->dataStartPtr + startPos, tLines, p->lineLen);
        break;
    case (NT_I16 | NT_UNSIGNED):
        SwapEvenT ((unsigned short*)p->dataStartPtr + startPos, tLines, p->lineLen);
        break;
    case NT_I32:
        SwapEvenT ((long*)p->dataStartPtr + startPos, tLines, p->lineLen);
        break;
    case (NT_I32| NT_UNSIGNED):
        SwapEvenT ((unsigned long*)p->dataStartPtr + startPos, tLines, p->lineLen);
        break;
    case NT_FP32:
        SwapEvenT ((float*)p->dataStartPtr + startPos, tLines, p->lineLen);
        break;
    case NT_FP64:
        SwapEvenT ((double*)p->dataStartPtr + startPos, tLines, p->lineLen);
        break;
        break;
    }
    return nullptr;
}

/* SwapEven XOP entry function
 SwapEvenParams:
 wave handle to start of data of input wave, which is overwrittten
 result =  0 or error code
Last Modified 2025/06/23 by Jamie Boyd */
extern "C" int SwapEven (SwapEvenParamsPtr p){
    waveHndl wavH = nullptr;        // handle to the input wave
    int waveType; //  Wavetypes numeric codes for things like 32 bit floating point, 16 bit int, etc
    int numDimensions;    // number of dimensions in input and output waves
    CountInt dimensionSizes[MAX_DIMENSIONS+1];    // an array used to hold the width, height, layers, and chunk sizes
    IndexInt dataOffset;    //offset in bytes from begnning of handle to a wave to the actual data - size of headers, units, etc.
    CountInt lineLen;    // The length of each line in the image
    CountInt numLines;    // The number of lines in the file
    int result;    // The error returned from various Wavemetrics functions
    char* dataStartPtr;    // Pointer to start of data in input wave. Need to use char for these to use WM function to get data offset
    // Multi threading
    UInt8 iThread, nThreads;
    SwapEvenThreadParamsPtr paramArrayPtr = nullptr;
    pthread_t* threadsPtr = nullptr;
    try {
        // Get handle to input wave. Make sure it exists.
        wavH = p->w1;
        if (wavH == nullptr) throw result = NON_EXISTENT_WAVE;
        // Get wave data type.
        waveType = WaveType(wavH);
        // Can't process text waves
        if (waveType == TEXT_WAVE_TYPE) throw result = NOTEXTWAVES;
        // Get number of used dimensions in wave.
        if (MDGetWaveDimensions(wavH, &numDimensions, dimensionSizes)) throw result = WAVEERROR_NOS;
        // Check that wave1 is 2D or 3D
        if ((numDimensions == 1) || (numDimensions == 4)) throw result = INPUTNEEDS_2D3D_WAVE;
        // Get dimension size info and calculate number of lines to process
        lineLen = dimensionSizes[0];
        if (numDimensions == 2){
            numLines = dimensionSizes[1];
        }else{
            numLines = dimensionSizes[1] * dimensionSizes[2];
        }
        // Get the offsets to the data in the wave
        if (MDAccessNumericWaveData(wavH, kMDWaveAccessMode0, &dataOffset)) throw result = WAVEERROR_NOS;
        dataStartPtr = (char*)(*wavH) + dataOffset;
        // Multi threading
        nThreads = gNumProcessors;
        paramArrayPtr = (SwapEvenThreadParamsPtr)WMNewPtr(nThreads * sizeof(SwapEvenThreadParams));
        if (paramArrayPtr == nullptr) throw result = MEMFAIL;
        // make an array of pthread_t
        threadsPtr =(pthread_t*)WMNewPtr(nThreads * sizeof(pthread_t));
        if (threadsPtr == nullptr) throw result = MEMFAIL;
    }catch (int result){
        if (paramArrayPtr != nullptr){
            WMDisposePtr ((Ptr)paramArrayPtr);
            if (threadsPtr != nullptr){
                WMDisposePtr ((Ptr)threadsPtr);
            }
        }
        p -> result = (double)(result - FIRST_XOP_ERR);
#ifdef NO_IGOR_ERR
        return (0);
#else
        return (result);
#endif
    }
    // fill threadParams array
    for (iThread = 0; iThread < nThreads; iThread++){
        paramArrayPtr[iThread].inPutWaveType = waveType;
        paramArrayPtr[iThread].dataStartPtr = dataStartPtr;
        paramArrayPtr[iThread].numLines = numLines;
        paramArrayPtr[iThread].lineLen = lineLen;
        paramArrayPtr[iThread].ti=iThread; // number of this thread, starting from 0
        paramArrayPtr[iThread].tN =nThreads; // total number of threads
    }
    // create the threads
    for (iThread = 0; iThread < nThreads; iThread++){
        pthread_create (&threadsPtr[iThread], NULL, SwapEvenThread, (void *) &paramArrayPtr[iThread]);
    }
    // Wait till all the threads are finished
    for (iThread = 0; iThread < nThreads; iThread++){
        pthread_join (threadsPtr[iThread], NULL);
    }
    // free memory for pThreads Array
    WMDisposePtr ((Ptr)threadsPtr);
    // Free paramaterArray memory
    WMDisposePtr ((Ptr)paramArrayPtr);
    // Inform Igor that we have changed the input wave.
    WaveHandleModified(wavH);
    p -> result = (0);
    return (0);
}

/* -------------------------------------DownSample--------------------------------------------------------
Down Sample takes a wave and resizes its X dimension, combining boxfactor points by taking the
 average, sum, max, or median value of boxFactor points. Used after "oversampling" during scanning,
 i.e., taking several A/D conversions per pixel
DownSample methods: 1 = average, 2 = sum, 3 = max, 4 = median
------------------------------------------------------------------------------------------------------------ */


/* Template functions for the various DownSample methods
Case 1 Kalman-style averaging.
Last Modified 2013/07/15 by Jamie Boyd */
template <typename T> void DSaverageT (T *dataStartPtr, CountInt points, UInt16 boxFactor) {
    // make Pointer to the end of the data
    T *dataEndPtr = dataStartPtr + points;
    T *outPutPtr; // output position varies slowly
    T *inPutPtr;  // input position varies quickly
    UInt16 ifactor;
    //iterate through the wave, averaging input range into output value
    for (outPutPtr = inPutPtr = dataStartPtr;inPutPtr < dataEndPtr; outPutPtr++){
        *outPutPtr = *inPutPtr;
        inPutPtr++;
        for (ifactor = 1; ifactor < boxFactor; ifactor++, inPutPtr++){
            *outPutPtr = (*outPutPtr * (ifactor/(ifactor + 1))) + *inPutPtr /(ifactor + 1);
        }
    }
}

// Case 2: Summing. Need different function for different bounds checking for floats, signed ints, and unsigned ints
//Floats: Don't check Max or Min - overflows should automatically go to INF and -INF
template <typename T> void DSsumFloatT (T *dataStartPtr, CountInt points, UInt16 boxFactor) {
    T *dataEndPtr = dataStartPtr + points;
    T *outPutPtr;
    T *inPutPtr;
    UInt16 ifactor;
    for (outPutPtr = inPutPtr = dataStartPtr;inPutPtr < dataEndPtr; outPutPtr++){
        for (*outPutPtr = *inPutPtr,inPutPtr++, ifactor = 1; ifactor < boxFactor; ifactor++, inPutPtr++)
            *outPutPtr = *outPutPtr + *inPutPtr;
    }
}

// Signed Ints: Check OverFlows at both ends (Max and Min)
template <typename T> void DSsumSignedIntT (T *dataStartPtr, CountInt points, UInt16 boxFactor, T Tmax, T Tmin){
    T *dataEndPtr = dataStartPtr + points;
    T *outPutPtr;
    T *inPutPtr;
    UInt16 ifactor;
    for (outPutPtr = inPutPtr = dataStartPtr;inPutPtr < dataEndPtr; outPutPtr++){
        for (*outPutPtr = *inPutPtr,inPutPtr++, ifactor = 1; ifactor < boxFactor; ifactor++, inPutPtr++){
            if ((*outPutPtr + *inPutPtr) > Tmax) //Check Max
                *outPutPtr = Tmax;
            else{
                if ((*outPutPtr + *inPutPtr) < Tmin) //Check Min
                    *outPutPtr = Tmin;
                else
                    *outPutPtr = *outPutPtr + *inPutPtr;
            }
        }
    }
}

// Unsigned integer types: Check OverFlows only for Max
template <typename T> void DSsumUnSignedIntT (T *dataStartPtr, CountInt points, UInt16 boxFactor, T Tmax){
    T *dataEndPtr = dataStartPtr + points;
    T *outPutPtr;
    T *inPutPtr;
    UInt16 ifactor;
    for (outPutPtr = inPutPtr = dataStartPtr;inPutPtr < dataEndPtr; outPutPtr++){
        for (*outPutPtr = *inPutPtr,inPutPtr++, ifactor = 1; ifactor < boxFactor; ifactor++, inPutPtr++){
            if ((*outPutPtr + *inPutPtr) > Tmax){  //check Max
                *outPutPtr = Tmax;
                break;    //can move on to the next input point - value can't get any higher
            }
            else
                *outPutPtr = *outPutPtr + *inPutPtr;
        }
    }
}

//Case 3: Take maximum Value
template <typename T> void DSMaxT (T *dataStartPtr, CountInt points, UInt16 boxFactor){
    // make Pointer to the end of the data
    T *dataEndPtr = dataStartPtr + points;
    T *outPutPtr;
    T *inPutPtr;
    UInt16 ifactor;

    for (outPutPtr = inPutPtr = dataStartPtr;inPutPtr < dataEndPtr; outPutPtr++){
        *outPutPtr = *inPutPtr;
        inPutPtr++;
        for (ifactor = 1; ifactor < boxFactor; ifactor++, inPutPtr++){
            if (*outPutPtr < *inPutPtr)
                *outPutPtr = *inPutPtr;
        }
    }
}

//case 4: take median value.
template <typename T> void DSMedianT (T *dataStartPtr, CountInt points, UInt16 boxFactor){
    // make Pointer to the end of the data
    T *dataEndPtr = dataStartPtr + points;
    T *outPutPtr;
    T *inPutPtr;
    // Median variables
    UInt16 i,right=boxFactor -1,j,left=0,mid;
    UInt16 k = (boxFactor)/2;
    T a,temp;
    for (outPutPtr = inPutPtr = dataStartPtr;inPutPtr < dataEndPtr; outPutPtr++, inPutPtr += boxFactor){
        for (;;) {
            if (right <= left + 1) {
                if (right == left + 1 && *(inPutPtr + right) < *(inPutPtr + left)) {
                    SWAP((*(inPutPtr + left)),(*(inPutPtr + right)));
                }
                *outPutPtr = *(inPutPtr + k);
                break;
            } else {
                mid=(left + right) >> 1;
                SWAP((*(inPutPtr + mid)),(*(inPutPtr + left + 1)));
                if (*(inPutPtr + left) > *(inPutPtr + right)) {
                    SWAP((*(inPutPtr + left)),(*(inPutPtr + right)));
                }
                if (*(inPutPtr + left + 1) > *(inPutPtr + right)) {
                    SWAP((*(inPutPtr + left + 1)), ( *(inPutPtr + right)));
                }
                if (*(inPutPtr + left) > *(inPutPtr + left+1)) {
                    SWAP((*(inPutPtr + left)),(*(inPutPtr + left+1)));
                }
                i=left+1;
                j=right;
                a=*(inPutPtr + left+1);
                for (;;) {
                    do i++; while (*(inPutPtr + i) < a);
                    do j--; while (*(inPutPtr + j) > a);
                    if (j < i) break;
                    SWAP((*(inPutPtr + i)),(*(inPutPtr + j)));
                }
                *(inPutPtr + left + 1) = *(inPutPtr + j);
                *(inPutPtr + j) = a;
                if (j >= k) right = j-1;
                if (j <= k) left=i;
            }
        }
    }
}


/* Structure to pass data to each DownSample thread
Last Modified 2013/07/15 by Jamie Boyd */
typedef struct DownSampleThreadParams{
    int inPutWaveType;
    char* dataStartPtr;
    CountInt points;
    CountInt boxFactor;
    UInt8 DSType;
    UInt8 ti; // number of this thread, starting from 0
    UInt8 tN; // total number of threads (255 "should be enough for anyone")
}DownSampleThreadParams, *DownSampleThreadParamsPtr;


/* Each thread to down sampe a range of points  starts with this function
Last Modified 2025/06/23 by Jamie Boyd */
void* DownSampleThread (void* threadarg){
    struct DownSampleThreadParams* p;
    p = (struct DownSampleThreadParams*) threadarg;
    CountInt tPoints = (p->points/p->tN);  // number of points to do per thread
    CountInt startPos = p->ti * tPoints;   // starting position for this thread
    if (p->ti == p->tN - 1) tPoints +=  (p->points % p->tN); // the last thread gets any left-over points
    CountInt boxFactor = p->boxFactor;
    UInt8 DSType = p->DSType;
    char* dataStartPtr = p->dataStartPtr;
    // call the right template function for the down sample type
    switch (p->inPutWaveType) {
        case NT_FP64:
            switch (DSType){
                case 1:    //average
                    DSaverageT ((double*) dataStartPtr + startPos , tPoints, boxFactor);
                    break;
                case 2: // sum
                    DSsumFloatT ((double*)dataStartPtr + startPos, tPoints, boxFactor);
                    break;
                case 3: //max
                    DSMaxT ((double*)dataStartPtr + startPos, tPoints, boxFactor);
                    break;
                case 4: // median
                    DSMedianT ((double*)dataStartPtr + startPos, tPoints, boxFactor);
                    break;
            }
            break;
        case NT_FP32:
            switch (DSType){
                case 1:    //average
                    DSaverageT ((float*) dataStartPtr + startPos, tPoints, boxFactor);
                    break;
                case 2: // sum
                    DSsumFloatT ((float*)dataStartPtr + startPos, tPoints, boxFactor);
                    break;
                case 3: //max
                    DSMaxT ((float*)dataStartPtr + startPos, tPoints, boxFactor);
                    break;
                case 4: // median
                    DSMedianT ((float*)dataStartPtr + startPos, tPoints, boxFactor);
                    break;
            }
            break;
        case (NT_I32 | NT_UNSIGNED):
            switch (DSType){
                case 1:    //average
                    DSaverageT ((unsigned long*) dataStartPtr + startPos, tPoints, boxFactor);
                    break;
                case 2: // sum
                    DSsumUnSignedIntT ((unsigned long*)dataStartPtr + startPos, tPoints, boxFactor,(unsigned long)ULONG_MAX);
                    break;
                case 3: //max
                    DSMaxT ((unsigned long*)dataStartPtr + startPos, tPoints, boxFactor);
                    break;
                case 4: // median
                    DSMedianT ((unsigned long*)dataStartPtr + startPos, tPoints, boxFactor);
                    break;
            }
            break;
        case NT_I32:
            switch (DSType){
                case 1:    //average
                    DSaverageT ((long*) dataStartPtr + startPos, tPoints, boxFactor);
                    break;
                case 2: // sum
                    DSsumSignedIntT ((long*)dataStartPtr + startPos, tPoints, boxFactor,(long)LONG_MAX, (long)LONG_MIN);
                    break;
                case 3: //max
                    DSMaxT ((long*)dataStartPtr + startPos, tPoints, boxFactor);
                    break;
                case 4: // median
                    DSMedianT ((long*)dataStartPtr + startPos, tPoints, boxFactor);
                    break;
            }
            break;
        case (NT_I16 | NT_UNSIGNED):
            switch (DSType){
                case 1:    //average
                    DSaverageT ((unsigned short*) dataStartPtr + startPos, tPoints, boxFactor);
                    break;
                case 2: // sum
                    DSsumUnSignedIntT ((unsigned short*)dataStartPtr + startPos, tPoints, boxFactor,(unsigned short)USHRT_MAX);
                    break;
                case 3: //max
                    DSMaxT ((unsigned short*)dataStartPtr + startPos, tPoints, boxFactor);
                    break;
                case 4: // median
                    DSMedianT ((unsigned short*)dataStartPtr + startPos, tPoints, boxFactor);
                    break;
            }
            break;
        case NT_I16:
            switch (DSType){
                case 1:    //average
                    DSaverageT ((short*) dataStartPtr + startPos, tPoints, boxFactor);
                    break;
                case 2: // sum
                    DSsumSignedIntT ((short*)dataStartPtr + startPos, tPoints, boxFactor,(short)SHRT_MAX, (short)SHRT_MIN);
                    break;
                case 3: //max
                    DSMaxT ((short*)dataStartPtr + startPos, tPoints, boxFactor);
                    break;
                case 4: // median
                    DSMedianT ((short*)dataStartPtr + startPos, tPoints, boxFactor);
                    break;
            }
            break;
        case (NT_I8 | NT_UNSIGNED):
            switch (DSType){
                case 1:    //average
                    DSaverageT ((unsigned char*) dataStartPtr + startPos, tPoints, boxFactor);
                    break;
                case 2: // sum
                    DSsumUnSignedIntT ((unsigned char*)dataStartPtr + startPos, tPoints, boxFactor,(unsigned char)UCHAR_MAX);
                    break;
                case 3: //max
                    DSMaxT ((unsigned char*)dataStartPtr + startPos, tPoints, boxFactor);
                    break;
                case 4: // median
                    DSMedianT ((unsigned char*)dataStartPtr + startPos, tPoints, boxFactor);
                    break;
            }
            break;
        case NT_I8:
            switch (DSType){
                case 1:    //average
                    DSaverageT ((char*) dataStartPtr + startPos, tPoints, boxFactor);
                    break;
                case 2: // sum
                    DSsumSignedIntT ((char*)dataStartPtr + startPos, tPoints, boxFactor,(char)SCHAR_MAX, (char)SCHAR_MIN);
                    break;
                case 3: //max
                    DSMaxT ((char*)dataStartPtr + startPos, tPoints, boxFactor);
                    break;
                case 4: // median
                    DSMedianT ((char*)dataStartPtr + startPos, tPoints, boxFactor);
                    break;
            }
            break;
    }    // end of switch
    return nullptr;
}

/* DownSample XOP entry function
 DownSampleParams:
 wavehandle to input wave, which is overwritten
 downSample type 1 = average, 2 = sum, 3 = max, 4 = median
 boxfactor = number of pixels boxed together
Last Modified 2025/06/13 by Jamie Boyd */
extern "C" int DownSample (DownSampleParamsPtr p) {
    int result = 0;    // The error returned from various Wavemetrics functions
    int waveType; //  Wavetypes numeric codes for things like 32 bit floating point, 16 bit int, etc
    int numDimensions;    // number of dimensions in input  waves
    CountInt dimensionSizes[MAX_DIMENSIONS+1];    // an array used to hold the width, height, layers, and chunk sizes
    BCInt dataOffset;    //offset in bytes from begnning of handle to a wave to the actual data - size of headers, units, etc.
    UInt16 boxFactor;
    CountInt xSize, ySize, zSize, points;
    waveHndl wavH;        // handle to the input wave
    char* dataStartPtr;    // Pointer to start of data in input wave. Need to use char for these to use WM function to get data offset
    UInt8 DSType;
    DownSampleThreadParamsPtr paramArrayPtr = nullptr;
    pthread_t* threadsPtr = nullptr;
    UInt8 iThread, nThreads;
    try{
        // Get handles to input wave. Make sure it exists.
        wavH = p->w1;
        if (wavH == NIL) throw result = NON_EXISTENT_WAVE;
        // Get wave data type.
        waveType = WaveType(wavH);
        // Can't process text waves
        if (waveType == TEXT_WAVE_TYPE) throw result = NOTEXTWAVES;
        // Get number of used dimensions in wave.
        if (MDGetWaveDimensions(wavH, &numDimensions, dimensionSizes)) throw result = WAVEERROR_NOS;
        // Check that input is 2D or 3D
        if ((numDimensions != 2) && (numDimensions!= 3)) throw result = INPUTNEEDS_2D3D_WAVE;
        // Check that Down sample type variable is o.k. 1 = average, 2 = sum, 3 = max, 4 = median
        DSType = (UInt8) p -> dsType;
        if ((DSType > 4) || (DSType < 1)) throw result = BADDSTYPE;
        // Get dimension size info
        xSize = dimensionSizes[0];
        boxFactor = (UInt16) p ->boxFactor;
        if (xSize % boxFactor) throw result = BADFACTOR;
        ySize = dimensionSizes [1];
        if (numDimensions == 2){
            points = xSize * ySize;
        }else{
            zSize = dimensionSizes [2];
            points = xSize * ySize * zSize;
        }
        // get wave handle
        if (MDAccessNumericWaveData(wavH, kMDWaveAccessMode0, &dataOffset)) throw result = WAVEERROR_NOS;
        // make pointer to start of wave data
        dataStartPtr = (char*)(*wavH) + dataOffset;
        // Get ready for Multi threading
        nThreads = gNumProcessors;
        // make an arrray of paramater pointers
        paramArrayPtr = (DownSampleThreadParamsPtr)WMNewPtr(nThreads * sizeof(DownSampleThreadParams));
        // make an array of pthread_t
       threadsPtr = (pthread_t*)WMNewPtr(nThreads * sizeof(pthread_t));
    }catch (int result){
        // dispose of any memory we may have allocated
        if (paramArrayPtr != nullptr) WMDisposePtr ((Ptr)paramArrayPtr);
        if (threadsPtr != nullptr) WMDisposePtr ((Ptr)threadsPtr);
        p -> result = (double)(result - FIRST_XOP_ERR);
#ifdef NO_IGOR_ERR
        return (0);
#else
        return (result);
#endif
    }
    // fill paramater array
    for (iThread = 0; iThread < nThreads; iThread++){
        paramArrayPtr[iThread].inPutWaveType = waveType;
        paramArrayPtr[iThread].dataStartPtr = dataStartPtr;
        paramArrayPtr[iThread].points = points;
        paramArrayPtr[iThread].boxFactor = boxFactor;
        paramArrayPtr[iThread].DSType = DSType;
        paramArrayPtr[iThread].ti = iThread; // number of this thread, starting from 0
        paramArrayPtr[iThread].tN = nThreads; // total number of threads
    }
    // create the threads
    for (iThread = 0; iThread < nThreads; iThread++){
        pthread_create (&threadsPtr[iThread], NULL, DownSampleThread, (void *) &paramArrayPtr[iThread]);
    }
    // Wait till all the threads are finished
    for (iThread = 0; iThread < nThreads; iThread++){
        pthread_join (threadsPtr[iThread], NULL);
    }
    // free memory for pThreads Array
    WMDisposePtr ((Ptr)threadsPtr);
    // Free paramaterArray memory
    WMDisposePtr ((Ptr)paramArrayPtr);
    dimensionSizes [0] = xSize/boxFactor;
    dimensionSizes [1] = ySize;
    dimensionSizes [2] = zSize;
    dimensionSizes [3] = 0;
    MDChangeWave2 (wavH, -1, dimensionSizes, 1);
    // Inform Igor that we have changed the input wave.
    WaveHandleModified (wavH);
    p -> result = (0);
    return (0);
}


/*  ---------------------------------------Transpose Frames------------------------------------------------------------
Transposes each frame in a 3D wave (equivalent to a horizontal flip and a 90 degree counter-clockwise rotation).
Useful because microscope may have an odd number of  mirrors in the light path
 ------------------------------------------------------------------------------------------------------------------- */

/* template function for transposeFrames that are not square
 Last Modified: 2025/06/23 by Jamie Boyd */
template <typename T> void TransposeFramesT (T *dataStartPtr, T *frameCopyStart, CountInt xSize, CountInt ySize, CountInt zSize) {
    CountInt frameSize = (xSize * ySize);
    CountInt frameBytes = (frameSize * sizeof(T));
    // make Pointer to end of the data
    T* dataEndPtr = dataStartPtr + (frameSize * zSize);
    // make pointer for progressing through the given data
    T *dataWavePtr;
    //pointers for iterating through x and y
    T *frameCopyPtr, *frameXend, *frameYend;
    CountInt toNextX = frameSize-1;
    for (dataWavePtr= dataStartPtr; dataWavePtr < dataEndPtr;){
        //make a copy of the current frame into the frame copy buffer
        frameCopyStart = (T*) memcpy (frameCopyStart, dataWavePtr, frameBytes);
        //iterate through x
        for (frameCopyPtr = frameCopyStart, frameXend = frameCopyStart + xSize;frameCopyPtr < frameXend;frameCopyPtr -= toNextX){
            // iterate through y
            for (frameYend = frameCopyPtr + frameSize; frameCopyPtr < frameYend ;frameCopyPtr += xSize, dataWavePtr+=1){
                *dataWavePtr = *frameCopyPtr;
            }
        }
    }
}

/* Template function for TransposeFrames functions with frames that are square, as they often are.
 Does not need to use a frame buffer, so should be a little faster
 Last Modified 2013/07/15 by Jamie Boyd */
template <typename T> void TransposeSquareFramesT (T *dataStartPtr, CountInt xySize, CountInt zSize){
    CountInt frameSize = xySize * xySize;
    // make Pointer to end of the data
    T* dataEndPtr = dataStartPtr + ((frameSize * zSize) - xySize);
    // make pointers for progressing through the data
    T *xPtr, *yPtr;
    // variables to point to the ends of the current line. x and y both the same, so we only need count X
    T *xEndLinePtr;
    // variable to count lines
    CountInt lines;
    // frames - 1, just because it's easier to read this way
    CountInt frameSizeMinus1 = frameSize - 1;
    // temp variable for SWAP macro
    T temp;
    //frames
    for (xPtr = dataStartPtr + 1, yPtr = dataStartPtr + xySize; xPtr < dataEndPtr; xPtr += 1,yPtr += 1){
        //columns/rows
        for (lines = 1; lines < xySize; lines += 1, xPtr += lines, yPtr -= (frameSizeMinus1 - lines*xySize)){
            // points in individual column/row
            for (xEndLinePtr = xPtr + (xySize - lines); xPtr < xEndLinePtr; xPtr += 1,yPtr += xySize){
                SWAP (*xPtr, *yPtr);
            }
        }
    }
}


/* Structure to pass data to each TransposeFrames thread
 Last Modified 2025/06/23 by Jamie Boyd */
typedef struct TransposeFramesThreadParams{
    int inPutWaveType;          // WaveMetrics code for waveType
    char* dataStartPtr;         // pointer to start of input wave, which is overwritten
    char* bufferPtr;            // pointer to start of a frame-sized buffer for temporary calculations, not used for square frames
    CountInt xSize;            // number of columns in each frame
    CountInt ySize;            // number of rows in each frame
    CountInt zSize;             // number of frames in wave
    UInt8 ti; // number of this thread, starting from 0
    UInt8 tN; // total number of threads (255 "should be enough for anyone")
} TransposeFramesThreadParams, *TransposeFramesThreadParamsPtr;


/* Each thread to transpose a range of frames starts with this function
Last Modified 2025/06/23 by Jamie Boyd */
void* TransposeFramesThread (void* threadarg){
    struct TransposeFramesThreadParams* p = (struct TransposeFramesThreadParams*) threadarg;
    CountInt frameSize = p->xSize * p->ySize;
    CountInt tFrames = (p->zSize/p->tN);  // number of frames to do per thread
    CountInt startOffset = p->ti * frameSize * tFrames ;   // starting position for this thread
    if (p->ti == p->tN - 1) tFrames += (tFrames % p->tN); // the last thread gets any left-over frames
    if (p->xSize == p->ySize){
        switch (p->inPutWaveType) {
            case NT_FP64:
                TransposeSquareFramesT ((double*) p->dataStartPtr + startOffset, p->xSize, tFrames);
                break;
            case NT_FP32:
                TransposeSquareFramesT (((float*) p->dataStartPtr) +startOffset, p->xSize, tFrames);
                break;
            case (NT_I32 | NT_UNSIGNED):
                TransposeSquareFramesT ((unsigned long*) p->dataStartPtr + startOffset, p->xSize, tFrames);
                break;
            case NT_I32:
                TransposeSquareFramesT ((long*) p->dataStartPtr + startOffset, p->xSize, tFrames);
                break;
            case (NT_I16 | NT_UNSIGNED):
                TransposeSquareFramesT ((unsigned short*) p->dataStartPtr + startOffset, p->xSize, tFrames);
                break;
            case NT_I16:
                TransposeSquareFramesT ((short*) p->dataStartPtr + startOffset, p->xSize, tFrames);
                break;
            case (NT_I8 | NT_UNSIGNED):
                TransposeSquareFramesT ((unsigned char*) p->dataStartPtr + startOffset, p->xSize, tFrames);
                break;
            case NT_I8:
                TransposeSquareFramesT ((char*) p->dataStartPtr + startOffset, p->xSize, tFrames);
                break;
        }
    }else{
        CountInt bufferOffset =p->ti * frameSize;
        switch (p->inPutWaveType) {
            case NT_FP64:
                TransposeFramesT ((double*) p->dataStartPtr + startOffset, (double*) p->bufferPtr + bufferOffset, p->xSize, p->ySize, tFrames);
                break;
            case NT_FP32:
                TransposeFramesT ((float*) p->dataStartPtr + startOffset, (float*) p->bufferPtr + bufferOffset, p->xSize, p->ySize, tFrames);
                break;
            case (NT_I32 | NT_UNSIGNED):
                TransposeFramesT ((unsigned long*) p->dataStartPtr + startOffset, (unsigned long*) p->bufferPtr + bufferOffset, p->xSize, p->ySize, tFrames);
                break;
            case NT_I32:
                TransposeFramesT ((long*) p->dataStartPtr + startOffset, (long*)p->bufferPtr + bufferOffset, p->xSize, p->ySize, tFrames);
                break;
            case (NT_I16 | NT_UNSIGNED):
                TransposeFramesT ((unsigned short*) p->dataStartPtr + startOffset, (unsigned short*)p->bufferPtr + bufferOffset, p->xSize, p->ySize, tFrames);
                break;
            case NT_I16:
                TransposeFramesT ((short*) p->dataStartPtr + startOffset, (short*)p->bufferPtr + bufferOffset, p->xSize, p->ySize, tFrames);
                break;
            case (NT_I8 | NT_UNSIGNED):
                TransposeFramesT ((unsigned char*) p->dataStartPtr + startOffset, (unsigned char*) p->bufferPtr + bufferOffset, p->xSize, p->ySize, tFrames);
                break;
            case NT_I8:
                TransposeFramesT ((char*) p->dataStartPtr + startOffset, (char*)p->bufferPtr + bufferOffset, p->xSize, p->ySize, tFrames);
                break;
        }
    }
    return nullptr;
 }

/* TransposeFrames XOP entry function
 TransposeFramesParams:
 wave handle to input wave, which is always overwritten
 result which is 0 or error code
 Last Modified 2025/06/23 by Jamie Boyd */
extern "C" int TransposeFrames (TransposeFramesParamsPtr p) {
    int result = 0;    // The error returned from various Wavemetrics functions
    int waveType; //  Wavetypes numeric codes for things like 32 bit floating point, 16 bit int, etc
    int numDimensions;    // number of dimensions in input and output waves
    CountInt dimensionSizes[MAX_DIMENSIONS+1];    // an array used to hold the width, height, layers, and chunk sizes
    CountInt dataOffset;    //offset in bytes from begnning of handle to a wave to the actual data - size of headers, units, etc.
    CountInt xSize;    // The length of each line in each frame
    CountInt ySize; // number of lines in each frame
    CountInt zSize; // number of frames in the stack
    UInt8 is3D; // need to know 2D from 3D when redimensioning
    waveHndl wavH;        // handle to the input wave
    char* dataStartPtr;    // Pointer to start of data in input wave. Need to use char for these to use WM function to get data offset
    // threading
    UInt8 iThread, nThreads;
    TransposeFramesThreadParamsPtr paramArrayPtr = nullptr;
    pthread_t* threadsPtr = nullptr; // pointer to threads array
    char* bufferPtr = nullptr;  // pointer to temp buffer for threads
    // try/catch block to allocate all memory and catch errors before starting threads
    try {
        // Get handle to input wave. Make sure it exists.
        wavH = p->w1;
        if (wavH == NIL) throw result = NON_EXISTENT_WAVE;
        // Get wave data type.
        waveType = WaveType(wavH);
        // Can't process text waves
        if (waveType == TEXT_WAVE_TYPE) throw result = NOTEXTWAVES;
        // Get number of used dimensions in wave.
        if (MDGetWaveDimensions(wavH, &numDimensions, dimensionSizes)) throw result = WAVEERROR_NOS;
        // Check that wave is 2D or 3D
        if ((numDimensions == 1) || (numDimensions == 4)) throw result = INPUTNEEDS_2D3D_WAVE;
        // Get dimension sizes info and do special case when number of frames is 1
        xSize = dimensionSizes[0];
        ySize = dimensionSizes[1];
        if (numDimensions == 2){
            zSize = 1;
            is3D = 0;
        }else{
            zSize = dimensionSizes[2];
            is3D = 1;
        }
        // Get the offset to the data in the wave
        if (MDAccessNumericWaveData(wavH, kMDWaveAccessMode0, &dataOffset)) throw result = WAVEERROR_NOS;
        dataStartPtr = (char*)(*wavH) + dataOffset;
        // get ready for Multi threading
        nThreads = gNumProcessors;
        if (zSize < gNumProcessors) nThreads = zSize;
        // make an array of threadPramsStruct
        paramArrayPtr = (TransposeFramesThreadParamsPtr)WMNewPtr(nThreads * sizeof(TransposeFramesThreadParams));
        if (paramArrayPtr == nullptr) throw result = MEMFAIL;
        // make an array of pthread_t
        threadsPtr = (pthread_t*)WMNewPtr(nThreads * sizeof(pthread_t));
        if (threadsPtr == nullptr) throw result = MEMFAIL;
        if (xSize != ySize){ // not square frames, so need to  make a frame sized buffer
            switch (waveType) {
                case NT_I64 | NT_UNSIGNED:
                case NT_I64:
                case NT_FP64:
                    bufferPtr = (char*)WMNewPtr (xSize * ySize * nThreads * 8);
                    break;
                case NT_I32 | NT_UNSIGNED:
                case NT_I32:
                case NT_FP32:
                    bufferPtr = (char*)WMNewPtr (xSize * ySize * nThreads * 4);
                    break;
                case NT_I16 | NT_UNSIGNED:
                case NT_I16:
                    bufferPtr = (char*)WMNewPtr (xSize * ySize * nThreads * 2);
                    break;
                case NT_I8 | NT_UNSIGNED:
                case NT_I8:
                    bufferPtr = (char*)WMNewPtr (xSize * ySize * nThreads * 1);
                    break;
                default:
                    throw result = NUMTYPE;
                    break;
            }
            if (bufferPtr == nullptr) throw result = NOMEM;
        }
    }catch (int result){ // free any memory we may have allocated so far
        if (bufferPtr != nullptr) WMDisposePtr ((Ptr)bufferPtr);
        if (threadsPtr != nullptr) WMDisposePtr ((Ptr)threadsPtr);
        if (paramArrayPtr != nullptr) WMDisposePtr ((Ptr)paramArrayPtr);
        p -> result = (double)(result - FIRST_XOP_ERR);
#ifdef NO_IGOR_ERR
        return (0);
#else
        return (result);
#endif
    }
    // fill the array of paramater structs
    for (iThread = 0; iThread < nThreads; iThread++){
        paramArrayPtr[iThread].inPutWaveType = waveType;
        paramArrayPtr[iThread].dataStartPtr = dataStartPtr;
        if (xSize != ySize) paramArrayPtr[iThread].bufferPtr = bufferPtr;
        paramArrayPtr[iThread].xSize = xSize;
        paramArrayPtr[iThread].ySize = ySize;
        paramArrayPtr[iThread].zSize = zSize;
        paramArrayPtr[iThread].ti=iThread; // number of this thread, starting from 0
        paramArrayPtr[iThread].tN =nThreads; // total number of threads
    }
    // create the threads
    for (iThread = 0; iThread < nThreads; iThread++){
        pthread_create (&threadsPtr[iThread], NULL, TransposeFramesThread, (void *) &paramArrayPtr[iThread]);
    }
    // Wait till all the threads are finished
    for (iThread = 0; iThread < nThreads; iThread++){
        pthread_join (threadsPtr[iThread], NULL);
    }
   // stuff for un-square frames
    if (xSize != ySize){
        WMDisposePtr ((Ptr)bufferPtr);
        // redimension, swapping X and Ys
        dimensionSizes [ROWS] = ySize;
        dimensionSizes [COLUMNS] = xSize;
        if (is3D){
            dimensionSizes [LAYERS] = zSize;
        }else{
            dimensionSizes [LAYERS] =0;
        }
        dimensionSizes [CHUNKS] =0;
        MDChangeWave2 (wavH, -1, dimensionSizes, 1);
    }
    // free memory for pThreads Array
    WMDisposePtr ((Ptr)threadsPtr);
    // Free paramaterArray memory
    WMDisposePtr ((Ptr)paramArrayPtr);
    // Inform Igor that we have changed the input wave.
    WaveHandleModified(wavH);
    p -> result = (0);
    return (0);
}

/* -------------------------------Decumulate Functions--------------------------------------------------
 For photon counting, the counter keeps a running total; i.e., it accumulates. To get counts for each pixel,
 we need to subtract from the from count at each pixel the count of the pixel before it, i.e., decumulate
 Also, the first pixel in each line contains photons from flyback time, so set it to 0
 -------------------------------------------------------------------------------------------------------- */


/* The following template is used to handle any one of the different types of wave data
 Last Modified 2013 by Jamie Boyd */
template <typename T> void DecumulateT(T *srcWaveStart, CountInt NumPnts, UInt32 epMax, UInt32 maxCount){
    T *srcWavePtr;
    // set srcWavePtr to last point in wave and work backwards
    for (srcWavePtr = srcWaveStart + NumPnts; srcWavePtr > srcWaveStart; srcWavePtr --){
        *srcWavePtr -= *(srcWavePtr -1);
    }
}


/****************************************************************************************************************
Decumulate.
typedef struct DecumulateParams
double bitSize;   //bitsize of the counter. expected to be either 24 or 32
waveHndl w1;  // input wave - is overwritten
Last Modified 2025/06/13 by Jamie Boyd
****************************************************************************************************************/
int Decumulate (DecumulateParamsPtr p) {

    int result = 0;    // The error returned from various Wavemetrics functions
    int waveType; //  Wavetypes numeric codes for things like 32 bit floating point, 16 bit int, etc
    int numDimensions;    // number of dimensions in input and output waves
    CountInt dimensionSizes[MAX_DIMENSIONS+1];    // an array used to hold the width, height, layers, and chunk sizes
    BCInt dataOffset;    //offset in bytes from begnning of handle to a wave to the actual data - size of headers, units, etc.
    waveHndl wavH;        // handle to the input wave
    CountInt numPnts;    // number of points in the wave
    char* srcWaveStart;    // Pointer to start of data in input wave. Need to use char for these to use WM function to get data offset
    UInt32 epMax;        // expected maximum counts per pixel
    UInt32 maxCnt;        // maximum value of counter before rollover

    try{
        // Get handle to input wave. Make sure it exists.
        wavH = p->w1;
        if (wavH == NIL) throw result = NON_EXISTENT_WAVE;
        //Get data type, no text waves
        waveType = WaveType(p->w1);
        if (waveType==TEXT_WAVE_TYPE) throw result = NOTEXTWAVES;
        // Get number of used dimensions in wave.
        if (MDGetWaveDimensions(wavH, &numDimensions, dimensionSizes)) throw result = WAVEERROR_NOS;
        // Get the offsets to the start of the data in the wave
        if (MDAccessNumericWaveData(wavH, kMDWaveAccessMode0, &dataOffset)) throw result=WAVEERROR_NOS;
        //make pointer to stat of data
        srcWaveStart = (char*)(*wavH) + dataOffset;
        //Figure out how many points are in the wave
        numPnts =WavePoints(wavH);
        //get expected maximum values per pixel
        epMax = p->expMax;
        maxCnt = pow (2, p->bitSize);
        //Do the decumulating
        switch(waveType){
        case NT_I8:
            DecumulateT((char*)srcWaveStart,numPnts,epMax, maxCnt);
            break;
        case (NT_I8 | NT_UNSIGNED):
            DecumulateT((unsigned char*)srcWaveStart,numPnts,epMax, maxCnt);
            break;
        case NT_I16:
            DecumulateT((short*)srcWaveStart,numPnts,epMax, maxCnt);
            break;
        case (NT_I16 | NT_UNSIGNED):
            DecumulateT((unsigned short*)srcWaveStart,numPnts,epMax, maxCnt);
            break;
        case NT_I32:
            DecumulateT((long*)srcWaveStart,numPnts,epMax, maxCnt);
            break;
        case (NT_I32| NT_UNSIGNED):
            DecumulateT((unsigned long*)srcWaveStart,numPnts,epMax, maxCnt);
            break;
        case NT_FP32:
            DecumulateT((float*)srcWaveStart,numPnts,epMax, maxCnt);
            break;
        case NT_FP64:
            DecumulateT((double*)srcWaveStart, numPnts,epMax, maxCnt);
            break;
        default:    // Unknown data type - possible in a future version of Igor.
            throw result = NT_FNOT_AVAIL;
            break;
        }
        WaveHandleModified(wavH);            // Inform Igor that we have changed the input wave.
        p->result = result;        // return 0 for success
        return result;
    }catch (int result){
        p -> result = result;
        return result;  // XFUNC error code.
    }
}
