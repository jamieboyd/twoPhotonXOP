#include "twoPhoton.h"
#include <math.h>

/* ------------------------------LSM Utilities --------------------------------------------------
utility functions specialized for Laser Scanning Microscope data acquisition
 Last Modified 2026/08/11 by Jamie Boyd - changed long to SInt32 and unsigned long to SInt32 and added support for 64 bit Integer waves and added FastIntCopy
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
    CountInt toNextStart = (2 * lineLen) - (lineLen/2);    //due to rounding of integer math, it doesn't matter if LineLen is even or odd
    CountInt toNextEnd = (2 * lineLen) + (lineLen/2);
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
Last Modified 2026/08/31 by Jamie Boyd */
typedef struct SwapEvenThreadParams{
    int inPutWaveType;
    char* dataStartPtr;     // pointer to start of data in the wave
    CountInt startOffset;     // point offset to the line this thread starts on
    CountInt numLines;      // number of lines for this thread to do
    CountInt lineLen;       // number of points in each line
}SwapEvenThreadParams, *SwapEvenThreadParamsPtr;


/* Each thread to swap a range of rows starts with this function
 Last Modified 2026/08/31 by Jamie Boyd */
void* SwapEvenThread (void* threadarg){
    struct SwapEvenThreadParams* p;
    p = (struct SwapEvenThreadParams*) threadarg;
    // call the right template function for the wave types
    switch (p->inPutWaveType) {
        case NT_I8:
            SwapEvenT ((char*) p->dataStartPtr + p->startOffset, p->numLines, p->lineLen);
            break;
        case (NT_I8 | NT_UNSIGNED):
            SwapEvenT ((unsigned char*) p->dataStartPtr + p->startOffset, p->numLines, p->lineLen);
            break;
        case NT_I16:
            SwapEvenT ((short*)p->dataStartPtr + p->startOffset, p->numLines, p->lineLen);
            break;
        case (NT_I16 | NT_UNSIGNED):
            SwapEvenT ((unsigned short*)p->dataStartPtr + p->startOffset, p->numLines, p->lineLen);
            break;
        case NT_I32:
            SwapEvenT ((SInt32*)p->dataStartPtr + p->startOffset, p->numLines, p->lineLen);
            break;
        case (NT_I32| NT_UNSIGNED):
            SwapEvenT ((UInt32*)p->dataStartPtr + p->startOffset, p->numLines, p->lineLen);
            break;
        case NT_I64:
            SwapEvenT((SInt64*)p->dataStartPtr + p->startOffset, p->numLines, p->lineLen);
            break;
        case (NT_I64 | NT_UNSIGNED):
            SwapEvenT((UInt64*)p->dataStartPtr + p->startOffset, p->numLines, p->lineLen);
            break;
        case NT_FP32:
            SwapEvenT ((float*)p->dataStartPtr + p->startOffset, p->numLines, p->lineLen);
            break;
        case NT_FP64:
            SwapEvenT ((double*)p->dataStartPtr + p->startOffset, p->numLines, p->lineLen);
            break;
    }
    return nullptr;
}

/* SwapEven XOP entry function
 SwapEvenParams:
 wave handle to start of data of input wave, which is overwrittten
 result =  0 or error code
 Last Modified 2026/08/31 by Jamie Boyd */
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
    int iThread, nThreads;
    CountInt linesPerThread, lastThreadLines;
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
        // Get the offset to the data in the wave
        if (MDAccessNumericWaveData(wavH, kMDWaveAccessMode0, &dataOffset)) throw result = WAVEERROR_NOS;
        dataStartPtr = (char*)(*wavH) + dataOffset;
        // Multi threading. Each thread must get an even number of lines and should have at least 10 lines, as a ball park figure
        nThreads = gNumProcessors;
        linesPerThread = numLines/nThreads;
        if (linesPerThread % 2)  linesPerThread --;
        if (linesPerThread < 10){
            nThreads = (int)(numLines/10);
            if (nThreads == 0){     // because of integer truncation if numLines < 10
                nThreads = 1;
                linesPerThread = numLines;
            } else{
                linesPerThread = numLines/nThreads;
                if (linesPerThread % 2) linesPerThread --;
            }
        }
        lastThreadLines = numLines % nThreads; // last thread gets any left over lines
        paramArrayPtr = (SwapEvenThreadParamsPtr)WMNewPtr(nThreads * sizeof(SwapEvenThreadParams));
        if (paramArrayPtr == nullptr) throw result = MEMFAIL;
        // make an array of pthread_t
        threadsPtr =(pthread_t*)WMNewPtr(nThreads * sizeof(pthread_t));
        if (threadsPtr == nullptr) throw result = MEMFAIL;
    }catch (int result){
        if (paramArrayPtr != nullptr) WMDisposePtr ((Ptr)paramArrayPtr);
        if (threadsPtr != nullptr) WMDisposePtr ((Ptr)threadsPtr);
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
        paramArrayPtr[iThread].startOffset = iThread * linesPerThread * lineLen;
        paramArrayPtr[iThread].numLines = linesPerThread;
        if (iThread == nThreads -1) paramArrayPtr[iThread].numLines += lastThreadLines;
        paramArrayPtr[iThread].lineLen = lineLen;
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

/* ------------------------------ SwapEvenReverseEven ------------------------------------------------------
 horizontally swaps every other line in an image or series of images plus vertically reverses every other frame
 Used after doing double back and forth scanning, where every other line is scanned from the opposite direction
 and every other frame is scanned from the oppoite direction in Y
 ----------------------------------------------------------------------------------------------------- */


//preprocessor macro to swap 4 values
#define SWAP4(a,b,c,d) temp=(d);(d)=(c);(c)=(b);(b)=(a);(a)=temp;


/* Function Template for SwapEvenReverseEven function
 every other frame needs to reversed in Y as well as every other line swapped in X
 The function expects to get a pointer (dataStartPtr) to the start of an odd numbered frame, i.e., not reversed
 CoutInt numFrames is number of frames to process, it may be even or odd
 CountInt numLines is the number of lines in each frame, not total number of lines as in swapFrames,
 because we need to track frames,
 Last Modified 2026/08/27 by Jamie Boyd */
template <typename T> void SwapEvenReverseEvenT (T *dataStartPtr, CountInt numFrames, CountInt numLines, CountInt lineWidth) {
    
    // variable used for swapping
    T temp;
    // calculate frame size
    CountInt frameSize = numLines * lineWidth;
    // calculate end of data
    T *dataEndPtr = dataStartPtr + (numFrames  * frameSize);
    // half values are used for y swapping
    CountInt halfWidth = lineWidth/2;
    CountInt halfFrame = frameSize/2;
    // calculate jumps to next Swap X frame for pSL and pEL after finishing a reverse Y frame
    CountInt toNextSwapFrame = halfFrame + lineWidth;
    // calculate distance to pSU and pEU from initialized pSL and PEL for a reverse Y frame
    CountInt fromLtoU = (numLines -1) * lineWidth;
    // calculate value to add/subtract to get to next line for pSL and pEU
    CountInt toNextpSLEU = halfWidth;
    CountInt toNextpSUEL = lineWidth + halfWidth;
    // values needed for X swapping
    CountInt toNextLineStart = (2 * lineWidth) - (lineWidth/2);    //due to rounding of integer math, it doesn't matter if LineWidth is even or odd
    CountInt toNextLineEnd = (2 * lineWidth) + (lineWidth/2);
    // pointers initialized to the start and end of the line at the upper and lower bounds of the frame (the 4 corners)
    T *pSL, *pSU, *pEL, *pEU;
    // will point to end of loop through a frame of data, when Y swapping it only needs to be half a frame
    T *frameEnd;
    for ( pSL = dataStartPtr + lineWidth - toNextSwapFrame, pEL = dataStartPtr + 2*lineWidth - 1 - toNextSwapFrame; ; ){
        // start of a swap X frame, outer loop goes by line by line, inner loop swaps values for a single line
        for (pSL += toNextSwapFrame, pEL += toNextSwapFrame, frameEnd = pSL + frameSize ; pSL < frameEnd && pSL < dataEndPtr ; pSL += toNextLineStart, pEL += toNextLineEnd){
            for (; pSL < pEL ; pSL++, pEL-- ){
                SWAP (*pSL, *pEL);
            }
        } // end of a swap X frame
        if (pSL > dataEndPtr) break; // number of frames may be odd, so we may end on a swap X frame
        // start of a Reverse Y frame
        for (pSL -= lineWidth,  pSU = pSL + fromLtoU, pEL -= lineWidth, pEU = pEL + fromLtoU, frameEnd = pSL + halfFrame ; pSL < frameEnd  ; pSL += toNextpSLEU, pSU -= toNextpSUEL, pEL += toNextpSUEL, pEU -= toNextpSLEU){
            // Y swap 2 lines where lower line is X swapped and upper line is not
            for (; pSL < pEL; pSL +=1, pSU +=1, pEL -=1, pEU -=1){
                SWAP4(*pSL, *pSU, *pEL, *pEU);
            }
            pSL += toNextpSLEU;
            pSU -= toNextpSUEL;
            pEL += toNextpSUEL;
            pEU -= toNextpSLEU;
            // if number of lines is not divisible by 4, we may have run out of lines, so check
            if (pSL >= frameEnd) break;
            // Y swap 2 lines where lower line is not X swapped and upper line is X swapped
            for (; pSL < pEL; pSL +=1, pSU +=1, pEL -=1, pEU -=1){
                SWAP4(*pSL, *pEU, *pEL, *pSU);
            }
        } // end of a Reverse Y frame
    } // end of frames
}
 
/* Structure to pass data to each SwapEvenReverseEven Thread
Last Modified 2026/08/31 by Jamie Boyd */
typedef struct SwapEvenReverseEvenThreadParams{
    int inPutWaveType;          // data type (integer, signed unsigned, floating point) for wave
    char* dataStartPtr;         // where this thread starts processing data
    CountInt startOffset;;        // point offset to frame where this thread starts processing data
    CountInt numFrames;         // number of frames for this thread to process
    CountInt numLines;          // number of lines in a frame (Y size)
    CountInt lineWidth;         // number of pixels in a line (X size)
}SwapEvenReverseEvenThreadParams, *SwapEvenReverseEvenThreadParamsPtr;


/* Each thread to swap a range of rows starts with this function
 Last Modified 2026/08/27 by Jamie Boyd */
void* SwapEvenReverseEvenThread (void* threadarg){
    struct SwapEvenReverseEvenThreadParams* p;
    p = (struct SwapEvenReverseEvenThreadParams*) threadarg;
    // call the right template function for the wave types

    switch (p->inPutWaveType) {
    case NT_I8:
        SwapEvenReverseEvenT ((char*) p->dataStartPtr + p->startOffset, p->numFrames, p->numLines, p->lineWidth);
        break;
    case (NT_I8 | NT_UNSIGNED):
            SwapEvenReverseEvenT ((unsigned char*)p->dataStartPtr + p->startOffset, p->numFrames, p->numLines, p->lineWidth);
        break;
    case NT_I16:
            SwapEvenReverseEvenT ((short*)p->dataStartPtr + p->startOffset, p->numFrames, p->numLines, p->lineWidth);
        break;
    case (NT_I16 | NT_UNSIGNED):
            SwapEvenReverseEvenT ((unsigned short*)p->dataStartPtr + p->startOffset, p->numFrames, p->numLines, p->lineWidth);
        break;
    case NT_I32:
            SwapEvenReverseEvenT ((SInt32*)p->dataStartPtr + p->startOffset, p->numFrames, p->numLines, p->lineWidth);
        break;
    case (NT_I32| NT_UNSIGNED):
            SwapEvenReverseEvenT ((UInt32*)p->dataStartPtr + p->startOffset, p->numFrames, p->numLines, p->lineWidth);
        break;
    case NT_I64:
            SwapEvenReverseEvenT((SInt64*)p->dataStartPtr + p->startOffset, p->numFrames, p->numLines, p->lineWidth);
        break;
    case (NT_I64 | NT_UNSIGNED):
            SwapEvenReverseEvenT((UInt64*)p->dataStartPtr + p->startOffset, p->numFrames, p->numLines, p->lineWidth);
        break;
    case NT_FP32:
            SwapEvenReverseEvenT ((float*)p->dataStartPtr + p->startOffset, p->numFrames, p->numLines, p->lineWidth);
        break;
    case NT_FP64:
            SwapEvenReverseEvenT ((double*)p->dataStartPtr + p->startOffset, p->numFrames, p->numLines, p->lineWidth);
        break;
    }
    return nullptr;
}

/* SwapEvenReverseEven XOP entry function
 SwapEvenParams:
 wave handle to start of data of input wave, which is overwrittten. For "double-turbo" mode where the Y galvo also reverses
 direction instead of "flying back" at the end of the frame. Every other frame is collected in reverse order. We need to swap
 even-numbered lines and reverse even-numbered frames.
 
 result =  0 or error code
 Last Modified 2026/08/27 by Jamie Boyd */
extern "C" int SwapEvenReverseEven (SwapEvenParamsPtr p){
    waveHndl wavH = nullptr;        // handle to the input wave
    int waveType; //  Wavetypes numeric codes for things like 32 bit floating point, 16 bit int, etc
    int numDimensions;    // number of dimensions in input and output waves
    CountInt dimensionSizes[MAX_DIMENSIONS+1];    // an array used to hold the width, height, layers, and chunk sizes
    IndexInt dataOffset;    //offset in bytes from begnning of handle to a wave to the actual data - size of headers, units, etc.
    CountInt lineWidth;    // The width of each line in the image
    CountInt numLines;    // The number of lines in each frame
    CountInt numFrames;     // the number of frames in the wave
    int result;    // The error returned from various Wavemetrics functions
    char* dataStartPtr;    // Pointer to start of data in input wave. Need to use char for these to use WM function to get data offset
    // Multi threading
    int iThread, nThreads;
    CountInt framesPerThread, lastThreadFrames;
    SwapEvenReverseEvenThreadParamsPtr paramArrayPtr = nullptr;
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
        // Check that wave1 is 3D
        if (numDimensions != 3) throw result = INPUTNEEDS_3D_WAVE;
        // Get dimension size info and calculate number of lines to process
        lineWidth = dimensionSizes[0];
        numLines = dimensionSizes[1];
        numFrames = dimensionSizes[2];
        // Get the offsets to the data in the wave
        if (MDAccessNumericWaveData(wavH, kMDWaveAccessMode0, &dataOffset)) throw result = WAVEERROR_NOS;
        dataStartPtr = (char*)(*wavH) + dataOffset;
        // Multi threading
        nThreads = gNumProcessors;
        // each thread must get an even number of frames (except possibly last frame)
        framesPerThread = numFrames/nThreads;
        if (framesPerThread % 2) framesPerThread --;
        if (framesPerThread < 2){
            nThreads = (int)(numFrames/2);
            if (nThreads == 0){
                nThreads = 1;
                framesPerThread = numFrames;
            } else{
                framesPerThread = numFrames/nThreads;
                if (framesPerThread % 2) framesPerThread --;
            }
        }
        lastThreadFrames = numFrames % nThreads; // last thread gets any leftovere frames
        paramArrayPtr = (SwapEvenReverseEvenThreadParamsPtr)WMNewPtr(nThreads * sizeof(SwapEvenReverseEvenThreadParams));
        if (paramArrayPtr == nullptr) throw result = MEMFAIL;
        // make an array of pthread_t
        threadsPtr =(pthread_t*)WMNewPtr(nThreads * sizeof(pthread_t));
        if (threadsPtr == nullptr) throw result = MEMFAIL;
    }catch (int result){
        if (paramArrayPtr != nullptr) WMDisposePtr ((Ptr)paramArrayPtr);
        if (threadsPtr != nullptr) WMDisposePtr ((Ptr)threadsPtr);
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
        paramArrayPtr[iThread].startOffset =  iThread * framesPerThread * numLines * lineWidth;
        paramArrayPtr[iThread].numFrames = framesPerThread;
        if (iThread == nThreads -1) paramArrayPtr[iThread].numFrames += lastThreadFrames;
        paramArrayPtr[iThread].numLines = numLines;
        paramArrayPtr[iThread].lineWidth = lineWidth;
    }
    // create the threads
    for (iThread = 0; iThread < nThreads; iThread++){
        pthread_create (&threadsPtr[iThread], NULL, SwapEvenReverseEvenThread, (void *) &paramArrayPtr[iThread]);
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


/* ------------------------------ ReverseEven ------------------------------------------------------
 horizontally swaps every other line in an image while vertically reversing the image
 Used when doing double back and forth scannings and canning continuously, but processing
 frame by frame, so odd number frames get SwapEven and even numbered frames get ReverseEven
 ----------------------------------------------------------------------------------------------------- */

/* Function Template for ReverseEven function
 When a single frame needs to be reversed in Y as well as every other line swapped in X
 The function expects to get a pointer (dataStartPtr) to the start of an odd numbered frame, i.e., not reversed
 CoutInt numFrames is number of frames to process, it may be even or odd
 CountInt numLines is the number of lines in each frame, not total number of lines as in swapFrames,
 because we need to track frames,
 Last Modified 2026/08/27 by Jamie Boyd */
template <typename T> void ReverseEvenT (T *dataStartPtr, CountInt numLines, CountInt lineWidth) {
    
    // variable used for swapping
    T temp;
    CountInt frameSize = numLines * lineWidth;     // size of the image frame
    // half values are used for y swapping
    CountInt halfWidth = lineWidth/2;
    T *halfFrame = dataStartPtr + frameSize/2;
    // calculate value to add/subtract to get to next line for pSL and pEU
    CountInt toNextpSLEU = halfWidth;
    CountInt toNextpSUEL = lineWidth + halfWidth;
    
    // pointers initialized to the start and end of the line at the upper and lower bounds of the frame (the 4 corners)
    T *pSL, *pSU, *pEL, *pEU;
    pSL = dataStartPtr;
    pEL  = dataStartPtr + lineWidth -1;
    pSU = dataStartPtr + frameSize - lineWidth;
    pEU = dataStartPtr + frameSize -1;
    for (;  pSL < halfFrame  ; pSL += toNextpSLEU, pSU -= toNextpSUEL, pEL += toNextpSUEL, pEU -= toNextpSLEU){
        for (; pSL < pEL; pSL +=1, pSU +=1, pEL -=1, pEU -=1){
            SWAP4(*pSL, *pSU, *pEL, *pEU);
        }
        pSL += toNextpSLEU;
        pEU -= toNextpSLEU;
        pSU -= toNextpSUEL;
        pEL += toNextpSUEL;
        // if number of lines is not divisible by 4, we may have run out of lines, so check
        if (pSL >= halfFrame) break ;
        for (; pSL < pEL; pSL +=1, pSU +=1, pEL -=1, pEU -=1){
            SWAP4(*pSL, *pEU, *pEL, *pSU);
        }
    }
}

/* ReverseEven XOP entry function
 SwapEvenParams:
 wave handle to start of data of input wave, which is overwrittten. For "double-turbo" mode where the Y galvo also reverses
 direction instead of "flying back" at the end of the frame. Every other frame is collected in reverse order. We need to swap
 even-numbered lines and reverse even-numbered frames.
 
 result =  0 or error code
 Last Modified 2026/08/27 by Jamie Boyd */
extern "C" int ReverseEven (SwapEvenParamsPtr p){
    waveHndl wavH = nullptr;        // handle to the input wave
    int waveType; //  Wavetypes numeric codes for things like 32 bit floating point, 16 bit int, etc
    int numDimensions;    // number of dimensions in input and output waves
    CountInt dimensionSizes[MAX_DIMENSIONS+1];    // an array used to hold the width, height, layers, and chunk sizes
    IndexInt dataOffset;    //offset in bytes from begnning of handle to a wave to the actual data - size of headers, units, etc.
    CountInt lineWidth;    // The width of each line in the image
    CountInt numLines;    // The number of lines in each frame
    int result;    // The error returned from various Wavemetrics functions
    char* dataStartPtr;    // Pointer to start of data in input wave. Need to use char for these to use WM function to get data offset
    // Multi threading ? - one frame, so no threads
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
        // Check that wave1 is 2D
        if (numDimensions != 2) throw result = INPUTNEEDS_2D_WAVE;
        // Get dimension size info and calculate number of lines to process
        lineWidth = dimensionSizes[0];
        numLines = dimensionSizes[1];
        // Get the offset to the data in the wave
        if (MDAccessNumericWaveData(wavH, kMDWaveAccessMode0, &dataOffset)) throw result = WAVEERROR_NOS;
        dataStartPtr = (char*)(*wavH) + dataOffset;
    }catch (int result){
        p -> result = (double)(result - FIRST_XOP_ERR);
#ifdef NO_IGOR_ERR
        return (0);
#else
        return (result);
#endif
    }
    switch (waveType) {
        case NT_I8:
            ReverseEvenT ((char*) dataStartPtr, numLines, lineWidth);
            break;
        case (NT_I8 | NT_UNSIGNED):
            ReverseEvenT ((unsigned char*)dataStartPtr, numLines, lineWidth);
            break;
        case NT_I16:
            ReverseEvenT ((short*)dataStartPtr, numLines, lineWidth);
            break;
        case (NT_I16 | NT_UNSIGNED):
            ReverseEvenT ((unsigned short*)dataStartPtr, numLines, lineWidth);
            break;
        case NT_I32:
            ReverseEvenT ((SInt32*)dataStartPtr, numLines, lineWidth);
            break;
        case (NT_I32| NT_UNSIGNED):
            ReverseEvenT ((UInt32*)dataStartPtr, numLines, lineWidth);
            break;
        case NT_I64:
            ReverseEvenT((SInt64*)dataStartPtr, numLines, lineWidth);
            break;
        case (NT_I64 | NT_UNSIGNED):
            ReverseEvenT((UInt64*)dataStartPtr, numLines, lineWidth);
            break;
        case NT_FP32:
            ReverseEvenT ((float*)dataStartPtr, numLines, lineWidth);
            break;
        case NT_FP64:
            ReverseEvenT ((double*)dataStartPtr, numLines, lineWidth);
            break;
    }
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

        case (NT_I64 | NT_UNSIGNED):
            switch (DSType) {
            case 1:    //average
                DSaverageT((UInt64*)dataStartPtr + startPos, tPoints, boxFactor);
                break;
            case 2: // sum
                DSsumUnSignedIntT((UInt64*)dataStartPtr + startPos, tPoints, boxFactor, (UInt64)ULLONG_MAX);
                break;
            case 3: //max
                DSMaxT((UInt64*)dataStartPtr + startPos, tPoints, boxFactor);
                break;
            case 4: // median
                DSMedianT((UInt64*)dataStartPtr + startPos, tPoints, boxFactor);
                break;
            }
            break;
        case NT_I64:
            switch (DSType) {
            case 1:    //average
                DSaverageT((SInt64*)dataStartPtr + startPos, tPoints, boxFactor);
                break;
            case 2: // sum
                DSsumSignedIntT((SInt64*)dataStartPtr + startPos, tPoints, boxFactor, (SInt64)LLONG_MAX, (SInt64)LLONG_MIN);
                break;
            case 3: //max
                DSMaxT((SInt64*)dataStartPtr + startPos, tPoints, boxFactor);
                break;
            case 4: // median
                DSMedianT((SInt64*)dataStartPtr + startPos, tPoints, boxFactor);
                break;
            }
            break;


        case (NT_I32 | NT_UNSIGNED):
            switch (DSType){
                case 1:    //average
                    DSaverageT ((UInt32*) dataStartPtr + startPos, tPoints, boxFactor);
                    break;
                case 2: // sum
                    DSsumUnSignedIntT ((UInt32*)dataStartPtr + startPos, tPoints, boxFactor,(UInt32)ULONG_MAX);
                    break;
                case 3: //max
                    DSMaxT ((UInt32*)dataStartPtr + startPos, tPoints, boxFactor);
                    break;
                case 4: // median
                    DSMedianT ((UInt32*)dataStartPtr + startPos, tPoints, boxFactor);
                    break;
            }
            break;
        case NT_I32:
            switch (DSType){
                case 1:    //average
                    DSaverageT ((SInt32*) dataStartPtr + startPos, tPoints, boxFactor);
                    break;
                case 2: // sum
                    DSsumSignedIntT ((SInt32*)dataStartPtr + startPos, tPoints, boxFactor,(SInt32)LONG_MAX, (SInt32)LONG_MIN);
                    break;
                case 3: //max
                    DSMaxT ((SInt32*)dataStartPtr + startPos, tPoints, boxFactor);
                    break;
                case 4: // median
                    DSMedianT ((SInt32*)dataStartPtr + startPos, tPoints, boxFactor);
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
        if (wavH == nullptr) throw result = NON_EXISTENT_WAVE;
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
            zSize = 1;
        } else {
            zSize = dimensionSizes[2];
        }
        points = xSize * ySize * zSize;
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
        frameCopyStart = (T*) memcpy ((void*)frameCopyStart, (void*)dataWavePtr, frameBytes);
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
    void* dataStartPtr;         // pointer to start of input wave, which is overwritten
    void* bufferPtr;            // pointer to start of a frame-sized buffer for temporary calculations, not used for square frames
    CountInt xSize;            // number of columns in each frame
    CountInt ySize;            // number of rows in each frame
    CountInt zSize;             // number of frames in wave
    UInt8 ti; // number of this thread, starting from 0
    UInt8 tN; // total number of threads (255 "should be enough for anyone")
} TransposeFramesThreadParams, *TransposeFramesThreadParamsPtr;


/* Each thread to transpose a range of frames starts with this function
Last Modified 2025/06/27 by Jamie Boyd */
void* TransposeFramesThread (void* threadarg){
    struct TransposeFramesThreadParams* p = (struct TransposeFramesThreadParams*) threadarg;
    CountInt frameSize = p->xSize * p->ySize;
    CountInt tFrames = (p->zSize/p->tN);  // number of frames to do per thread
    CountInt startOffset = p->ti * frameSize * tFrames ;   // starting position for this thread
    if (p->ti == p->tN - 1) tFrames += (p->zSize % p->tN); // the last thread gets any left-over frames
    if (p->xSize == p->ySize){
        switch (p->inPutWaveType) {
            case NT_FP64:
                TransposeSquareFramesT (((double*) p->dataStartPtr) + startOffset, p->xSize, tFrames);
                break;
            case NT_FP32:
                TransposeSquareFramesT (((float*) p->dataStartPtr) + startOffset, p->xSize, tFrames);
                break;
            case (NT_I32 | NT_UNSIGNED):
                TransposeSquareFramesT (((UInt32*) p->dataStartPtr) + startOffset, p->xSize, tFrames);
                break;
            case NT_I32:
                TransposeSquareFramesT (((UInt32*) p->dataStartPtr) + startOffset, p->xSize, tFrames);
                break;
            case (NT_I64 | NT_UNSIGNED):
                TransposeSquareFramesT(((UInt64*)p->dataStartPtr) + startOffset, p->xSize, tFrames);
                break;
            case NT_I64:
                TransposeSquareFramesT(((SInt64*)p->dataStartPtr) + startOffset, p->xSize, tFrames);
                break;
            case (NT_I16 | NT_UNSIGNED):
                TransposeSquareFramesT (((unsigned short*) p->dataStartPtr) + startOffset, p->xSize, tFrames);
                break;
            case NT_I16:
                TransposeSquareFramesT (((short*) p->dataStartPtr) + startOffset, p->xSize, tFrames);
                break;
            case (NT_I8 | NT_UNSIGNED):
                TransposeSquareFramesT (((unsigned char*) p->dataStartPtr) + startOffset, p->xSize, tFrames);
                break;
            case NT_I8:
                TransposeSquareFramesT (((char*) p->dataStartPtr) + startOffset, p->xSize, tFrames);
                break;
        }
    }else{
        CountInt bufferOffset =p->ti * frameSize; // this is in points, not bytes.
        switch (p->inPutWaveType) {
            case NT_FP64:
                TransposeFramesT (((double*) p->dataStartPtr) + startOffset, ((double*) p->bufferPtr) + bufferOffset, p->xSize, p->ySize, tFrames);
                break;
            case NT_FP32:
                TransposeFramesT (((float*) p->dataStartPtr) + startOffset, ((float*) p->bufferPtr) + bufferOffset, p->xSize, p->ySize, tFrames);
                break;
            case (NT_I64 | NT_UNSIGNED):
                TransposeFramesT(((UInt64*)p->dataStartPtr) + startOffset, ((UInt64*)p->bufferPtr) + bufferOffset, p->xSize, p->ySize, tFrames);
                break;
            case NT_I64:
                TransposeFramesT(((SInt64*)p->dataStartPtr) + startOffset, ((SInt64*)p->bufferPtr) + bufferOffset, p->xSize, p->ySize, tFrames);
                break;

            case (NT_I32 | NT_UNSIGNED):
                TransposeFramesT (((UInt32*) p->dataStartPtr) + startOffset, ((UInt32*) p->bufferPtr) + bufferOffset, p->xSize, p->ySize, tFrames);
                break;
            case NT_I32:
                TransposeFramesT (((SInt32*) p->dataStartPtr) + startOffset, ((SInt32*)p->bufferPtr) + bufferOffset, p->xSize, p->ySize, tFrames);
                break;
            case (NT_I16 | NT_UNSIGNED):
                TransposeFramesT (((unsigned short*) p->dataStartPtr) + startOffset, ((unsigned short*)p->bufferPtr) + bufferOffset, p->xSize, p->ySize, tFrames);
                break;
            case NT_I16:
                TransposeFramesT (((SInt16*) p->dataStartPtr) + startOffset, ((SInt16*)p->bufferPtr) + bufferOffset, p->xSize, p->ySize, tFrames);
                break;
            case (NT_I8 | NT_UNSIGNED):
                TransposeFramesT (((unsigned char*) p->dataStartPtr) + startOffset, ((unsigned char*) p->bufferPtr) + bufferOffset, p->xSize, p->ySize, tFrames);
                break;
            case NT_I8:
                TransposeFramesT (((char*) p->dataStartPtr) + startOffset, ((char*)p->bufferPtr) + bufferOffset, p->xSize, p->ySize, tFrames);
                break;
        }
    }
    return nullptr;
 }

/* TransposeFrames XOP entry function
 TransposeFramesParams:
 wave handle to input wave, which is always overwritten
 result which is 0 or error code
 Last Modified 2025/06/27 by Jamie Boyd */
extern "C" int TransposeFrames (TransposeFramesParamsPtr p) {
    int result = 0;    // The error returned from various Wavemetrics functions
    int waveType; //  Wavetypes numeric codes for things like 32 bit floating point, 16 bit int, etc
    int numDimensions;    // number of dimensions in input and output waves
    CountInt dimensionSizes[MAX_DIMENSIONS+1];    // an array used to hold the width, height, layers, and chunk sizes
    CountInt dataOffset;    //offset in bytes from begnning of handle to a wave to the actual data - size of headers, units, etc.
    CountInt xSize;    // The length of each line in each frame
    CountInt ySize; // number of lines in each frame
    CountInt zSize; // number of frames in the stack
    waveHndl wavH;        // handle to the input wave
    char* dataStartPtr;    // Pointer to start of data in input wave. Need to use char for these to use WM function to get data offset
    // threading
    UInt8 iThread, nThreads;
    TransposeFramesThreadParamsPtr paramArrayPtr = nullptr;
    pthread_t* threadsPtr = nullptr; // pointer to threads array
    void* bufferPtr = nullptr;  // pointer to temp buffer for threads
    // try/catch block to allocate all memory and catch errors before starting threads
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
        // Check that wave is 2D or 3D
        if ((numDimensions == 1) || (numDimensions == 4)) throw result = INPUTNEEDS_2D3D_WAVE;
        // Get dimension sizes info and do special case when number of frames is 1
        xSize = dimensionSizes[0];
        ySize = dimensionSizes[1];
        if (numDimensions == 2){
            zSize = 1;
        }else{
            zSize = dimensionSizes[2];
        }
        // Get the offset to the data in the wave
        if (MDAccessNumericWaveData(wavH, kMDWaveAccessMode0, &dataOffset)) throw result = WAVEERROR_NOS;
        dataStartPtr =  ((char*)*wavH) + dataOffset;
        // get ready for Multi threading
        nThreads = gNumProcessors;
        if (zSize < gNumProcessors) nThreads = (UInt8)zSize;
        // make an array of threadPramsStruct
        paramArrayPtr = (TransposeFramesThreadParamsPtr)WMNewPtr(nThreads * sizeof(TransposeFramesThreadParams));
        if (paramArrayPtr == nullptr) throw result = MEMFAIL;
        // make an array of pthread_t
        threadsPtr = (pthread_t*)WMNewPtr(nThreads * sizeof(pthread_t));
        if (threadsPtr == nullptr) throw result = MEMFAIL;
        if (xSize != ySize){ // not square frames, so need to  make a frame sized buffer for each thread
            switch (waveType) {
                case NT_I64 | NT_UNSIGNED:
                    bufferPtr = WMNewPtr(xSize * ySize * nThreads * sizeof(UInt64));
                    break;
                case NT_I64:
                    bufferPtr = WMNewPtr(xSize * ySize * nThreads * sizeof(SInt64));
                    break;
                case NT_FP64:
                    bufferPtr = WMNewPtr (xSize * ySize * nThreads * sizeof(double));
                    break;
                case NT_I32 | NT_UNSIGNED:
                    bufferPtr = WMNewPtr(xSize * ySize * nThreads * sizeof(UInt32));
                    break;
                case NT_I32:
                    bufferPtr = WMNewPtr(xSize * ySize * nThreads * sizeof(SInt32));
                    break;
                case NT_FP32:
                    bufferPtr = WMNewPtr (xSize * ySize * nThreads * sizeof(float));
                    break;
                case NT_I16 | NT_UNSIGNED:
                    bufferPtr = WMNewPtr(xSize * ySize * nThreads * sizeof(UInt16));
                    break;
                case NT_I16:
                    bufferPtr = WMNewPtr(xSize * ySize * nThreads * sizeof(SInt16));
                    break;
                case NT_I8 | NT_UNSIGNED:
                    bufferPtr = WMNewPtr(xSize * ySize * nThreads * sizeof(UInt8));
                    break;
                case NT_I8:
                    bufferPtr = WMNewPtr (xSize * ySize * nThreads * sizeof(SInt8));
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
        paramArrayPtr[iThread].dataStartPtr = (void*)dataStartPtr;
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
 The National Instrument boards used by twoPhoton have 24 or 32 bit counters, but we make a templated
 version in case someone has a different use case.
 -------------------------------------------------------------------------------------------------------- */

/* If you have a 32 bit counter into an unsigned 32 bit wave, counter overflows overflows are handled automatically
Last Modified: 2025/06/27 by Jamie Boyd*/
void Decumulate32_32(UInt32* dataStart, CountInt numPnts, UInt32 firstValue) {
    UInt32* srcWavePtr;
    for (srcWavePtr = dataStart + numPnts - 1; srcWavePtr > dataStart; srcWavePtr --) {
        *srcWavePtr -= *(srcWavePtr - 1);
    } // first point in thread subtracts saved value - otherwise we would subtract a value maybe already decumulated by other thread
    *srcWavePtr -= firstValue;
}

/* The following template is used to handle any one of the other different types of wave data/counter size
combos but your wave type must be wider than your counter size.
 Last Modified 2025/06/27 by Jamie Boyd */
template <typename T> void DecumulateT (T *dataStart, CountInt NumPnts, UInt32 maxCount, UInt32 firstValue){
    T *srcWavePtr;
    // set srcWavePtr to last point in wave and work backwards
    for (srcWavePtr = dataStart + NumPnts-1; srcWavePtr > dataStart; srcWavePtr --){
        if (*(srcWavePtr - 1) > *srcWavePtr)              // counter rollover occurred
            *srcWavePtr += maxCount;
        *srcWavePtr -= *(srcWavePtr - 1);
    } // first point in thread subtracts saved value - otherwise we would subtract a value maybe already decumulated by other thread
    if (firstValue > *srcWavePtr)
        *srcWavePtr += maxCount;
    *srcWavePtr -= firstValue;
}


/* Structure to pass data to each Decumulate thread
 Last Modified 2025/06/27 by Jamie Boyd */
typedef struct DecumulateThreadParams {
    int inPutWaveType;          // WaveMetrics code for waveType
    void* dataStartPtr;         // pointer to start of input wave, which is overwritten
    CountInt numPoints;         // number of points of data in wave total
    UInt32 maxCount;            // max count before rollover, 2^counterSize -1
    UInt32 savedPrevValue;      // first point is in another thread, so needs to be preloaded
    UInt8 ti;                   // number of this thread, starting from 0
    UInt8 tN;                   // total number of threads (255 "should be enough for anyone")
} DecumulateThreadParams, * DecumulateThreadParamsPtr;

/* Each thread to decumulate a wave starts with this function
Last Modified 2025/06/28 by Jamie Boyd */
void* DecumulateThread (void* threadarg){
    struct DecumulateThreadParams* p = (struct DecumulateThreadParams*) threadarg;
    CountInt tPoints = p->numPoints/p->tN;
    CountInt startOffset = 1 + (p->ti * tPoints); // need the 1 point offset because subtracting from previous point
    if (p->ti == p->tN - 1) tPoints += (p->numPoints % p->tN); // the last thread gets any left-over points
   
    switch (p->inPutWaveType) {
        case NT_FP64:
            DecumulateT (((double*) p->dataStartPtr) + startOffset, tPoints, p->maxCount, p->savedPrevValue);
            break;
        case NT_FP32:
            DecumulateT (((float*) p->dataStartPtr) + startOffset, tPoints, p->maxCount, p->savedPrevValue);
            break;
        case NT_I64:
            DecumulateT (((SInt64*) p->dataStartPtr) + startOffset, tPoints, p->maxCount, p->savedPrevValue);
            break;
        case (NT_I64 | NT_UNSIGNED):
            DecumulateT (((UInt64*) p->dataStartPtr) + startOffset, tPoints, p->maxCount, p->savedPrevValue);
            break;
        case (NT_I32 | NT_UNSIGNED):
            if (p->maxCount == pow (2, (double)32) -1)  // 32-bit counter into 32 bit unsigned wave.
                Decumulate32_32 (((UInt32*) p->dataStartPtr) + startOffset, tPoints, p->savedPrevValue);
            else
                DecumulateT (((UInt32*) p->dataStartPtr) + startOffset, tPoints, p->maxCount, p->savedPrevValue);
            break;
        case (NT_I16):
            DecumulateT (((SInt16*) p->dataStartPtr) + startOffset, tPoints, p->maxCount, p->savedPrevValue);
            break;
        case (NT_I16 | NT_UNSIGNED):
            DecumulateT (((UInt16*) p->dataStartPtr) + startOffset, tPoints, p->maxCount, p->savedPrevValue);
            break;
        case (NT_I8):
            DecumulateT (((SInt8*) p->dataStartPtr) + startOffset, tPoints, p->maxCount, p->savedPrevValue);
            break;
        case (NT_I8 | NT_UNSIGNED):
            DecumulateT (((UInt8*) p->dataStartPtr) + startOffset, tPoints, p->maxCount, p->savedPrevValue);
            break;
    }
    return nullptr;
}


/* Decumulate XOP entry function
typedef struct DecumulateParams
double bitSize      bitsize of the counter.  either 24 or 32 for NI boards, but could be any integer <= 32
waveHndl w1         input wave - is overwritten. Could make an option to decumulate into an output wave?
Last Modified 2025/06/28 by Jamie Boyd */
int Decumulate (DecumulateParamsPtr p) {
    int result = 0;    // The error returned from various Wavemetrics functions
    int waveType; //  Wavetypes numeric codes for things like 32 bit floating point, 16 bit int, etc
    BCInt dataOffset;    //offset in bytes from begnning of handle to a wave to the actual data - size of headers, units, etc.
    waveHndl wavH;        // handle to the input wave
    CountInt numPnts;    // number of points in the wave
    char* srcWaveStart;    // Pointer to start of data in input wave. Need to use char for these to use WM function to get data offset
    UInt32 maxCount;        // maximum value of counter before rollover.
    CountInt tPoints;       // points to process per thread. We need this here because we need to know value of "edge" points
    // threading
    UInt8 iThread, nThreads;
    DecumulateThreadParamsPtr paramArrayPtr = nullptr;
    pthread_t* threadsPtr = nullptr; // pointer to threads array
    try {
        // Get handle to input wave. Make sure it exists.
        wavH = p->w1;
        if (wavH == nullptr) throw result = NON_EXISTENT_WAVE;
        //Get data type, no text waves
        waveType = WaveType(p->w1);
        if (waveType == TEXT_WAVE_TYPE) throw result = NOTEXTWAVES;
        // Get the offsets to the start of the data in the wave
        if (MDAccessNumericWaveData(wavH, kMDWaveAccessMode0, &dataOffset)) throw result = WAVEERROR_NOS;
        //make pointer to stat of data
        srcWaveStart = ((char*)(*wavH)) + dataOffset;
        // how many points are in the wave
        numPnts = WavePoints(wavH);
        maxCount = pow(2, p->bitSize) - 1;
        // threading
        nThreads = gNumProcessors;
        // make an array of threadPramsStruct
        paramArrayPtr = (DecumulateThreadParamsPtr)WMNewPtr(nThreads * sizeof(DecumulateThreadParams));
        if (paramArrayPtr == nullptr) throw result = MEMFAIL;
        // make an array of pthread_t
        threadsPtr = (pthread_t*)WMNewPtr(nThreads * sizeof(pthread_t));
        if (threadsPtr == nullptr) throw result = MEMFAIL;
        // do some checks that wave type is wider than counter size. Largest supported counter size is 32 bit.
        // Find value of points on "edge" of chunks for each thread. Needed ahead of time beacuse the first
        // point belonging to one thread subtracts the last value in another thread. Can't guarantee order of threads
        tPoints = numPnts/nThreads;
        for (iThread = 0; iThread < nThreads; iThread++) {
            paramArrayPtr[iThread].inPutWaveType = waveType;
            paramArrayPtr[iThread].dataStartPtr = (void*)srcWaveStart;
            paramArrayPtr[iThread].numPoints = numPnts; // total number of points in the wave
            paramArrayPtr[iThread].ti = iThread;
            paramArrayPtr[iThread].tN = nThreads;
            paramArrayPtr[iThread].maxCount = maxCount; // 2^counterBits -1. maximum count before counter rollover
            switch (waveType) {
                case NT_I8:
                    if (maxCount > SCHAR_MAX) throw result = NUMTYPE;
                    paramArrayPtr[iThread].savedPrevValue = *(((SInt8*)srcWaveStart) + (iThread * tPoints));
                    break;
                case (NT_I8 | NT_UNSIGNED):
                    if (maxCount > UCHAR_MAX) throw result = NUMTYPE;
                    paramArrayPtr[iThread].savedPrevValue = *(((UInt8*)srcWaveStart) + (iThread * tPoints));
                    break;
                case NT_I16:
                    if (maxCount > SHRT_MAX) throw result = NUMTYPE;
                    paramArrayPtr[iThread].savedPrevValue = *(((SInt16*)srcWaveStart) + (iThread * tPoints));
                    break;
                case (NT_I16 | NT_UNSIGNED):
                    if (maxCount > USHRT_MAX) throw result = NUMTYPE;
                    paramArrayPtr[iThread].savedPrevValue = *(((UInt16*)srcWaveStart) + (iThread * tPoints));
                    break;
                case (NT_I32):
                    paramArrayPtr[iThread].savedPrevValue = *(((SInt32*)srcWaveStart) + (iThread * tPoints));
                    break;
                case (NT_I32 | NT_UNSIGNED):
                    paramArrayPtr[iThread].savedPrevValue = *(((UInt32*)srcWaveStart) + (iThread * tPoints));
                    break;
                case (NT_I64):
                    paramArrayPtr[iThread].savedPrevValue = (UInt32)*(((SInt64*)srcWaveStart) + (iThread * tPoints));
                    break;
                case (NT_I64 | NT_UNSIGNED):
                    paramArrayPtr[iThread].savedPrevValue = (UInt32)*(((UInt64*)srcWaveStart) + (iThread * tPoints));
                    break;
                case NT_FP32:
                    if (maxCount > 16777216)  throw result = NUMTYPE;
                    paramArrayPtr[iThread].savedPrevValue = (UInt32)*(((float*)srcWaveStart) + (iThread * tPoints));
                    break;
                case NT_FP64:
                    paramArrayPtr[iThread].savedPrevValue = (UInt32)*(((double*)srcWaveStart) + (iThread * tPoints));
                    break;
                default:
                    throw result = NUMTYPE;
                    break;
            }
        }
    } catch (int result) { // free any memory we may have allocated so far
        if (threadsPtr != nullptr) WMDisposePtr((Ptr)threadsPtr);
        if (paramArrayPtr != nullptr) WMDisposePtr((Ptr)paramArrayPtr);
        p->result = (double)(result - FIRST_XOP_ERR);
#ifdef NO_IGOR_ERR
        return (0);
#else
        return (result);
#endif
    }
    // create the threads
    for (iThread = 0; iThread < nThreads; iThread++) {
        pthread_create(&threadsPtr[iThread], NULL, DecumulateThread, (void*)&paramArrayPtr[iThread]);
    }
    // Wait till all the threads are finished
    for (iThread = 0; iThread < nThreads; iThread++) {
        pthread_join(threadsPtr[iThread], NULL);
    }
    WaveHandleModified(wavH);            // Inform Igor that we have changed the input wave.
    p->result = (0);        // return 0 for success
    return (0);
 }

/* *********************************** FastIntCopy ***********************************************************
 Copies data between 2 Integer waves. Unsigned/signed mismatch is handled by either rotating
 the values in the dest wave, or by clipping them.
 
 For rotating, all of the data values are preserved, but an offset is added, either positive
 for copying unsigned to signed, or negative for copying from signed to unsigned.
 
 For clipping, if copying data from a signed src wave to an unsigned dest wave,
 values below zero are clipped to zero, but values above zero retain the same absolute
 value as in the signed src wave. If copying data from an unsigned wave to a signed wave,
 values above half of max value are clipped to half max value, but values below half max value
 retain the same absolute value as in the unsigned src wave.
 
 **********************************************************************************************
 Structure to pass data to any of the FastIntCopyThread variants
 Last Modified 2026/08/07 by Jamie Boyd */
typedef struct FastIntCopyThreadParams{
    char* srcDataStartPtr;      // pointer to start of data in src wave
    IndexInt srcPointOffset;     // offset to where we start copying data from
    char* destDataStartPtr;     // pointer to start of data in dest wave
    IndexInt destPointOffset;    // offset to where we start copying data to
    int srcWaveType;            // dest wave type can be inferred from srcWaveType
    CountInt numCopyPoints;     // number of points to copy
    UInt8 ti;                   // number of this thread, starting from 0
    UInt8 tN;                   // total number of threads
} FastIntCopyThreadParams, *FastIntCopyThreadParamsPtr;


/* **********************************************************************************************
 fast integer copy when the signed/unsigned statuses of source and destination waves are the same.
 Last Modified 2026/08/07 by Jamie Boyd */
template <typename TI, typename TO> void FastIntCopyT (TI* srcPtr, TO* destPtr, CountInt numCopyPoints){
    TO* endPtr;
    for (endPtr = destPtr + numCopyPoints; destPtr < endPtr; srcPtr++ , destPtr++){
        *destPtr = *srcPtr;
    }
}

/* **********************************************************************************************
 FastIntCopyThread for straight copy when the signed/unsigned statuses of source and destination waves are the same.
 Last Modified 2026/08/08 by Jamie Boyd */
void* FastIntCopyThread (void* threadarg){
    struct FastIntCopyThreadParams* p = (struct FastIntCopyThreadParams*) threadarg;
    CountInt numCopyPoints = p->numCopyPoints;
    CountInt srcPointOffset = p->srcPointOffset;
    CountInt destPointOffset = p->destPointOffset;
    UInt8 ti = p->ti;
    UInt8 tN = p->tN;
    CountInt pixPerThread = numCopyPoints/tN;
    CountInt startPos = ti * pixPerThread; // which pixel to start this thread on depends on thread number * points per thread. ti is 0 based
    if (ti == (tN - 1)) pixPerThread += (numCopyPoints % tN); // last thread gets any extra pixels
    switch (p->srcWaveType) {
        case (NT_I8 | NT_UNSIGNED):
            FastIntCopyT ((unsigned char*)p->srcDataStartPtr + startPos + srcPointOffset, (unsigned char*)p->destDataStartPtr + startPos + destPointOffset, pixPerThread);
            break;
        case (NT_I16 | NT_UNSIGNED):
            FastIntCopyT ((unsigned short*)p->srcDataStartPtr + p->srcPointOffset + startPos, (unsigned short*) p->destDataStartPtr + p->destPointOffset + startPos, pixPerThread);
            break;
        case (NT_I32| NT_UNSIGNED):
            FastIntCopyT ((UInt32*)p->srcDataStartPtr + startPos + srcPointOffset, (UInt32*) p->destDataStartPtr + startPos + destPointOffset, pixPerThread);
            break;
        case (NT_I64| NT_UNSIGNED):
            FastIntCopyT ((UInt64*)p->srcDataStartPtr + startPos + srcPointOffset, (UInt64*) p->destDataStartPtr + startPos + destPointOffset, pixPerThread);
            break;
        case (NT_I8):
           FastIntCopyT ((char*)p->srcDataStartPtr  + startPos + srcPointOffset, (char*)p->destDataStartPtr + startPos + destPointOffset,  pixPerThread);
            break;
        case (NT_I16):
            FastIntCopyT ((short*)p->srcDataStartPtr  + startPos + srcPointOffset, (short*) p->destDataStartPtr + startPos + destPointOffset, pixPerThread);
            break;
        case (NT_I32):
            FastIntCopyT ((SInt32*)p->srcDataStartPtr + startPos + srcPointOffset, (SInt32*) p->destDataStartPtr + startPos + destPointOffset, pixPerThread);
            break;
        case (NT_I64):
            FastIntCopyT ((SInt64*)p->srcDataStartPtr + startPos + srcPointOffset, (SInt64*) p->destDataStartPtr + startPos + destPointOffset, pixPerThread);
            break;
    }
    return nullptr;
}

/* **********************************************************************************************
 fast integer copy with rotation, for when the signed/unsigned statuses of source and destination waves are opposite
 Last Modified 2026/08/07 by Jamie Boyd */
template <typename TI, typename TO> void FastIntCopyRotT (TI* srcPtr, TO* destPtr, CountInt numCopyPoints, TO(rotatePnts)){
    TO* endPtr;
    for (endPtr = destPtr + numCopyPoints; destPtr < endPtr; srcPtr++ , destPtr++){
        *destPtr = *srcPtr + rotatePnts;
    }
}

/* **********************************************************************************************
 FastIntCopyThread with rotation for when the signed/unsigned statuses of source and destination waves are opposite.
 Last Modified 2026/08/08 by Jamie Boyd */
void* FastIntCopyRotThread (void* threadarg){
    struct FastIntCopyThreadParams* p;
    p = (struct FastIntCopyThreadParams*) threadarg;
    CountInt srcPointOffset = p->srcPointOffset;
    CountInt destPointOffset = p->destPointOffset;
    CountInt numCopyPoints = p->numCopyPoints;
    UInt8 ti = p->ti;
    UInt8 tN = p->tN;
    CountInt pixPerThread = numCopyPoints/tN;
    CountInt startPos = ti * pixPerThread; // which pixel to start this thread on depends on thread number * points per thread. ti is 0 based
    if (ti == (tN - 1)) pixPerThread += (numCopyPoints % tN); // last thread gets any extra pixels
    switch (p->srcWaveType) {
        case (NT_I8 | NT_UNSIGNED):
            FastIntCopyRotT ((unsigned char*)p->srcDataStartPtr  + startPos + srcPointOffset, (char*)p->destDataStartPtr + startPos + destPointOffset, pixPerThread, (char)SCHAR_MIN);
            break;
        case (NT_I16 | NT_UNSIGNED):
            FastIntCopyRotT ((unsigned short*)p->srcDataStartPtr  + startPos + srcPointOffset, (short*) p->destDataStartPtr + startPos + destPointOffset, pixPerThread, (short)SHRT_MIN);
            break;
        case (NT_I32| NT_UNSIGNED):
           FastIntCopyRotT ((UInt32*)p->srcDataStartPtr + startPos + srcPointOffset, (SInt32*) p->destDataStartPtr + startPos + destPointOffset, pixPerThread, (SInt32)LONG_MIN);
            break;
        case (NT_I64| NT_UNSIGNED):
            FastIntCopyRotT ((UInt64*)p->srcDataStartPtr + startPos + srcPointOffset, (SInt64*) p->destDataStartPtr + startPos + destPointOffset, pixPerThread, (SInt64)LLONG_MIN);
            break;
        case (NT_I8):
            FastIntCopyRotT ((char*)p->srcDataStartPtr  + startPos + srcPointOffset, (unsigned char*)p->destDataStartPtr + startPos + destPointOffset, pixPerThread, (unsigned char)UCHAR_MAX);
            break;
        case (NT_I16):
            FastIntCopyRotT ((short*)p->srcDataStartPtr  + startPos + srcPointOffset, (unsigned short*) p->destDataStartPtr + startPos + destPointOffset, pixPerThread, (unsigned short)USHRT_MAX);
            break;
        case (NT_I32):
            FastIntCopyRotT ((SInt32*)p->srcDataStartPtr + startPos + srcPointOffset, (UInt32*) p->destDataStartPtr + startPos + destPointOffset, pixPerThread, (UInt32)LONG_MAX);
            break;
        case (NT_I64):
            FastIntCopyRotT ((SInt64*)p->srcDataStartPtr + startPos, (UInt64*) p->destDataStartPtr + startPos, pixPerThread, (UInt64)LLONG_MAX);
            break;
    }
    return nullptr;
}

/* **********************************************************************************************
 fast integer copy with clipping of negative values to zero, for when the source is signed and the destination is unsigned
Last Modified 2026/08/07 by Jamie Boyd */
template <typename TI, typename TO> void FastIntCopyStoUclipT (TI* srcPtr, TO* destPtr, CountInt numCopyPoints, TO(maxSigned)){
    TO* endPtr;
    for (endPtr = destPtr + numCopyPoints; destPtr < endPtr; srcPtr++ , destPtr++){
        *destPtr = *srcPtr;
        if (*destPtr > maxSigned){
            *destPtr  = 0;
        }
    }
}

/* **********************************************************************************************
 FastIntCopyStoUclipThread Thread for fast integer copy with clipping of negative values to zero, for when the source is signed and the destination is unsigned
Last Modified 2026/08/07 by Jamie Boyd */
void* FastIntCopyStoUclipThread (void* threadarg){
    struct FastIntCopyThreadParams* p;
    p = (struct FastIntCopyThreadParams*) threadarg;
    CountInt numCopyPoints = p->numCopyPoints;
    CountInt srcPointOffset = p->srcPointOffset;
    CountInt destPointOffset = p->destPointOffset;
    UInt8 ti = p->ti;
    UInt8 tN = p->tN;
    CountInt pixPerThread = numCopyPoints/tN;
    CountInt startPos = ti * pixPerThread; // which pixel to start this thread on depends on thread number * points per thread. ti is 0 based
    if (ti == (tN - 1)) pixPerThread += (numCopyPoints % tN); // last thread gets any extra pixels
    switch (p->srcWaveType) {
        case (NT_I8):
            FastIntCopyStoUclipT ((char*)p->srcDataStartPtr + startPos + srcPointOffset, (unsigned char*)p->destDataStartPtr + startPos + destPointOffset, pixPerThread, (unsigned char)SCHAR_MAX);
            break;
        case (NT_I16):
            FastIntCopyStoUclipT ((short*)p->srcDataStartPtr  + startPos + srcPointOffset, (unsigned short*) p->destDataStartPtr + startPos + destPointOffset, pixPerThread, (unsigned short)SHRT_MAX);
            break;
        case (NT_I32):
            FastIntCopyStoUclipT ((SInt32*)p->srcDataStartPtr + startPos + srcPointOffset, (UInt32*) p->destDataStartPtr + startPos + destPointOffset, pixPerThread, (UInt32)LONG_MAX);
            break;
        case (NT_I64):
              FastIntCopyStoUclipT ((SInt64*)p->srcDataStartPtr + startPos + srcPointOffset, (UInt64*) p->destDataStartPtr + startPos + destPointOffset, pixPerThread, (UInt64)LLONG_MAX);
              break;
    }
    return nullptr;
}


/* **********************************************************************************************
 template fast integer copy with clipping of large values to the maximum signed value, for when the source is unsigned and the destination is signed
Last Modified 2026/08/07 by Jamie Boyd */
template <typename TI, typename TO> void FastIntCopyUtoSclipT (TI* srcPtr, TO* destPtr, CountInt numCopyPoints, TO(maxSigned)){
    TO* endPtr;
    for (endPtr = destPtr + numCopyPoints; destPtr < endPtr; srcPtr++ , destPtr++){
        *destPtr = *srcPtr;
        if (*destPtr < 0){
            *destPtr  = maxSigned;
        }
    }
}

/* **********************************************************************************************
 thread for fast integer copy with clipping of large values to the maximum signed value, when the source is unsigned and the destination is signed
Last Modified 2026/08/07 by Jamie Boyd */
void* FastIntCopyUtoSclipThread (void* threadarg){
    struct FastIntCopyThreadParams* p;
    p = (struct FastIntCopyThreadParams*) threadarg;
    CountInt numCopyPoints = p->numCopyPoints;
    CountInt srcPointOffset = p->srcPointOffset;
    CountInt destPointOffset = p->destPointOffset;
    UInt8 ti = p->ti;
    UInt8 tN = p->tN;
    
    CountInt pixPerThread = numCopyPoints/tN;
    CountInt startPos = ti * pixPerThread; // which pixel to start this thread on depends on thread number * points per thread. ti is 0 based
    if (ti == (tN - 1)) pixPerThread += (numCopyPoints % tN); // last thread gets any extra pixels
    switch (p->srcWaveType) {
        case (NT_I8 | NT_UNSIGNED):
            FastIntCopyUtoSclipT ((unsigned char*)p->srcDataStartPtr  + startPos + srcPointOffset, (char*)p->destDataStartPtr + startPos + destPointOffset, pixPerThread, (char)SCHAR_MAX);
            break;
        case (NT_I16 | NT_UNSIGNED):
           FastIntCopyUtoSclipT ((unsigned short*)p->srcDataStartPtr  + startPos + srcPointOffset, (short*) p->destDataStartPtr + startPos + destPointOffset, pixPerThread, (short)SHRT_MAX);
            break;
        case (NT_I32| NT_UNSIGNED):
            FastIntCopyUtoSclipT ((UInt32*)p->srcDataStartPtr + startPos + srcPointOffset, (SInt32*) p->destDataStartPtr + startPos + destPointOffset, pixPerThread, (SInt32)LONG_MAX);
            break;
        case (NT_I64| NT_UNSIGNED):
            FastIntCopyUtoSclipT ((UInt64*)p->srcDataStartPtr + startPos + srcPointOffset, (SInt64*) p->destDataStartPtr + startPos + destPointOffset, pixPerThread, (SInt64)LLONG_MAX);
            break;
    }
    return nullptr;
}

/* **************************************** FastIntCopy XOP entry function **************************************************
 waveHndl srcWaveH;
 double srcOffset;
 waveHndl destWaveH;
 double destOffset;
 double numCopyPoints;
 double ClipNotRotate;
 Last Modified 2026/09 by Jamie Boyd */
extern "C" int FastIntCopy(FastIntCopyParamsPtr p) {
    int result = 0;    // The error returned from various Wavemetrics functions
    waveHndl srcWaveH = NULL, destWaveH = NULL;    // Handles to the input and output waves
    IndexInt srcOffset, destOffset;     // offsets to points in src and dest wave to copy to/from
    IndexInt numCopyPoints;       // the number of points to be copied
    UInt8 ClipNotRotate; // set to request clipping for signed/unsigned mismatch between dest and src
    int srcWaveType, destWaveType;    // Wavemetrics numeric code for data type of wave
    CountInt srcOffsetToData, destOffsetToData; //offset in bytes from begnning of handle to a wave to the actual data - size of headers, units
    IndexInt numSrcPoints, numDestPoints;       // total number of points in src and dest waves
    char *srcDataStartPtr, *destDataStartPtr; // pointers to start of data, will be cast to correct type for wave
    UInt8 iThread, nThreads;                 // for threading
    pthread_t* threadsPtr = nullptr;       // array of pthreads
    FastIntCopyThreadParamsPtr paramArrayPtr = nullptr; // for passing data to threads
    UInt8 method;       // code for which fast copy thread to use
    char XOPbuffer [256]; // string used for XOPAlert for debugging
    try {
        // copy data from input struct.
        srcWaveH = p->srcWaveH;
        destWaveH = p->destWaveH;
        if (destWaveH == nullptr){
            sprintf(XOPbuffer, "The destination wave handle is bad.\r");
            XOPNotice (XOPbuffer);
        }
        srcOffset = (CountInt)p->srcOffset;
        destOffset = (CountInt)p->destOffset;
        numCopyPoints = (IndexInt)p->numCopyPoints;
        ClipNotRotate = (UInt8)p->ClipNotRotate;
        // Check handles to src and dest waves make sure they exist.
        if ((srcWaveH == nullptr) || (destWaveH == nullptr)) throw result= NON_EXISTENT_WAVE;
        // make sure offsets and number of points to copy do not exceed wave size
        numSrcPoints = WavePoints(srcWaveH);
        numDestPoints = WavePoints(destWaveH);
        // short cut case where numCopyPnts is 0
        if (numCopyPoints < 1){
            if (srcOffset == 0){
                numCopyPoints = numSrcPoints;
            }else if (destOffset == 0){
                numCopyPoints = numDestPoints;
            }
        }
        if (((srcOffset + numCopyPoints) > numSrcPoints) || ((destOffset + numCopyPoints) > numDestPoints)) throw result = INPUT_RANGE;
        // Get waves data types and check they are integers of same bit width, give or take the unsigned flag
        srcWaveType = WaveType(srcWaveH);
        destWaveType =  WaveType(destWaveH);
        if (srcWaveType==TEXT_WAVE_TYPE) throw result = NOTEXTWAVES;
        if ((srcWaveType == NT_FP32 ) || (srcWaveType == NT_FP64)) throw result = WAVEERROR_NOS;
        if (!((srcWaveType == destWaveType) || (abs (srcWaveType - destWaveType) == NT_UNSIGNED))) throw result = WAVEERROR_NOS;
        //Get offsets to wave data for src and dest
        if (MDAccessNumericWaveData(srcWaveH, kMDWaveAccessMode0, &srcOffsetToData) != 0) throw result = WAVEERROR_NOS;
        srcDataStartPtr = (char*)(*srcWaveH) + srcOffsetToData;
        if (MDAccessNumericWaveData(destWaveH, kMDWaveAccessMode0, &destOffsetToData) != 0) throw result = WAVEERROR_NOS;
        destDataStartPtr = (char*)(*destWaveH) + destOffsetToData;
        // set fastIntCopy method as appropriate for wave types and clipnotRotate choice
        method = 0;
        if (srcWaveType & NT_UNSIGNED) method += 1;
        if (destWaveType & NT_UNSIGNED) method += 2;
        if (ClipNotRotate) method += 4;
        // get ready for thread initiaializaton
        nThreads = gNumProcessors;
        // make array of parameter structures
        paramArrayPtr = (FastIntCopyThreadParamsPtr)WMNewPtr (nThreads * sizeof(FastIntCopyThreadParams));
        if (paramArrayPtr == NULL) throw result = MEMFAIL;
        // make an array of pthread_t or HANDLE
        threadsPtr =(pthread_t*)WMNewPtr(nThreads * sizeof(pthread_t));
        if (threadsPtr == NULL) throw result = MEMFAIL;
        // catch error before starting threads
    }catch (int result){
        if (paramArrayPtr != NULL) WMDisposePtr ((Ptr)paramArrayPtr);
        if (threadsPtr != NULL) WMDisposePtr ((Ptr)threadsPtr);
        p -> result = (double)(result - FIRST_XOP_ERR);
#ifdef NO_IGOR_ERR
        return (0);
#else
        return (result);
#endif
    }
    // fill paramater array
    for (iThread = 0; iThread < nThreads; iThread++){
        paramArrayPtr[iThread].srcDataStartPtr = srcDataStartPtr;
        paramArrayPtr[iThread].srcPointOffset = srcOffset;
        paramArrayPtr[iThread].destDataStartPtr = destDataStartPtr;
        paramArrayPtr[iThread].destPointOffset = destOffset;
        paramArrayPtr[iThread].srcWaveType = srcWaveType;
    paramArrayPtr[iThread].numCopyPoints = numCopyPoints;
        paramArrayPtr[iThread].ti=iThread; // number of this thread, starting from 0
        paramArrayPtr[iThread].tN =nThreads; // total number of threads
    }
    // create the threads according to needed thread function
    switch (method){
        case 0:     //no rotation is needed, where the src and dest waves have the same type. request for Clipping is ignored
        case 3:
        case 4:
        case 7:
            for (iThread = 0; iThread < nThreads; iThread++){
                pthread_create (&threadsPtr[iThread], NULL, FastIntCopyThread, (void *) &paramArrayPtr[iThread]);
            }
            break;
           
        case 1:      // No cliping wanted, but the src and dest waves have signed/unsigned mis-match, so rotation is needed
        case 2:
            for (iThread = 0; iThread < nThreads; iThread++){
                pthread_create (&threadsPtr[iThread], NULL, FastIntCopyRotThread, (void *) &paramArrayPtr[iThread]);
            }
            break;
            
        case 5:     // for copying from TI Unsigned src wave to TO signed dest wave, with clipping
            for (iThread = 0; iThread < nThreads; iThread++){
                pthread_create (&threadsPtr[iThread], NULL, FastIntCopyUtoSclipThread, (void *) &paramArrayPtr[iThread]);
            }
            break;
            
        case 6:     //for copying from TI Signed src wave to TO Unsigned dest wave, with clipping
            for (iThread = 0; iThread < nThreads; iThread++){
                pthread_create (&threadsPtr[iThread], NULL, FastIntCopyStoUclipThread, (void *) &paramArrayPtr[iThread]);
            }
            break;
    }
    // Wait till all the threads are finished
    for (iThread = 0; iThread < nThreads; iThread++){
        pthread_join (threadsPtr[iThread], NULL);
    }
    // free memory for pThreads Array and paramaterArray memory
    if (paramArrayPtr != NULL) WMDisposePtr ((Ptr)paramArrayPtr);
    if (threadsPtr != NULL) WMDisposePtr ((Ptr)threadsPtr);
    // Inform Igor that we have changed output wave
    WaveHandleModified(destWaveH);
    p -> result = (0);
    return (0);
}
