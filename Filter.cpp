#include "twoPhoton.h"
/* ----------------------------------------Filter ------------------------------------------------------------------
 Code for filtering waves with 2D convolution kernels - arbitrary, symetrical, or median
 Last Modified 2026/10/15 by Jamie Boyd - adding Windows threading, changing buffer allocation to array
 Last Modified 2026/08/11 by Jamie Boyd - changed long to SInt32 and unsigned long to SInt32 and added support for 64 bitInteger waves
--------------------------------------------------------------------------------------------------------------------*/


/* ------------------------------------ convolution with a 2D kernel ------------------------------------------
 convolves each frame in input wave with an arbitrary sized 2D kernel
 --------------------------------------------------------------------------------------------------------------*/


/* *************************** ConvolveMakeKernelTable ********************************************************
 Makes a table of weights for a 2D kernel
 Last Modified 2014/01/20 by Jamie Boyd */
float* ConvolveMakeKernelTable (float * kernel, int kWidth, int kHeight){
    int kRadH = (kHeight -1)/2;
    int kRadW = (kWidth-1)/2;
    int kSize = kWidth * kHeight;
    // allocate memory for kernel table
    float *kernelTablePtr =(float*) WMNewPtr (kSize * sizeof (float));
    if (kernelTablePtr == NULL) return NULL;
    int iConvo, convoEnd, convoToNextRow;
    int kx, ky;
    // get value for first point in kernel table
    int iKT =0;
     kernelTablePtr [iKT] = 0;
    convoToNextRow=kRadW;
    //yLoop
    for (iConvo = kSize -1; iConvo >= (kSize/2); iConvo -= convoToNextRow){
        // xLoop
        for (convoEnd = iConvo - kRadW ; iConvo >= convoEnd; iConvo -=1){
            kernelTablePtr [iKT]  += kernel [iConvo];
        }
    }
    // get values for first column in kernelTable- adding a row to kernel each time
    for (iKT =kWidth; iKT <  (kSize/2) + 1;  iKT +=kWidth, iConvo -= convoToNextRow){
        kernelTablePtr [iKT] = kernelTablePtr [iKT -kWidth];
        // convo X
        for (convoEnd = iConvo - kRadW ; iConvo >= convoEnd; iConvo -=1){
            kernelTablePtr [iKT] += kernel [iConvo];
        }
    }
    // left edge bottom - removing a row from kernel each time
    for (iConvo = kSize -1; iKT < kSize; iKT += kWidth, iConvo -= convoToNextRow){
        kernelTablePtr [iKT] = kernelTablePtr [iKT -kWidth];
        // convo X
        for (convoEnd = iConvo - kRadW ; iConvo >= convoEnd; iConvo -=1){
            kernelTablePtr [iKT] -= kernel [iConvo];
        }
    }
    // Now do each row based on column 0 of that row
    // top
    for (ky =0, iKT =1; ky <= kRadH; ky +=1, iKT +=1){
        // top left - adding a new column at left each time
        for (kx =1; kx <= kRadW; kx +=1, iKT +=1){
            kernelTablePtr [iKT] = kernelTablePtr [iKT - 1];
            for (iConvo = (kRadH - ky)* kWidth + kRadW -kx; iConvo < kSize; iConvo += kWidth){
                kernelTablePtr [iKT] += kernel [iConvo];
            }
        }
        // top right - removing a column at right eadh time
        for (kx-=1; kx > 0; kx -=1, iKT +=1){
            kernelTablePtr [iKT] = kernelTablePtr [iKT - 1];
            for (iConvo =(kRadH - ky)* kWidth + kRadW + kx ; iConvo < kSize; iConvo += kWidth){
                kernelTablePtr [iKT] -= kernel [iConvo];
            }
        }
    }
    // bottom
    for (ky -=1; ky > 0; ky -=1, iKT +=1){
        // bottom left - adding a new column at left each time
        for (kx =1; kx <= kRadW; kx +=1, iKT +=1){
            kernelTablePtr [iKT] = kernelTablePtr [iKT - 1];
            for (iConvo = kRadW -kx, convoEnd = (kRadH + ky) * kWidth ; iConvo < convoEnd; iConvo += kWidth){
                kernelTablePtr [iKT] += kernel [iConvo];
            }
        }
        // bottom right - removing a  column at right each time
        for (kx -=1; kx > 0; kx -=1, iKT +=1){
            kernelTablePtr [iKT] = kernelTablePtr [iKT - 1];
            for (iConvo = kRadW + kx, convoEnd = (kRadH + ky) * kWidth; iConvo < convoEnd; iConvo += kWidth){
                kernelTablePtr [iKT] -= kernel [iConvo];
            }
        }
    }
    return kernelTablePtr;
}


/* ******************************************** ConvolveX *****************************************************************
 Function template to convolve a single row in an image.  Doesn't need to explicitly know image Y size/position within Y or
 Kernel Y size. Just needs to know start-Y and end-Y position in the kernel. So the same function can be called for any Y
 position in an image. Note that destWave is passed by reference
 Last modified 2014/01/24 by Jamie Boyd */
template <typename TI, typename TO> void ConvolveX (TI *srcWave, TO *&destWave, float *kernel, float *kernelTable, UInt16 nKernelX, UInt16 radKernelX, CountInt nWaveX, UInt16 endKernelY, UInt16 &startKernel, CountInt &startConvo, UInt16 &ikernelTable, UInt16 &toNextKernelX, CountInt &nConvoX, CountInt &toNextConvoX){
    
    CountInt iWaveX;
    CountInt iConvo;
    double outVal;
    UInt16 iKernel;
    UInt16 endKernelX;
    
    // X Loop for LEFT
    for (iWaveX= 0; iWaveX < radKernelX; iWaveX+=1){
        // Y Loop for Covolution
        outVal = 0;
        endKernelX = startKernel + nConvoX;
        for (iKernel = startKernel, iConvo = startConvo; iKernel < endKernelY; iKernel += toNextKernelX, iConvo += toNextConvoX){
            // X Loop for Convolution
            for (; iKernel <  endKernelX; iKernel += 1, iConvo += 1)
                outVal += srcWave [iConvo] * kernel [iKernel];
            endKernelX += nKernelX;
        } // end of Y loop for Convolution
        // set output value
        *destWave = outVal/kernelTable[ikernelTable];
        // get ready for next pixel in this row of TOP LEFT
        startKernel -= 1; // start convo 1  point earlier than for last pixel
        nConvoX +=1; // one more convo point per row than for last pixel
        toNextKernelX -=1;
        toNextConvoX -=1;
        ikernelTable += 1; // move tp next position in convo table
        destWave +=1;
    } // end of X loop for LEFT
    // X Loop for CENTRE
    for (; iWaveX < (nWaveX - radKernelX); iWaveX +=1){
        outVal = 0;
        endKernelX = startKernel + nConvoX;
        for (iKernel = startKernel, iConvo = startConvo; iKernel < endKernelY; iKernel += toNextKernelX, iConvo += toNextConvoX){
            // X Loop for Convolution
            for (; iKernel < endKernelX; iKernel += 1, iConvo += 1)
                outVal += srcWave [iConvo] * kernel [iKernel];
            endKernelX += nKernelX;
        } // end of Y loop for Convolution
        // set output value
        *destWave = outVal/kernelTable[ikernelTable];
        // get ready for next pixel in this row
        destWave +=1;
        startConvo += 1;
    }	// end of X loop for CENTRE
    // X Loop for RIGHT
    ikernelTable += 1; // move to next position in kernel table
    for (; iWaveX < nWaveX; iWaveX +=1){
        // get ready for next pixel in this row of TOP RIGHT
        nConvoX -=1; // one less convo point per row than for last pixel
        toNextKernelX +=1;
        toNextConvoX +=1;
        endKernelX = startKernel + nConvoX;
        outVal = 0;
        // Y loop for convolution
        for (iKernel = startKernel, iConvo = startConvo; iKernel < endKernelY; iKernel += toNextKernelX, iConvo += toNextConvoX){
            // X Loop for Convolution
            for (; iKernel < endKernelX; iKernel += 1, iConvo += 1)
                outVal += srcWave[iConvo] * kernel [iKernel];
            endKernelX += nKernelX;
        } // end of Y loop for Convo
        // set output value
        *destWave = outVal/kernelTable[ikernelTable];
        destWave += 1;
        startConvo += 1;
        ikernelTable += 1; // move to next position in kernel table
    } // End of X Loop for RIGHT
}


/* **************************************** ConvolveT *************************************************************
 function template for convolving one wave with another and putting results in an output wave.
 Input wave can be 2 or 3D, but each plane is done as a separate 2D image.
 Last Modified 2014/01/24 by Jamie Boyd */
template <typename TI, typename TO> void ConvolveT (TI* srcWave, TO* destWave, TO* frameBuffer, CountInt nWaveX, CountInt nWaveY, CountInt nWaveZ, float* kernel, float * kernelTable, UInt16 nKernelX, UInt16 nKernelY){
    UInt16 radKernelX = (nKernelX - 1)/2; //radius of the kernel width, not including the central pixel
    UInt16 radKernelY = (nKernelY - 1)/2; //radius of the kernel height, not including the central pixel
    UInt16 nKernel = nKernelX * nKernelY;
    CountInt fSize = (nWaveX * nWaveY);  //frame size
    // iterating through kernel and kernel table
    UInt16 ikernelTable; // position in kernelTable
    UInt16 startKernel; // position of kernel at start of convolution of first pixel of the row
    CountInt startFrameKernel = (radKernelY * nKernelX) + radKernelX;
    UInt16 endKernelY; // end point in kernel (determined by starting position and number of rows
    // Convolution variables
    CountInt startConvo; // position in the src wave for the convolution for first pixel in the row to convolve)
    CountInt toNextConvoX = nWaveX - radKernelX -1; // amount to add to iConvo to get to next row for convolution
	UInt16 toNextKernelX= radKernelX; // amount to add to iKernel to get to next row of kernel
	CountInt nConvoX = radKernelX + 1; // number of X to convolve at start of each line
    // Pointers to Push
    TI *srcFramePtr=srcWave; // will point to start of each frame in input image
    TO *destPixPtr; // will point to each pixel in turn in output image
    TI *srcEndPtr = srcWave + fSize * nWaveZ; // end of data to process
    UInt8 inPlace = 0;
    if ((TI*)destWave == (TI*)srcWave){
		inPlace = 1;
        destPixPtr = frameBuffer;
    }else{
        destPixPtr = destWave;
    }
	// Loop through all frames
    CountInt iWaveY;
    for (;srcFramePtr < srcEndPtr; srcFramePtr += fSize){
        // Y Loop for TOP (first radKernelY rows)
        ikernelTable = 0; // at start of frame, kernel table pos = 0
        startKernel=startFrameKernel;
        endKernelY=nKernel;
        for (iWaveY= 0;  iWaveY < radKernelY; iWaveY+=1){
            startConvo =0;
            ConvolveX (srcFramePtr, destPixPtr, kernel, kernelTable, nKernelX, radKernelX, nWaveX, endKernelY, startKernel, startConvo, ikernelTable, toNextKernelX, nConvoX, toNextConvoX);
            startKernel -= nKernelX -radKernelX;
        }
        // Y Loop for MIDDLE (up to nWaveY - radKernelY rows)
        startConvo =0;
        for (; iWaveY < (nWaveY - radKernelY); iWaveY +=1){
            startKernel = radKernelX;
            ConvolveX (srcFramePtr, destPixPtr, kernel, kernelTable, nKernelX, radKernelX, nWaveX, endKernelY, startKernel, startConvo, ikernelTable, toNextKernelX, nConvoX, toNextConvoX);
            ikernelTable -= nKernelX;
            startConvo += radKernelX;
        }
        // Loop for Y BOTTOM
        ikernelTable += nKernelX;
        for (;iWaveY < nWaveY;iWaveY +=1){
            startKernel = radKernelX;
            endKernelY -= nKernelX;
            ConvolveX (srcFramePtr, destPixPtr, kernel, kernelTable, nKernelX, radKernelX, nWaveX, endKernelY, startKernel, startConvo, ikernelTable, toNextKernelX, nConvoX, toNextConvoX);
            startConvo += radKernelX;
        }
        // if convolving in place copy buffer back on top of src wave  at end of frame
        if (inPlace){
            memcpy ((void*)srcFramePtr, (void*)frameBuffer, fSize * sizeof (TO));
            destPixPtr = frameBuffer; // reset destination pixel to start of buffer
        }
    }
}


/* **************************************** ConvolveFramesThreadParams ***************************************************************
 Structure to pass data to each ConvolveFrames thread or ConvolveSymFrames thread
 Last Modified 2026/10/15 by Jamie Boyd */
typedef struct ConvolveFramesThreadParams{
    int inPutWaveType;          // WaveMetrics code for waveType
    char* inPutDataPtr;         // pointer to start of input wave
    char* outPutDataPtr;        // pointer to start of output wave
    char* frameBufferPtr;       // pointer to frame sized buffer for symConvolveFrames, or when overwriting with ConvolveFrames
    CountInt xSize;             // number of columns in each frame
    CountInt ySize;             // number of rows in each frame
    CountInt startFrame;        // frame to start on for this threrad
    CountInt threadFrames;      // number of frames to do for this thread
    float * kernelDataPtr;      // pointer to start of kernel
    float * kernelTablePtr;     // pointer to start of kernel table
    UInt16 kWidth;              // number of columns in kernel
    UInt16 kHeight;             // number of rows in kernel (ignored by convolveSymFrames)
    UInt8 isFloat;              // waveType of outPut wave. 0 for same type as input wave, 1 for floating point wave
} ConvolveFramesThreadParams, *ConvolveFramesThreadParamsPtr;


/* ************************************ ConvolveFramesThread ******************************************************************
 Each thread to colvolve a range of frames starts with this function
 Last Modified 2014/01/20 by Jamie Boyd */
void *ConvolveFramesThread (void *threadarg){
    ConvolveFramesThreadParams *p = (ConvolveFramesThreadParams*) threadarg;
    CountInt frameSize =p->xSize * p->ySize;
    CountInt startPos = p->startFrame * frameSize; //change start position from frames to data points by multiplying by frame size
	if (p->isFloat){
		switch (p->inPutWaveType) {
            case NT_I8:
                ConvolveT (((char*)p->inPutDataPtr) + startPos, (float*)p->outPutDataPtr + startPos, (float*)p->frameBufferPtr, p->xSize, p->ySize, p->threadFrames, p->kernelDataPtr, p->kernelTablePtr, p->kWidth, p->kHeight);
                break;
            case (NT_I8 | NT_UNSIGNED):
                ConvolveT ((unsigned char*)p->inPutDataPtr + startPos, (float*)p->outPutDataPtr + startPos, (float*)p->frameBufferPtr, p->xSize, p->ySize, p->threadFrames, p->kernelDataPtr, p->kernelTablePtr, p->kWidth, p->kHeight);
                break;
            case NT_I16:
                ConvolveT ((short*)p->inPutDataPtr + startPos, (float*)p->outPutDataPtr + startPos, (float*)p->frameBufferPtr, p->xSize, p->ySize, p->threadFrames, p->kernelDataPtr, p->kernelTablePtr, p->kWidth, p->kHeight);
                break;
            case (NT_I16 | NT_UNSIGNED):
                ConvolveT ((unsigned short*)p->inPutDataPtr + startPos,(float*)p->outPutDataPtr + startPos, (float*)p->frameBufferPtr, p->xSize, p->ySize, p->threadFrames, p->kernelDataPtr, p->kernelTablePtr, p->kWidth, p->kHeight);
                break;
            case NT_I32:
                ConvolveT ((SInt32*)p->inPutDataPtr + startPos,(float*)p->outPutDataPtr + startPos, (float*)p->frameBufferPtr, p->xSize, p->ySize, p->threadFrames, p->kernelDataPtr, p->kernelTablePtr, p->kWidth, p->kHeight);
                break;
            case (NT_I32| NT_UNSIGNED):
                ConvolveT ((UInt32*)p->inPutDataPtr + startPos,(float*)p->outPutDataPtr + startPos, (float*)p->frameBufferPtr, p->xSize, p->ySize, p->threadFrames, p->kernelDataPtr, p->kernelTablePtr, p->kWidth, p->kHeight);
                break;
            case NT_I64:
                ConvolveT((SInt64*)p->inPutDataPtr + startPos, (float*)p->outPutDataPtr + startPos, (float*)p->frameBufferPtr, p->xSize, p->ySize, p->threadFrames, p->kernelDataPtr, p->kernelTablePtr, p->kWidth, p->kHeight);
                break;
            case (NT_I64 | NT_UNSIGNED):
                ConvolveT((UInt64*)p->inPutDataPtr + startPos, (float*)p->outPutDataPtr + startPos, (float*)p->frameBufferPtr, p->xSize, p->ySize, p->threadFrames, p->kernelDataPtr, p->kernelTablePtr, p->kWidth, p->kHeight);
                break;
            case NT_FP32:
                ConvolveT ((float*)p->inPutDataPtr + startPos,(float*)p->outPutDataPtr + startPos, (float*)p->frameBufferPtr, p->xSize, p->ySize, p->threadFrames, p->kernelDataPtr, p->kernelTablePtr, p->kWidth, p->kHeight);
                break;
            case NT_FP64:
                ConvolveT ((double*)p->inPutDataPtr + startPos,(float*)p->outPutDataPtr + startPos, (float*)p->frameBufferPtr, p->xSize, p->ySize, p->threadFrames, p->kernelDataPtr, p->kernelTablePtr, p->kWidth, p->kHeight);
                break;
		}
	}else {
        switch (p->inPutWaveType) {
            case NT_I8:
                ConvolveT ((char*)p->inPutDataPtr + startPos,(char*)p->outPutDataPtr + startPos, (char*)p->frameBufferPtr, p->xSize, p->ySize, p->threadFrames, p->kernelDataPtr, p->kernelTablePtr,p-> kWidth, p->kHeight);
                break;
            case (NT_I8 | NT_UNSIGNED):
                ConvolveT ((unsigned char*)p->inPutDataPtr + startPos,(unsigned char*)p->outPutDataPtr + startPos, (unsigned char*)p->frameBufferPtr, p->xSize, p->ySize, p->threadFrames, p->kernelDataPtr, p->kernelTablePtr, p->kWidth, p->kHeight);
                break;
            case NT_I16:
                ConvolveT ((short*)p->inPutDataPtr + startPos,(short*)p->outPutDataPtr + startPos, (short*)p->frameBufferPtr, p->xSize, p->ySize, p->threadFrames, p->kernelDataPtr, p->kernelTablePtr, p->kWidth, p->kHeight);
                break;
            case (NT_I16 | NT_UNSIGNED):
                ConvolveT ((unsigned short*)p->inPutDataPtr + startPos,(unsigned short*)p->outPutDataPtr + startPos, (unsigned short*)p->frameBufferPtr, p->xSize, p->ySize, p->threadFrames, p->kernelDataPtr, p->kernelTablePtr, p->kWidth, p->kHeight);
                break;
            case NT_I32:
                ConvolveT ((SInt32*)p->inPutDataPtr + startPos,(SInt32*)p->outPutDataPtr + startPos, (SInt32*)p->frameBufferPtr, p->xSize, p->ySize, p->threadFrames, p->kernelDataPtr, p->kernelTablePtr, p->kWidth, p->kHeight);
                break;
            case (NT_I32| NT_UNSIGNED):
                ConvolveT ((UInt32*)p->inPutDataPtr + startPos,(UInt32*)p->outPutDataPtr + startPos, (UInt32*)p->frameBufferPtr, p->xSize, p->ySize, p->threadFrames, p->kernelDataPtr, p->kernelTablePtr, p->kWidth, p->kHeight);
                break;
            case NT_I64:
                ConvolveT((SInt64*)p->inPutDataPtr + startPos, (SInt64*)p->outPutDataPtr + startPos, (SInt64*)p->frameBufferPtr, p->xSize, p->ySize, p->threadFrames, p->kernelDataPtr, p->kernelTablePtr, p->kWidth, p->kHeight);
                break;
            case (NT_I64 | NT_UNSIGNED):
                ConvolveT((UInt64*)p->inPutDataPtr + startPos, (UInt64*)p->outPutDataPtr + startPos, (UInt64*)p->frameBufferPtr, p->xSize, p->ySize, p->threadFrames, p->kernelDataPtr, p->kernelTablePtr, p->kWidth, p->kHeight);
                break;
            case NT_FP32:
                ConvolveT ((float*)p->inPutDataPtr + startPos,(float*)p->outPutDataPtr + startPos, (float*)p->frameBufferPtr, p->xSize, p->ySize, p->threadFrames, p->kernelDataPtr, p->kernelTablePtr, p->kWidth, p->kHeight);
                break;
            case NT_FP64:
                ConvolveT ((double*)p->inPutDataPtr + startPos,(double*)p->outPutDataPtr + startPos, (double*)p->frameBufferPtr, p->xSize, p->ySize, p->threadFrames, p->kernelDataPtr, p->kernelTablePtr, p->kWidth, p->kHeight);
                break;
		}
	}
	return nullptr;
}

/* *********************************************** ConvolveFrames *******************************************************
 ConvolveFrames XOP entry function
 Convolves a 2D or 3D wave with a smaller 2D wave (the kernel), and sends the output to an output wave.
 Treats each plane in a 3D wave as a separate image
 Convolves any type of input wave and outputs to either the same type of wave, or to a 32 bit floating point wave
 Last modified 2025/06/23 by Jamie Boyd
 
ConvolveFramesParams
 double overWrite; // 1 if it is o.k. to overwrite existing waves, 0 to exit with error
 waveHndl kernelH;	// convolution kernel, a 2D wave odd number of pixels high and wide
 double outPutType; // 0 for same type as input wave, non-zero for floating point wave
 Handle outPutPath;	// A handle to a string containing path to output wave we want to make, or empty string to overwrite existing wave
 waveHndl inPutWaveH;//input wave. needs to be 3D wave
 double result; */
extern "C" int ConvolveFrames(ConvolveFramesParamsPtr p){
	int result = 0;                                 // The error returned from various Wavemetrics functions
	waveHndl inPutWaveH, outPutWaveH, kernelH;		// handles to the input and output waves and the kernel
	int inPutWaveType, kernelWaveType;              //  Wavetypes numeric codes for things like 32 bit floating point, 16 bit int, etc
    int inPutDimensions, kernelDimensions;	        // number of dimensions in input and kernel waves
   
    CountInt inPutDimSizes[MAX_DIMENSIONS+1];       // array used to hold the width, height, layers, and chunk sizes of input wave
    CountInt outPutDimSizes[MAX_DIMENSIONS+1];
    CountInt kernelDimSizes[MAX_DIMENSIONS+1];	    // array used to hold the width, height, layers, and chunk sizes of kernel wave
	CountInt zSize;                                 // we use this separate from array for size calculation
    BCInt inPutOffset, outPutOffset, kernelOffset;	// offset in bytes from begnning of handle to a wave to the actual data - size of headers, units, etc.
	DataFolderHandle outPutDFHandle;	// Handle to the datafolder where we will put the output wave
	DFPATH outPutPath;                   // strings to hold data folder paths of input and outPut waves
	WVNAME outPutWaveName;           // C strings to hold names of input and output waves
	UInt8 overWriteOK = (UInt8)(p->overWrite);	    // 0 to not overwrite output wave if it already exists, 1 to overwrite old waves
	UInt8 isFloat = (UInt8)(p-> outPutType);        // 0 to use input type for output, non-zero to use make bit floating point output
	UInt8 isOverWritingInput =0;                    // will be set to 1 if we are overwriting the input wave
    int outPutWaveType;
    int outPutDimensions;
    int requestedOutPutType;                           // will be inputWaveType ot NT_32
    float *kernelTablePtr;                          // kernel table (calculated weighting for truncation of convolution at edges)
	char *inPutDataStartPtr, *outPutDataStartPtr, *kernelDataStartPtr;
    // for threads
    UInt8 iThread, nThreads;
    CountInt threadFrames;
    char **bufferPtr = nullptr;  // pointer to array of pointers
    ConvolveFramesThreadParamsPtr paramArrayPtr = nullptr;
#ifdef __GNUC__
    pthread_t* threadsPtr = nullptr;   // pointer to an array of pThread_t structures we will make
#endif
#ifdef  _WINDOWS_
    HANDLE *threadsPtr = nullptr;      // pointer to an array of Windows handles, which will point to thread structures we make
    DWORD *threadIDsPtr = nullptr;     // pointer to an array of DWORD (32-bit unsigned integer) that Windows uses for thread IDs
#endif
    try{
        // Get handles to input wave and kernel. Make sure both waves exist.
        inPutWaveH = p->inPutWaveH;
        kernelH = p ->kernelH;
        if ((inPutWaveH == nullptr)||(kernelH == nullptr)) throw result = NON_EXISTENT_WAVE;
        // Get wave data type
        inPutWaveType = WaveType(inPutWaveH);
        if (inPutWaveType==TEXT_WAVE_TYPE) throw result = NOTEXTWAVES;
        // Get number of used dimensions in waves.
        if (MDGetWaveDimensions(inPutWaveH, &inPutDimensions, inPutDimSizes)) throw result = WAVEERROR_NOS;
        // Check that inputwave is 2D or 3D
        if ((inPutDimensions == 1) || (inPutDimensions == 4)) throw result = INPUTNEEDS_2D3D_WAVE;
        // if z size in array is 0 (for a 2D wave), set zSize to 1 for purposes of calculating size
        if (inPutDimSizes [2] == 0)
            zSize = 1;
        else
            zSize=inPutDimSizes [2];
        //	Check that kernel is 2D, and of odd size in both x and y
        if (MDGetWaveDimensions(kernelH, &kernelDimensions, kernelDimSizes)) throw result = WAVEERROR_NOS;
        if ((((kernelDimensions != 2) || ((kernelDimSizes[0] % 2) == 0)) || ((kernelDimSizes[1] % 2) == 0))) throw result = BADKERNEL;
        // make sure kernel is 32 bit float - change if necessary
        kernelWaveType = WaveType (kernelH);
        if (kernelWaveType !=  (NT_FP32)){
            if (MDChangeWave (kernelH, NT_FP32, kernelDimSizes))throw result = WAVEERROR_NOS;
        }
        // If outPutPath is empty string, we are being asked to overwrite the input wave, which may or may not be OK
        if (WMGetHandleSize (p->outPutPath) == 0){
            if (!(overWriteOK)) throw result = OVERWRITEALERT;
            if ((isFloat) && (inPutWaveType != NT_FP32)){
                if (MDChangeWave (inPutWaveH, NT_FP32, inPutDimSizes)) throw result = WAVEERROR_NOS;
                inPutWaveType = NT_FP32;
            }
            outPutWaveH = inPutWaveH;
            isOverWritingInput = 1;
        }else{ // we have a passed in string we need to parse
            // Parse passed in outPut path into dataFolder path and waveName
            ParseWavePath (p->outPutPath, outPutPath, outPutWaveName);
            //No liberal wave names for output wave, maybe we don't really need to enforce non-liberal names?
            CleanupName (0, outPutWaveName, MAX_OBJ_NAME);
            // Get a handle to the output datafolder
            if (GetNamedDataFolder (NULL, outPutPath, &outPutDFHandle))throw result = WAVEERROR_NOS;
            // try to get a handle to the output wave specified in the string
            outPutWaveH = FetchWaveFromDataFolder(outPutDFHandle, outPutWaveName);
            if (outPutWaveH == nullptr){        // output wave specified in wave name does not exist, so make it
                if (isFloat){
                    if (MDMakeWave (&outPutWaveH, outPutWaveName, outPutDFHandle, inPutDimSizes, NT_FP32, overWriteOK)) throw result = WAVEERROR_NOS;
                }else{
                    if (MDMakeWave (&outPutWaveH, outPutWaveName, outPutDFHandle, inPutDimSizes, inPutWaveType, overWriteOK)) throw result = WAVEERROR_NOS;
                }
            }else{ // wave specified in the string already exists and outPutWaveH now points to it
                if (!(overWriteOK)) throw result = OVERWRITEALERT;  // because output wave already exists
                if (outPutWaveH == inPutWaveH){ // overwriting input wave
                    if ((isFloat) && (inPutWaveType != NT_FP32)){
                        if (MDChangeWave (inPutWaveH, NT_FP32, inPutDimSizes)) throw result = WAVEERROR_NOS;
                        inPutWaveType = NT_FP32;
                    }
                    isOverWritingInput = 1;
                }else{  // outPut wave is an existing wave, and is not the input wave, redimension and change type if needed
                    // get wave data type of output wave and check that we don't have a text wave
                    outPutWaveType = WaveType(outPutWaveH);
                    if (outPutWaveType==TEXT_WAVE_TYPE) throw result = NOTEXTWAVES;
                    //Get number of dimensions and dimension sizes in output wave
                    if (MDGetWaveDimensions(outPutWaveH, &outPutDimensions, outPutDimSizes))throw result = WAVEERROR_NOS;
                    // redimension output wave and change wave type if necessary. it needs to be 2D,
                    // same X and Y size as input, and float data type
                    if (isFloat)
                        requestedOutPutType = NT_FP32;
                    else
                        requestedOutPutType = inPutWaveType;
                    if ((outPutWaveType |= requestedOutPutType) || (outPutDimSizes [0] |= inPutDimSizes [0]) || (outPutDimSizes [1] |= inPutDimSizes [1]) || (outPutDimSizes [2] |= inPutDimSizes [2])){
                        if (MDChangeWave(outPutWaveH, requestedOutPutType, inPutDimSizes)) throw result = WAVEERROR_NOS;
                    }
                }
            }
        }
        //Get data offsets for the 3 waves (2 waves, if overwriting)
        if (MDAccessNumericWaveData(inPutWaveH, kMDWaveAccessMode0, &inPutOffset)) throw result = WAVEERROR_NOS;
		if (isOverWritingInput){
			outPutOffset = inPutOffset;
		}else{
			if (MDAccessNumericWaveData(outPutWaveH, kMDWaveAccessMode0, &outPutOffset)) throw result = WAVEERROR_NOS;
		}
        if (MDAccessNumericWaveData(kernelH, kMDWaveAccessMode0, &kernelOffset)) throw result = WAVEERROR_NOS;
        inPutDataStartPtr = (char*)(*inPutWaveH) + inPutOffset;
		outPutDataStartPtr =  (char*)(*outPutWaveH) + outPutOffset;
        kernelDataStartPtr = (char*)(*kernelH) + kernelOffset;
        //make kernel table
		kernelTablePtr = ConvolveMakeKernelTable ((float*)kernelDataStartPtr, (int)kernelDimSizes[0], (int) kernelDimSizes[1]);
		if (kernelTablePtr == nullptr) throw result = NOMEM;
        // multiprocessor init
        nThreads = gNumProcessors;
        if (zSize < nThreads) nThreads = zSize;
        threadFrames = zSize/nThreads;
        // make frame sized buffers if overwriting input
        if (isOverWritingInput){ // input = output wave, so need to  make a frame sized buffer
            bufferPtr = (char**)WMNewPtr(sizeof(char*) * nThreads);
            if (bufferPtr == nullptr) throw result = MEMFAIL;
            for (iThread = 0 ; iThread < nThreads; iThread +=1){
                switch (inPutWaveType) {
                    case NT_I64 | NT_UNSIGNED:
                    case NT_I64:
                    case NT_FP64:
                        bufferPtr [iThread] = (char*)WMNewPtr (inPutDimSizes[ROWS] * inPutDimSizes[COLUMNS]  * 8);
                        if (bufferPtr [iThread] == nullptr) throw result = MEMFAIL;
                        break;
                    case NT_I32 | NT_UNSIGNED:
                    case NT_I32:
                    case NT_FP32:
                        bufferPtr[iThread] = (char*)WMNewPtr (inPutDimSizes[ROWS] * inPutDimSizes[COLUMNS]  * 4);
                        if (bufferPtr [iThread] == nullptr) throw result = MEMFAIL;
                        break;
                    case NT_I16 | NT_UNSIGNED:
                    case NT_I16:
                        bufferPtr [iThread] = (char*)WMNewPtr (inPutDimSizes[ROWS] * inPutDimSizes[COLUMNS]   * 2);
                        if (bufferPtr [iThread] == nullptr) throw result = MEMFAIL;
                        break;
                    case NT_I8 | NT_UNSIGNED:
                    case NT_I8:
                        bufferPtr [iThread] = (char*)WMNewPtr (inPutDimSizes[ROWS] * inPutDimSizes[COLUMNS]   * 1);
                        if (bufferPtr [iThread] == nullptr) throw result = MEMFAIL;
                        break;
                    default:
                        throw result = NUMTYPE;
                        break;
                }
            }
        }
        // make array with data for threads
        paramArrayPtr = (ConvolveFramesThreadParamsPtr)WMNewPtr (nThreads * sizeof(ConvolveFramesThreadParams));
        if (paramArrayPtr == nullptr) throw result = NOMEM;
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
 
    }catch (int (result)) { // catch before starting threads
        if (bufferPtr != nullptr){
            for (iThread = 0 ; iThread < nThreads; iThread +=1){
                if (bufferPtr [iThread] == nullptr)
                    break;
                else
                    WMDisposePtr ((Ptr)bufferPtr [iThread]);
            }
            WMDisposePtr ((Ptr)bufferPtr);
        }
        if (paramArrayPtr != nullptr) WMDisposePtr ((Ptr)paramArrayPtr);
        if (threadsPtr != nullptr) WMDisposePtr ((Ptr)threadsPtr);
#ifdef _WINDOWS_
        if (threadIDsPtr != nullptr) WMDisposePtr((Ptr)threadIDsPtr);
#endif
        if (kernelTablePtr !=nullptr) WMDisposePtr ((Ptr)kernelTablePtr);
        WMDisposeHandle (p->outPutPath);    // free input string for output path
        // set result
        p -> result = (double)(result - FIRST_XOP_ERR);
        #ifdef NO_IGOR_ERR
            return (0);
        #else
            return (result);
        #endif
    } // end of try-catch block
    // fill array of parameter structures
    for (iThread = 0; iThread < nThreads; iThread++){
        paramArrayPtr[iThread].inPutWaveType = inPutWaveType;
        paramArrayPtr[iThread].inPutDataPtr = inPutDataStartPtr;
        paramArrayPtr[iThread].outPutDataPtr = outPutDataStartPtr;
        paramArrayPtr[iThread].xSize = inPutDimSizes [0];
        paramArrayPtr[iThread].ySize = inPutDimSizes [1];
        paramArrayPtr[iThread].startFrame = iThread * threadFrames;
        paramArrayPtr[iThread].threadFrames = threadFrames;
        if (isOverWritingInput){
            paramArrayPtr[iThread].frameBufferPtr = bufferPtr[iThread];
        }
        paramArrayPtr[iThread].kernelDataPtr=(float*)kernelDataStartPtr;
        paramArrayPtr[iThread].kernelTablePtr=kernelTablePtr;
        paramArrayPtr[iThread].kWidth = kernelDimSizes [0];
        paramArrayPtr[iThread].kHeight= kernelDimSizes [1];
        paramArrayPtr[iThread].isFloat = isFloat; // 0 for same type as input wave, non-zero for floating point wave
    }
    paramArrayPtr[nThreads -1].threadFrames += zSize % threadFrames;    // last thread gets any left-over frames
#ifdef __GNUC__
    // create the threads
    for (iThread = 0; iThread < nThreads; iThread++){
        pthread_create(&threadsPtr[iThread], NULL, ConvolveFramesThread, (void*)&paramArrayPtr[iThread]);
    }
    // wait for threads to finish
    for (iThread = 0; iThread < nThreads; iThread++) {
        pthread_join(threadsPtr[iThread], NULL);
    }
#endif
#ifdef _WINDOWS_
    // create the threads with Default security attributes, Default stack size, KalmanThread function, a pointer to a parameter structure, Default creation flags, and a pointer to a threadId
    for (iThread = 0; iThread < nThreads; iThread++){
        threadsPtr[iThread] = CreateThread( NULL, 0, (LPTHREAD_START_ROUTINE)ConvolveFramesThread, (LPVOID)&paramArrayPtr[iThread], 0, &threadIDsPtr[iThread]);
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

    // free frameBuffer, if made
    if (isOverWritingInput){
        for (iThread = 0 ; iThread < nThreads; iThread +=1){
            WMDisposePtr ((Ptr)bufferPtr [iThread]);
        }
        WMDisposePtr ((Ptr)bufferPtr);
    }
    //free kernel table memory
    WMDisposePtr ((Ptr)kernelTablePtr);
    WMDisposeHandle (p->outPutPath);
    // Inform Igor that we have changed the output wave.
    WaveHandleModified(outPutWaveH);
    p -> result = (0);
    return (0);
}

/* -------------------------------- convolution with a symetrical 2D kernel ------------------------------------------
 A symetrical 2D kernel applied as the same 1D kernel done first in X and then in Y. Each frame is done as 2D plane
 ------------------------------------------------------------------------------------------------------------------- */
 

/* ************************************** SymConvolveMakeKernelTable ************************************************
 Makes and returns a table of weights for a 1D kernel, first Checking that the kernel is, indeed, symetrical
 returns a null pointer if kernel is not symetrical
 Last Modified 2014/01/21 by Jamie Boyd */
float * SymConvolveMakeKernelTable (float * kernel, UInt16 kWidth){
    // kernel has to be odd
    if((kWidth % 2) == 0) return NIL;
    // check symetry
    // Get sum of entire kernel for making kernel table while we are at it
    UInt16 kRadW = (kWidth-1)/2; // radius of kernel, not counting central point
    UInt16 iKT, iConvo;
    float kSum=0;
    for (iKT =0, iConvo = kWidth -1; iKT < kRadW; iKT +=1, iConvo -=1){
        if (kernel [iKT] != kernel [iConvo]) return NIL;
        kSum += kernel [iKT] + kernel [iConvo];
    }
    // If we get to here, kernel is symetrical. Sum of kernel lacks central point
    kSum += kernel [kRadW];
    // make kernel table. It's symetrical, so we only need to make size kRadW + 1
    float *kernelTable =(float*) WMNewPtr ((kRadW + 1) * sizeof (float));
    if (kernelTable == NIL) return NIL;
    // last point = kernel sum
    kernelTable [kRadW] = kSum;
    // each earlier point removes a value from left side
    iKT = kRadW-1;
    for (iConvo = 0; iKT >0; iKT -=1, iConvo +=1)
        kernelTable [iKT] = kernelTable [iKT + 1] - kernel [iConvo];
    kernelTable [iKT] = kernelTable [iKT + 1] - kernel [iConvo];
    return kernelTable;
}


/* ************************************** SymConvolveT **********************************************************************
 template function for convolving one wave with a 1D symetrical kernel (e.g., Guassian) and putting results in an output wave.
 Input wave can be 2 or 3D, but each plane is done as a separate 2D image.
 Last Modified 2014/01/21 by Jamie Boyd */
template <typename TI, typename TO> void SymConvolveT (TI* srcWave, TO* destWave, TO* buffer, CountInt nWaveX, CountInt nWaveY, CountInt nWaveZ, float* kernel, float * kernelTable, UInt16 kWidth){
    CountInt frameSize = nWaveX * nWaveY;
    //variables for iterating through waves
    UInt16 kRad = (kWidth-1)/2; // radius of kernel, not including the central point
    //CountInt nWave = frameSize * nWaveZ;
    //Position of pixel in x,y,z
    CountInt iWaveY,  iWaveX;
    //CountInt nEdgeY = kRad * nWaveX;
    // iterating through kernel and convo positions, linear
    CountInt iKernel, iConvo, iKernelTable=0, convoStart, kernelStart, kernelEnd;
    // temporary value to calculate value for each pixel
    double outVal;
    CountInt iOut; // position of output value, in buffer, then in output wave
    TI* srcFramePtr = srcWave;
    TI* srcEndPtr = srcWave + frameSize * nWaveZ;
    TO* destFramePtr = destWave;
    for (; srcFramePtr < srcEndPtr; srcFramePtr += frameSize, destFramePtr += frameSize){
        // Do Y filtering first, sending output to frame buffer. each Y starts at startFrameP + iWaveX
        for (iWaveX =0; iWaveX < nWaveX; iWaveX +=1){
            // TOP of each row.
            iOut = iWaveX;
            convoStart = iWaveX;
            kernelStart = kRad;
            kernelEnd = kWidth;
            for (iWaveY=0; iWaveY < kRad; iWaveY += 1){
                // convolution in Y
                outVal =0;
                for (iKernel =kernelStart, iConvo =convoStart; iKernel < kernelEnd; iKernel +=1, iConvo += nWaveX)
                    outVal += srcFramePtr [iConvo] * kernel [iKernel];
                buffer [iOut] = outVal / kernelTable [iKernelTable];
                kernelStart -=1;
                iKernelTable +=1;
                iOut += nWaveX;
            }
            // Middle of each row
            for (; iWaveY < (nWaveY - kRad); iWaveY += 1){
                outVal =0;
                // convolution in Y
                for (iKernel =kernelStart, iConvo =convoStart; iKernel < kernelEnd; iKernel +=1, iConvo += nWaveX)
                    outVal += srcFramePtr [iConvo] * kernel [iKernel];
                buffer [iOut] = outVal / kernelTable [iKernelTable];
                convoStart += nWaveX;
                iOut += nWaveX;
            }
            // bottom of each row
            for (; iWaveY < nWaveY; iWaveY += 1){
                kernelEnd -=1;
                iKernelTable -=1;
                // convolution in Y
                outVal =0;
                for (iKernel =kernelStart, iConvo =convoStart; iKernel < kernelEnd; iKernel +=1, iConvo += nWaveX)
                    outVal += srcFramePtr [iConvo] * kernel [iKernel];
                buffer [iOut] = outVal / kernelTable [iKernelTable];
                convoStart += nWaveX;
                iOut += nWaveX;
            }
        }
        // Do X filtering second, sending output to output wave.
        iOut =0;
        convoStart = 0;
        for (iWaveY=0; iWaveY < nWaveY; iWaveY +=1){
            kernelStart = kRad;
            kernelEnd = kWidth;
            // left edge
            for (iWaveX =0; iWaveX < kRad; iWaveX +=1){
                // convolution in X
                outVal =0;
                for (iKernel =kernelStart, iConvo =convoStart; iKernel < kernelEnd; iKernel +=1, iConvo += 1)
                    outVal += buffer [iConvo] * kernel [iKernel];
                destFramePtr [iOut] = outVal / kernelTable [iKernelTable];
                kernelStart -=1;
                iKernelTable +=1;
                iOut += 1;
            }
            // Middle
            for (; iWaveX < (nWaveX - kRad); iWaveX+= 1){
                outVal =0;
                // convolution in X
                for (iKernel =kernelStart, iConvo =convoStart; iKernel < kernelEnd; iKernel +=1, iConvo += 1)
                    outVal += buffer [iConvo] * kernel [iKernel];
                destFramePtr [iOut] = outVal / kernelTable [iKernelTable];
                convoStart += 1;
                iOut += 1;
            }
            // right edge
            for (; iWaveX< nWaveX; iWaveX += 1){
                kernelEnd -=1;
                iKernelTable -=1;
                // convolution in X
                outVal =0;
                for (iKernel =kernelStart, iConvo =convoStart; iKernel < kernelEnd; iKernel +=1, iConvo += 1)
                    outVal += buffer [iConvo] * kernel [iKernel];
                destFramePtr [iOut] = outVal/ kernelTable [iKernelTable];
                convoStart += 1;
                iOut += 1;
            }
            convoStart += kRad; // get to next row for x filtering
        }
    }
}

/* ******************************************* SymConvolveFramesThread ************************************************************
 Each thread to symetrically convolve a range of frames starts with this function
 Last Modified 2025/06/23 by Jamie Boyd */
void* SymConvolveFramesThread (void* threadarg){
    struct ConvolveFramesThreadParams* p;
    p = (struct ConvolveFramesThreadParams*) threadarg;
    CountInt frameSize = p->xSize * p->ySize;
    CountInt startPos = p->startFrame * frameSize;
    //CountInt tFrames = p->zSize/p->tN; //frames per thread = number of frames / number of threads, truncated to an integer
    //CountInt startPos = p->ti * tFrames; // whaich frame to start this thread on depends on thread number * frames per thread. ti is 0 based
    //if (p->ti == p->tN - 1) tFrames +=  (p->zSize % p->tN); // the last thread gets any left-over frames
    //startPos *= (p->xSize * p->ySize); //change start position from frames to data points by multiplying by frame size
    //CountInt bufferOffset =p->ti * p->xSize * p->ySize;
    
    if (p->isFloat){
        switch (p->inPutWaveType) {
            case NT_I8:
                SymConvolveT ((char*)p->inPutDataPtr + startPos,(float*)p->outPutDataPtr + startPos, (float*)p->frameBufferPtr, p->xSize, p->ySize, p->threadFrames, p->kernelDataPtr, p->kernelTablePtr, p-> kWidth);
                break;
            case (NT_I8 | NT_UNSIGNED):
                SymConvolveT ((unsigned char*)p->inPutDataPtr + startPos,(float*)p->outPutDataPtr + startPos, (float*)p->frameBufferPtr, p->xSize, p->ySize, p->threadFrames, p->kernelDataPtr, p->kernelTablePtr, p-> kWidth);
                break;
            case NT_I16:
                SymConvolveT ((short*)p->inPutDataPtr + startPos,(float*)p->outPutDataPtr + startPos, (float*)p->frameBufferPtr, p->xSize, p->ySize, p->threadFrames, p->kernelDataPtr, p->kernelTablePtr, p-> kWidth);
                break;
            case (NT_I16 | NT_UNSIGNED):
                SymConvolveT ((unsigned short*)p->inPutDataPtr + startPos,(float*)p->outPutDataPtr + startPos, (float*)p->frameBufferPtr, p->xSize, p->ySize, p->threadFrames, p->kernelDataPtr, p->kernelTablePtr, p-> kWidth);
                break;
            case NT_I32:
                SymConvolveT ((SInt32*)p->inPutDataPtr + startPos,(float*)p->outPutDataPtr + startPos, (float*)p->frameBufferPtr, p->xSize, p->ySize, p->threadFrames, p->kernelDataPtr, p->kernelTablePtr, p-> kWidth);
                break;
            case (NT_I32| NT_UNSIGNED):
                SymConvolveT ((UInt32*)p->inPutDataPtr + startPos,(float*)p->outPutDataPtr + startPos, (float*)p->frameBufferPtr, p->xSize, p->ySize, p->threadFrames, p->kernelDataPtr, p->kernelTablePtr, p-> kWidth);
                break;
            case NT_I64:
                SymConvolveT((SInt64*)p->inPutDataPtr + startPos, (float*)p->outPutDataPtr + startPos, (float*)p->frameBufferPtr, p->xSize, p->ySize, p->threadFrames, p->kernelDataPtr, p->kernelTablePtr, p->kWidth);
                break;
            case (NT_I64 | NT_UNSIGNED):
                SymConvolveT((UInt64*)p->inPutDataPtr + startPos, (float*)p->outPutDataPtr + startPos, (float*)p->frameBufferPtr, p->xSize, p->ySize, p->threadFrames, p->kernelDataPtr, p->kernelTablePtr, p->kWidth);
                break;
            case NT_FP32:
                SymConvolveT ((float*)p->inPutDataPtr + startPos,(float*)p->outPutDataPtr + startPos, (float*)p->frameBufferPtr, p->xSize, p->ySize, p->threadFrames, p->kernelDataPtr, p->kernelTablePtr, p-> kWidth);
                break;
            case NT_FP64:
                SymConvolveT ((double*)p->inPutDataPtr + startPos,(float*)p->outPutDataPtr + startPos, (float*)p->frameBufferPtr, p->xSize, p->ySize, p->threadFrames, p->kernelDataPtr, p->kernelTablePtr, p-> kWidth);
                break;
        }
    }else {
        switch (p->inPutWaveType) {
            case NT_I8:
                SymConvolveT ((char*)p->inPutDataPtr + startPos,(char*)p->outPutDataPtr + startPos, (char*)p->frameBufferPtr, p->xSize, p->ySize, p->threadFrames, p->kernelDataPtr, p->kernelTablePtr, p-> kWidth);
                break;
            case (NT_I8 | NT_UNSIGNED):
                SymConvolveT ((unsigned char*)p->inPutDataPtr + startPos,(unsigned char*)p->outPutDataPtr + startPos, (unsigned char*)p->frameBufferPtr, p->xSize, p->ySize, p->threadFrames, p->kernelDataPtr, p->kernelTablePtr, p-> kWidth);
                break;
            case NT_I16:
                SymConvolveT ((short*)p->inPutDataPtr + startPos,(short*)p->outPutDataPtr + startPos,  (short*)p->frameBufferPtr, p->xSize, p->ySize, p->threadFrames, p->kernelDataPtr, p->kernelTablePtr, p-> kWidth);
                break;
            case (NT_I16 | NT_UNSIGNED):
                SymConvolveT ((unsigned short*)p->inPutDataPtr + startPos,(unsigned short*)p->outPutDataPtr + startPos, (unsigned short*)p->frameBufferPtr, p->xSize, p->ySize, p->threadFrames, p->kernelDataPtr, p->kernelTablePtr, p-> kWidth);
                break;
            case NT_I32:
                SymConvolveT ((SInt32*)p->inPutDataPtr + startPos,(SInt32*)p->outPutDataPtr + startPos, (SInt32*)p->frameBufferPtr, p->xSize, p->ySize, p->threadFrames, p->kernelDataPtr, p->kernelTablePtr, p-> kWidth);
                break;
            case (NT_I32| NT_UNSIGNED):
                SymConvolveT ((UInt32*)p->inPutDataPtr + startPos,(UInt32*)p->outPutDataPtr + startPos, (UInt32*)p->frameBufferPtr, p->xSize, p->ySize, p->threadFrames, p->kernelDataPtr, p->kernelTablePtr, p-> kWidth);
                break;
            case NT_I64:
                SymConvolveT((SInt64*)p->inPutDataPtr + startPos, (SInt64*)p->outPutDataPtr + startPos, (SInt64*)p->frameBufferPtr, p->xSize, p->ySize, p->threadFrames, p->kernelDataPtr, p->kernelTablePtr, p->kWidth);
                break;
            case (NT_I64 | NT_UNSIGNED):
                SymConvolveT((UInt64*)p->inPutDataPtr + startPos, (UInt64*)p->outPutDataPtr + startPos, (UInt64*)p->frameBufferPtr, p->xSize, p->ySize, p->threadFrames, p->kernelDataPtr, p->kernelTablePtr, p->kWidth);
                break;
            case NT_FP32:
                SymConvolveT ((float*)p->inPutDataPtr + startPos,(float*)p->outPutDataPtr + startPos, (float*)p->frameBufferPtr, p->xSize, p->ySize, p->threadFrames, p->kernelDataPtr, p->kernelTablePtr, p-> kWidth);
                break;
            case NT_FP64:
                SymConvolveT ((double*)p->inPutDataPtr + startPos,(double*)p->outPutDataPtr + startPos, (double*)p->frameBufferPtr, p->xSize, p->ySize, p->threadFrames, p->kernelDataPtr, p->kernelTablePtr, p-> kWidth);
                break;
        }
    }
    return nullptr;
}

/* SymConvolveFrames XOP entry function
 Convolves a 2D or 3D wave with a 1D symetrical kernel in X and Y, and sends the output to an output wave.
 Treats each plane in a 3D wave as a separate image
 Convolves any type of input wave and outputs to either the same type of wave,  or to a 32 bit floating point wave
 Last modified 2025/10/20 by Jamie Boyd
 
 typedef struct ConvolveFramesParams{
 double overWrite; // 1 if it is o.k. to overwrite existing waves, 0 to exit with error if overwriting will occur
 waveHndl kernelH; // convolution kernel, a 2D wave odd number of pixels high and wide, or 1D odd number wave for Sym
 double outPutType; // 0 for same type as input wave, non-zero for floating point wave
 Handle outPutPath;	// A handle to a string containing path to output wave we want to make, or empty string to overwrite existing wave
 waveHndl inPutWaveH; //input wave. needs to be 2D or 3D wave
 double result; */
extern "C" int SymConvolveFrames(ConvolveFramesParamsPtr p) {
    int result = 0;	                                // The error returned from various Wavemetrics functions
    waveHndl inPutWaveH;		                    // handle to the input wave
    int inPutWaveType;                              // Wavetypes numeric codes for things like 32 bit floating point, 16 bit int, etc
    int inPutDimensions;	                        // number of dimensions in input wave
    CountInt inPutDimSizes[MAX_DIMENSIONS+1];	    // an array used to hold the width, height, layers, and chunk sizes
    BCInt inPutOffset;                              // offset to data in input wave, after headers, etc.
    char *inPutDataStartPtr;                        // pointer to start of data in input wave
    waveHndl outPutWaveH;                           //handle to output wave, which we wil create
    int outPutWaveType;                             // wavetype for output wave, needed if wave already exists
    int outPutDimensions;                           // number of dimensions in output wave, needed if wave already exists
    BCInt outPutOffset;                             // offset to data in output wave
    char *outPutDataStartPtr;                       // pointer to start of data in output wave
    CountInt outPutDimSizes[MAX_DIMENSIONS+1];      // an array used to hold the width, height, layers, and chunk sizes
    DataFolderHandle outPutDFHandle;                // Handle to the datafolder where we will put the output wave
    DFPATH outPutPath;                              // string to hold data folder paths of outPut wave
    WVNAME outPutWaveName;                          // C string to hold name of output wave
    waveHndl kernelH;                               // handle to wave used as convolution kernel
    int kernelType;                                 // wave type - we will make it float
    int kernelDimensions;                           // number of deimensions in kernal wave - needs to be 2
    CountInt kernelDimSizes[MAX_DIMENSIONS+1];      // an array used to hold the width, height, layers, and chunk sizes
    BCInt kernelOffset;                             //offset in bytes from begnning of handle to a wave to the actual data
    char *kernelDataStartPtr;                       // pointer to start of data in kernel wave
    UInt8 overWriteOK = (UInt8)(p->overWrite);        // 0 to not overwrite output wave if it already exists, 1 to overwrite old waves
    UInt8 isFloat = (UInt8)(p-> outPutType);        // 0 to use input type for output, non-zero to use make bit floating point output
    UInt8 isOverWritingInput = 0;                    // will be set to 1 if we are overwriting the input wave
   
    int requestedOutPutType;                         // will be floating point or same as input wave
    CountInt zSize;
    char **bufferPtr = nullptr;                     // pointer to array of pointers for buffers
    float *kernelTable = nullptr;                   // pointer to kernel Table we make
    UInt8 iThread, nThreads;                        // for iterating through thread data
    CountInt threadFrames;                          // number of frames each thread processes
    ConvolveFramesThreadParamsPtr paramArrayPtr;
#ifdef __GNUC__
    pthread_t* threadsPtr = nullptr;   // pointer to an array of pThread_t structures we will make
#endif
#ifdef  _WINDOWS_
    HANDLE *threadsPtr = nullptr;      // pointer to an array of Windows handles, which will point to thread structures we make
    DWORD *threadIDsPtr = nullptr;     // pointer to an array of DWORD (32-bit unsigned integer) that Windows uses for thread IDs
#endif

    try{
        // Get handles to input wave and kernel.
        inPutWaveH = p->inPutWaveH;
        kernelH = p->kernelH;
        if ((kernelH == nullptr) || (inPutWaveH == nullptr))throw result = NON_EXISTENT_WAVE;
        // Get wave data type
        inPutWaveType = WaveType(inPutWaveH);
        if (inPutWaveType==TEXT_WAVE_TYPE) throw result = NOTEXTWAVES;
        // Get number of used dimensions in waves.
        if (MDGetWaveDimensions(inPutWaveH, &inPutDimensions, inPutDimSizes)) throw result = WAVEERROR_NOS;
        // Check that inputwave is 2D or 3D
        if ((inPutDimensions == 1) || (inPutDimensions == 4)) throw result = INPUTNEEDS_2D3D_WAVE;
        // if layers is 0, make it 1 for purposes of calculating size
        if (inPutDimSizes [LAYERS] == 0)
            zSize = 1;
        else
            zSize=inPutDimSizes [LAYERS];        //	make sure kernel exists and is 1D and of odd dimensions
        // Get number of used dimensions in kernel. Make sure it is 1D and of odd length
        if (MDGetWaveDimensions(kernelH, &kernelDimensions, kernelDimSizes)) throw result = WAVEERROR_NOS;
        if ((kernelDimensions > 1) || (kernelDimSizes [0] % 2 == 0))  throw result = BADSYMKERNEL;
        // Get kernel wave data type, and make other numeric types float
        kernelType = WaveType(kernelH);
        if (kernelType==TEXT_WAVE_TYPE) throw result = NOTEXTWAVES;
        if (kernelType != NT_FP32){
            if  (MDChangeWave(kernelH, NT_FP32, kernelDimSizes)) throw result = WAVEERROR_NOS;
        }
        // make output wave
        // If outPutPath is empty string, we are being asked to overwrite the input wave, which may or may not be OK
        if (WMGetHandleSize (p->outPutPath) == 0){
            if (!(overWriteOK)) throw result = OVERWRITEALERT;
            if ((isFloat) && (inPutWaveType != NT_FP32)){
                if (MDChangeWave (inPutWaveH, NT_FP32, inPutDimSizes)) throw result = WAVEERROR_NOS;
            }
            outPutWaveH = inPutWaveH;
            isOverWritingInput = 1;
        }else{ // we have a passed in string we need to parse
            // Parse passed in outPut path into dataFolder path and waveName
            ParseWavePath (p->outPutPath, outPutPath, outPutWaveName);
            //No liberal wave names for output wave, maybe we don't really need to enforce non-liberal names?
            CleanupName (0, outPutWaveName, MAX_OBJ_NAME);
            // Get a handle to the output datafolder
            if (GetNamedDataFolder (NULL, outPutPath, &outPutDFHandle))throw result = WAVEERROR_NOS;
            // try to get a handle to the output wave specified in the string
            outPutWaveH = FetchWaveFromDataFolder(outPutDFHandle, outPutWaveName);
            if (outPutWaveH == nullptr){        // output wave specified in wave name does not exist, so make it
                if (isFloat){
                    if (MDMakeWave (&outPutWaveH, outPutWaveName, outPutDFHandle, inPutDimSizes, NT_FP32, overWriteOK)) throw result = WAVEERROR_NOS;
                }else{
                    if (MDMakeWave (&outPutWaveH, outPutWaveName, outPutDFHandle, inPutDimSizes, inPutWaveType, overWriteOK)) throw result = WAVEERROR_NOS;
                }
            }else{ // wave specified in the string already exists and outPutWaveH now points to it
                if (!(overWriteOK)) throw result = OVERWRITEALERT;
                if (outPutWaveH == inPutWaveH){ // overwriting input wave
                    if ((isFloat) && (inPutWaveType != NT_FP32)){
                        if (MDChangeWave (inPutWaveH, NT_FP32, inPutDimSizes)) throw result = WAVEERROR_NOS;
                        inPutWaveType = NT_FP32;
                    }
                    isOverWritingInput = 1;
                }else{  // outPut wave is an existing wave, and is not the input wave, redimension and change type if needed
                    // get wave data type of output wave and check that we don't have a text wave
                    outPutWaveType = WaveType(outPutWaveH);
                    if (outPutWaveType==TEXT_WAVE_TYPE) throw result = NOTEXTWAVES;
                    //Get number of dimensions and dimension sizes in output wave
                    if (MDGetWaveDimensions(outPutWaveH, &outPutDimensions, outPutDimSizes))throw result = WAVEERROR_NOS;
                    if (isFloat)
                        requestedOutPutType = NT_FP32;
                    else
                        requestedOutPutType = inPutWaveType;
                    // redimension output wave and change wave type if necessary. it needs to be 2D, same X and Y size as input, and float data type
                    if ((outPutWaveType |= requestedOutPutType) || (outPutDimSizes [0] |= inPutDimSizes [0]) || (outPutDimSizes [1] |= inPutDimSizes [1]) || (outPutDimSizes [2] |= inPutDimSizes [2])){
                        if (MDChangeWave(outPutWaveH, requestedOutPutType, inPutDimSizes)) throw result = WAVEERROR_NOS;
                    }
                }
            }
        }
        //Get data offsets for the 3 waves (2 waves, if overwriting)
        if (MDAccessNumericWaveData(inPutWaveH, kMDWaveAccessMode0, &inPutOffset)) throw result = WAVEERROR_NOS;
        if (isOverWritingInput){
            outPutOffset = inPutOffset;
        }else{
            if (MDAccessNumericWaveData(outPutWaveH, kMDWaveAccessMode0, &outPutOffset)) throw result = WAVEERROR_NOS;
        }
        inPutDataStartPtr = (char*)(*inPutWaveH) + inPutOffset;
        outPutDataStartPtr =  (char*)(*outPutWaveH) + outPutOffset;
        if (MDAccessNumericWaveData(kernelH, kMDWaveAccessMode0, &kernelOffset)) throw result = WAVEERROR_NOS;
        kernelDataStartPtr = (char*)(*kernelH) + kernelOffset;
        // make kernel table
        kernelTable = SymConvolveMakeKernelTable ((float*)kernelDataStartPtr, kernelDimSizes [0]);
        if (kernelTable == nullptr) throw result = MEMFAIL;
        // multiprocessor initialization
        nThreads = gNumProcessors;
        if (zSize < nThreads) nThreads = zSize;
        threadFrames = zSize/nThreads;
        // make frame-sized buffers for threads
        bufferPtr = (char**)WMNewPtr(sizeof(char*) * nThreads);
        if (bufferPtr == nullptr) throw result = MEMFAIL;
        for (iThread = 0 ; iThread < nThreads; iThread +=1){
            switch (inPutWaveType) {
                case NT_I64 | NT_UNSIGNED:
                case NT_I64:
                case NT_FP64:
                    bufferPtr [iThread] = (char*)WMNewPtr (inPutDimSizes[ROWS] * inPutDimSizes[COLUMNS]  * 8);
                    if (bufferPtr [iThread] == nullptr) throw result = MEMFAIL;
                    break;
                case NT_I32 | NT_UNSIGNED:
                case NT_I32:
                case NT_FP32:
                    bufferPtr[iThread] = (char*)WMNewPtr (inPutDimSizes[ROWS] * inPutDimSizes[COLUMNS] * 4);
                    if (bufferPtr [iThread] == nullptr) throw result = MEMFAIL;
                    break;
                case NT_I16 | NT_UNSIGNED:
                case NT_I16:
                    bufferPtr [iThread] = (char*)WMNewPtr (inPutDimSizes[ROWS] * inPutDimSizes[COLUMNS]   * 2);
                    if (bufferPtr [iThread] == nullptr) throw result = MEMFAIL;
                    break;
                case NT_I8 | NT_UNSIGNED:
                case NT_I8:
                    bufferPtr [iThread] = (char*)WMNewPtr (inPutDimSizes[ROWS] * inPutDimSizes[COLUMNS]  * 1);
                    if (bufferPtr [iThread] == nullptr) throw result = MEMFAIL;
                    break;
                default:
                    throw result = NUMTYPE;
                    break;
            }
        }
        // make an array of thread parameter structures
        paramArrayPtr= (ConvolveFramesThreadParamsPtr)WMNewPtr (nThreads * sizeof(ConvolveFramesThreadParams));
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
            
    }catch (int (result)) { // catch errors before starting threads
        if (bufferPtr != nullptr){
            for (iThread = 0 ; iThread < nThreads; iThread +=1){
                if (bufferPtr [iThread] == nullptr)
                    break;
                else
                    WMDisposePtr (bufferPtr [iThread]);
            }
            WMDisposePtr ((Ptr)bufferPtr);
        }
        if (threadsPtr != nullptr) WMDisposePtr ((Ptr)threadsPtr);
        if (paramArrayPtr != nullptr) WMDisposePtr ((Ptr)paramArrayPtr);
#ifdef _WINDOWS_
        if (threadIDsPtr != nullptr) WMDisposePtr((Ptr)threadIDsPtr);
#endif
        if (kernelTable !=nullptr) WMDisposePtr ((Ptr)kernelTable);
        WMDisposeHandle (p->outPutPath);    // free input string for output path
        p -> result = (double)(result - FIRST_XOP_ERR);
        #ifdef NO_IGOR_ERR
            return (0);
        #else
            return (result);
        #endif
    } // end of catch block
    // fill paramater array
    for (iThread = 0; iThread < nThreads; iThread++){
        paramArrayPtr[iThread].inPutWaveType = inPutWaveType;
        paramArrayPtr[iThread].inPutDataPtr = inPutDataStartPtr;
        paramArrayPtr[iThread].outPutDataPtr = outPutDataStartPtr;
        paramArrayPtr[iThread].frameBufferPtr = bufferPtr[iThread];
        paramArrayPtr[iThread].xSize = inPutDimSizes [0];
        paramArrayPtr[iThread].ySize = inPutDimSizes [1];
        paramArrayPtr[iThread].threadFrames = threadFrames;
        paramArrayPtr[iThread].startFrame = iThread * threadFrames;
        paramArrayPtr[iThread].kernelDataPtr=(float*)kernelDataStartPtr; //pointer to start of kernel
        paramArrayPtr[iThread].kWidth = kernelDimSizes[0]; //width of kernel
        paramArrayPtr[iThread].kernelTablePtr =kernelTable; // sum of kernel at start of column or line
        paramArrayPtr[iThread].isFloat = isFloat; // 0 for same type as input wave, non-zero for floating point wave
    }
    paramArrayPtr[nThreads -1].threadFrames += zSize % threadFrames;    // last thread gets any left-over frames
    
#ifdef __GNUC__
    // create the threads
    for (iThread = 0; iThread < nThreads; iThread++){
        pthread_create(&threadsPtr[iThread], NULL, SymConvolveFramesThread, (void*)&paramArrayPtr[iThread]);
    }
    // wait for threads to finish
    for (iThread = 0; iThread < nThreads; iThread++) {
        pthread_join(threadsPtr[iThread], NULL);
    }
#endif
#ifdef _WINDOWS_
    // create the threads with Default security attributes, Default stack size, KalmanThread function, a pointer to a parameter structure, Default creation flags, and a pointer to a threadId
    for (iThread = 0; iThread < nThreads; iThread++){
        threadsPtr[iThread] = CreateThread( NULL, 0, (LPTHREAD_START_ROUTINE)SymConvolveFramesThread, (LPVOID)&paramArrayPtr[iThread], 0, &threadIDsPtr[iThread]);
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

    // free frameBuffers
    for (iThread = 0 ; iThread < nThreads; iThread +=1){
        WMDisposePtr ((Ptr)bufferPtr [iThread]);
    }
    WMDisposePtr ((Ptr)bufferPtr);
    //free kernel table memory
    WMDisposePtr ((Ptr)kernelTable);
    WMDisposeHandle (p->outPutPath);
    // Inform Igor that we have changed the output wave.
    WaveHandleModified(outPutWaveH);
    p -> result = (0);
    return (0);
}


/* -------------------------------------- MedianFrames-------------------------------------------------------
 Applies a median filter to a 2D or 3D wave, treating each plane as a separate 2D image
 Does multiple frames in a single 3D wave
 -------------------------------------------------------------------------------------------------------------*/

/* ********************************************** medianT **********************************************
 template for a stand-alone median function for any data type
 n = number of data points, dataStrtPtr is pointer to start of data
 Last Modified 2013/07/16 by Jamie Boyd */
template <typename T> T medianT (UInt32 n, T* dataStrtPtr) {
    UInt32 medianOffset = n/2; // point offset to location of median value when function returns
    UInt32 i,right = n -1,j,left =0, mid; // point offsets used for partitioning
    T a; //
    T temp; // used for SWAP macro
    for (;;) {
        if (right <= left+1) {
            if (right == left+1 && *(dataStrtPtr + right) < *(dataStrtPtr + left)) {
                SWAP((*(dataStrtPtr + left)),(*(dataStrtPtr + right)));
            }
            return *(dataStrtPtr + medianOffset);
        } else {
            mid=(left + right) >> 1;
            SWAP((*(dataStrtPtr + mid)),(*(dataStrtPtr + left + 1)));
            if (*(dataStrtPtr + left) > *(dataStrtPtr + right)) {
                SWAP((*(dataStrtPtr + left)),(*(dataStrtPtr + right)));
            }
            if (*(dataStrtPtr + left + 1) > *(dataStrtPtr + right)) {
                SWAP((*(dataStrtPtr + left + 1)), ( *(dataStrtPtr + right)));
            }
            if (*(dataStrtPtr + left) > *(dataStrtPtr + left+1)) {
                SWAP((*(dataStrtPtr + left)),(*(dataStrtPtr + left+1)));
            }
            i=left+1;
            j=right;
            a=*(dataStrtPtr + left+1);
            for (;;) {
                do i++; while (*(dataStrtPtr + i) < a);
                do j--; while (*(dataStrtPtr + j) > a);
                if (j < i) break;
                SWAP((*(dataStrtPtr + i)),(*(dataStrtPtr + j)));
            }
            *(dataStrtPtr + left + 1) = *(dataStrtPtr + j);
            *(dataStrtPtr + j) = a;
            if (j >= medianOffset) right=j-1;
            if (j <= medianOffset) left=i;
        }
    }
}


 /* ********************************************** MedianFramesT **********************************************
  template for applying a median filter and putting the results in an output wave.
  Input wave can be 2 or 3D, but each plane is done as a separate 2D image.
 Last Modified 2025/06/24 by Jamie Boyd */
template <typename TI> void MedianFramesT (TI* srcWave, TI* destWave, TI* bufferStartPtr, CountInt xSize, CountInt ySize, CountInt zSize, UInt16 kWidth){
    CountInt fSize = (xSize * ySize);  //frame size
    CountInt bufferSize =(fSize *  sizeof (TI));
    UInt32 kSize = (kWidth * kWidth); // number of pixels in the kernel
    TI* kernelCopyStart; // pointer to start of the buffer to store copied data to be medianed in a destructive fashion
    UInt32  kBufferSize = (kSize *  sizeof (TI));
    UInt32 kRadW = (kWidth - 1)/2; //radius of the kernel width, including the central pixel
    UInt8 isOverWriting = 0; // will be set to 1 if overwriting existing wave
    TI* inPutPtr; // Center pixel input pixel. Offsets to kernel input are calculated relative to this. Value progresses monotonically through wave
    TI* convoPtr; // pointer to pixel of the image used when copying data
    TI* outPutPtr; // output pixel to get the calculated median value. Value progresses monotonically through the wave
    TI* kPtr;  // pointer to iterate through the copied data
    CountInt kX, kY, kXend, kYend; // variables for iterating through pieces of kernel
    CountInt wToNextRow; //amount to add to wToFirstRow to get to next row in input wave when iterating through a kernel
    CountInt wX, wY, wZ, wXend, wYend; // to keep track of progress through X and Y in input wave
    // overwriting source if dest == src
    if ((TI*)destWave == (TI*)srcWave){
        isOverWriting = 1;
    }
    // Make a buffer big enough to hold a "kernel's worth" of pixels - allocate in thread because it is too small to fail
    kernelCopyStart = (TI*) WMNewPtr (kBufferSize);
    // Loop through all frames
    for (wZ=0, outPutPtr = destWave, inPutPtr =srcWave; wZ < zSize; wZ++){
        if (isOverWriting){  // copy next frame to buffer and position input pointer at start of buffer
            memcpy ((void*) bufferStartPtr, (void*) outPutPtr, bufferSize);
            inPutPtr = (TI*) bufferStartPtr;
        } // If not copying to a buffer, inputPtr will be left in correct position at end of Z loop
        // Loop for TOP
        for (wY = 0; wY < kRadW; wY++){
            // loop for TOP LEFT
            for (wX =0; wX < kRadW; wX++, inPutPtr++, outPutPtr++){
                convoPtr = inPutPtr - (wY  * xSize) - wX;
                wToNextRow =  (UInt32)(xSize - kRadW - wX - 1);
                kPtr = kernelCopyStart;
                // loop through kernel at each location
                for (kY = kRadW - wY; kY < kWidth; kY++ , convoPtr +=wToNextRow){
                    for (kX = kRadW - wX; kX < kWidth; kX++, kPtr++, convoPtr++) *kPtr = *convoPtr;
                }
                *outPutPtr = medianT ((CountInt)(kPtr - kernelCopyStart), kernelCopyStart);
            } // End of Loop for TOP LEFT
            // Loop for TOP CENTER
            for (wXend =(xSize - kRadW); wX < wXend; wX++, inPutPtr++, outPutPtr++){
                // loop through kernel at each location
                convoPtr = inPutPtr - (wY  * xSize) - kRadW;
                wToNextRow = xSize - kWidth;
                kPtr = kernelCopyStart;
                for (kY = kRadW - wY; kY < kWidth; kY++, convoPtr += wToNextRow){
                    for (kX =0; kX < kWidth; kX++, kPtr++, convoPtr++) *kPtr = *convoPtr;
                }
                *outPutPtr = medianT ((CountInt)(kPtr - kernelCopyStart), kernelCopyStart);
            } // End of Loop for TOP CENTER
            // Loop for TOP RIGHT
            for (;wX < xSize; wX++, inPutPtr++, outPutPtr++){
                // loop through kernel at each location
                convoPtr =  inPutPtr - (wY  * xSize) - kRadW;
                wToNextRow = xSize - kRadW - (xSize - wX);
                kPtr = kernelCopyStart;
                for (kY = kRadW - wY; kY < kWidth; kY++, convoPtr += wToNextRow){
                    for (kX =0, kXend = kRadW + (xSize - wX); kX < kXend; kX++, kPtr++, convoPtr++) *kPtr = *convoPtr;
                }
                *outPutPtr = medianT ((CountInt)(kPtr - kernelCopyStart), kernelCopyStart);
            } // End of Loop for TOP RIGHT
        } // End of Loop for TOP
        // Loop for MIDDLE
        for (wYend =(ySize - kRadW) ; wY < wYend ; wY++){
            // Loop for MIDDLE LEFT
            for (wX =0; wX < kRadW; wX++, inPutPtr++, outPutPtr++){
                convoPtr = inPutPtr - (kRadW  * xSize) - wX;
                wToNextRow =  xSize - kRadW - wX - 1;
                kPtr = kernelCopyStart;
                // loop through kernel at each location
                for (kY = 0; kY < kWidth; kY++, convoPtr +=wToNextRow){
                    for (kX = kRadW - wX; kX < kWidth; kX++, kPtr++, convoPtr++) *kPtr = *convoPtr;
                }
                *outPutPtr = medianT ((CountInt)(kPtr - kernelCopyStart), kernelCopyStart);
            } // End of Loop for MIDDLE LEFT
            // Loop for MIDDLE CENTER
            for (wXend = (xSize - kRadW); wX < wXend ; wX++, inPutPtr++, outPutPtr++){
                convoPtr = inPutPtr - (kRadW  * xSize) - kRadW;
                wToNextRow = xSize - kWidth;
                kPtr = kernelCopyStart;
                // loop through kernel at each location
                for (kY = 0; kY < kWidth; kY++, convoPtr += wToNextRow){
                    for (kX =0; kX < kWidth; kX++, kPtr++, convoPtr++) *kPtr = *convoPtr;
                }
                *outPutPtr = medianT ((CountInt)(kPtr - kernelCopyStart), kernelCopyStart);
            } // End of Loop for MIDDLE CENTER
            // Loop for MIDDLE RIGHT
            for (;wX < xSize; wX++, inPutPtr++, outPutPtr++){
                // loop through kernel at each location
                convoPtr =  inPutPtr - (kRadW  * xSize) - kRadW;
                wToNextRow = xSize - kRadW - (xSize - wX);
                kPtr = kernelCopyStart;
                for (kY = 0; kY < kWidth; kY++, convoPtr += wToNextRow){
                    for (kX = 0, kXend = kRadW + (xSize - wX); kX < kXend; kX++, kPtr++, convoPtr++) *kPtr = *convoPtr;
                }
                *outPutPtr = medianT ((CountInt)(kPtr - kernelCopyStart), kernelCopyStart);
            } // End of Loop for MIDDLE RIGHT
        } // End of Loop for MIDDLE
        // Loop for BOTTOM
        for (; wY < ySize; wY++){
            // Loop for BOTTOM LEFT
            for (wX = 0; wX < kRadW; wX ++, inPutPtr++, outPutPtr++){
                // Loop through kernel at each location
                convoPtr = inPutPtr - (kRadW  * xSize) - wX;
                wToNextRow =  xSize - kRadW - wX - 1;
                kPtr = kernelCopyStart;
                for (kY = 0, kYend = kRadW + (ySize - wY); kY < kYend; kY++, convoPtr += wToNextRow){
                    for (kX =kRadW - wX; kX < kWidth; kX++, kPtr++, convoPtr++) *kPtr = *convoPtr;
                }
                *outPutPtr = medianT ((CountInt)(kPtr - kernelCopyStart), kernelCopyStart);
            } // End ofLoop for BOTTOM LEFT
            // Loop for BOTTOM CENTER
            for (; wX < (xSize - kRadW); wX++, inPutPtr++, outPutPtr++){
                // Loop through kernel at each location
                convoPtr = inPutPtr - (kRadW  * xSize) - kRadW;
                wToNextRow = xSize - kWidth;
                kPtr = kernelCopyStart;
                for (kY = 0, kYend = kRadW + (ySize - wY) ; kY < kYend; kY++, convoPtr += wToNextRow){
                    for (kX =0; kX < kWidth; kX++, convoPtr++, kPtr++) *kPtr = *convoPtr;
                }
                *outPutPtr = medianT ((CountInt)(kPtr - kernelCopyStart), kernelCopyStart);
            } // End of Loop for BOTTOM CENTER
            // Loop for BOTTOM RIGHT
            for (; wX < xSize; wX++, inPutPtr++, outPutPtr++){
                //Loop through kernel at each location
                convoPtr =  inPutPtr - (kRadW * xSize) - kRadW;
                wToNextRow = xSize - kRadW - (xSize - wX);
                kPtr = kernelCopyStart;
                for (kYend = kRadW + (ySize - wY), kY =0; kY < kYend; kY++, convoPtr += wToNextRow){
                    for (kX =0, kXend = kRadW + (xSize -wX); kX < kXend; kX++, convoPtr++, kPtr++) *kPtr = *convoPtr;
                }
                *outPutPtr = medianT ((CountInt)(kPtr - kernelCopyStart), kernelCopyStart);
            } // End of Loop for BOTTOM RIGHT
        } // End of Loop for BOTTOM
    } // End of Loop for Each Frame
    if (kernelCopyStart != nullptr) WMDisposePtr ((Ptr)kernelCopyStart);
}

/* ***************************** MedianFramesThreadParams *****************************
 Structure to pass data to each MedianFramesThread
 Last Modified 2013/07/16 by Jamie Boyd */
typedef struct MedianFramesThreadParams{
    int inPutWaveType;
    char* inPutDataStartPtr;
    char* outPutDataStartPtr;
    char* bufferPtr;
    CountInt xSize;
    CountInt ySize;
    CountInt zSize;
    CountInt startFrame;        // frame to start on for this threrad
    CountInt threadFrames;      // number of frames to do for this thread
    UInt16 kWidth;
}MedianFramesThreadParams, *MedianFramesThreadParamsPtr;


/* ***************************** MedianFramesThread *****************************
 Each thread to median filter a range of frames starts with this function
 Last Modified 2013/07/16 by Jamie Boyd */
void* MedianFramesThread (void* threadarg){
    struct MedianFramesThreadParams* p;
    p = (struct MedianFramesThreadParams*) threadarg;
    CountInt frameSize = p->xSize * p->ySize;
    CountInt startPos  = p->startFrame * frameSize; //change start position from frames to data points by multiplying by frame size
    // call the right template function for the wave types
    switch (p->inPutWaveType) {
        case NT_I8:
            MedianFramesT ((char*)p->inPutDataStartPtr + startPos, (char*)p->outPutDataStartPtr + startPos, (char*)p->bufferPtr, p->xSize, p->ySize, p->threadFrames, p->kWidth);
            break;
        case (NT_I8 | NT_UNSIGNED):
            MedianFramesT ((unsigned char*)p->inPutDataStartPtr + startPos,(unsigned char*)p->outPutDataStartPtr + startPos, (unsigned char*)p->bufferPtr, p->xSize, p->ySize, p->threadFrames, p->kWidth);
            break;
        case NT_I16:
            MedianFramesT ((short*)p->inPutDataStartPtr + startPos,(short*)p->outPutDataStartPtr + startPos, (short*)p->bufferPtr, p->xSize, p->ySize, p->threadFrames, p->kWidth);
            break;
        case (NT_I16 | NT_UNSIGNED):
            MedianFramesT ((unsigned short*)p->inPutDataStartPtr + startPos,(unsigned short*)p->outPutDataStartPtr + startPos, (unsigned short*)p->bufferPtr, p->xSize, p->ySize, p->threadFrames, p->kWidth);
            break;
        case NT_I32:
            MedianFramesT ((SInt32*)p->inPutDataStartPtr + startPos,(SInt32*)p->outPutDataStartPtr + startPos, (SInt32*)p->bufferPtr, p->xSize, p->ySize, p->threadFrames, p->kWidth);
            break;
        case (NT_I32| NT_UNSIGNED):
            MedianFramesT ((UInt32*)p->inPutDataStartPtr + startPos,(UInt32*)p->outPutDataStartPtr + startPos, (UInt32*)p->bufferPtr, p->xSize, p->ySize, p->threadFrames, p->kWidth);
            break;
        case NT_I64:
            MedianFramesT((SInt64*)p->inPutDataStartPtr + startPos, (SInt64*)p->outPutDataStartPtr + startPos, (SInt64*)p->bufferPtr, p->xSize, p->ySize, p->threadFrames, p->kWidth);
            break;
        case (NT_I64 | NT_UNSIGNED):
            MedianFramesT((UInt64*)p->inPutDataStartPtr + startPos, (UInt64*)p->outPutDataStartPtr + startPos, (UInt64*)p->bufferPtr, p->xSize, p->ySize, p->threadFrames, p->kWidth);
            break;
        case NT_FP32:
            MedianFramesT ((float*)p->inPutDataStartPtr + startPos,(float*)p->outPutDataStartPtr + startPos, (float*)p->bufferPtr, p->xSize, p->ySize, p->threadFrames, p->kWidth);
            break;
        case NT_FP64:
            MedianFramesT ((double*)p->inPutDataStartPtr + startPos,(double*)p->outPutDataStartPtr + startPos, (double*)p->bufferPtr, p->xSize, p->ySize, p->threadFrames, p->kWidth);
            break;
    }
    return nullptr;
}

/* ***************************** MedianFrames *****************************
 MedianFrames XOP entry function
 Applies a median filter to a 2D or 3D wave, treating each plane in a 3D wave as a separate 2D image
 typedef struct MedianFramesParams{
 double overWrite; // 1 if it is o.k. to overwrite existing waves, 0 to exit with error
 double kWidth; //width of the area over which to calculate the median. Must be an odd number
 Handle outPutPath;	// A handle to a string containing path to output wave we want to make, or empty string to overwrite existing wave
 waveHndl inPutWaveH;//input wave. needs to be 2D or 3D wave
 double result;
 } MedianFramesParams, *MedianFramesParamsPtr;
 Last modified 2013/07/16 by Jamie Boyd */
extern "C" int MedianFrames (MedianFramesParamsPtr p){
    int result = 0;                             // The error returned from various Wavemetrics functions
    waveHndl inPutWaveH, outPutWaveH;		    // handles to the input and output waves
	int inPutWaveType;                          //  Wavetypes numeric codes for things like 32 bit floating point, 16 bit int, etc
    int inPutDimensions;	                    // number of dimensions in input wave
    CountInt inPutDimSizes[MAX_DIMENSIONS+1];	// an array used to hold the width, height, layers, and chunk sizes
    int outPutWaveType;
    CountInt outPutDimSizes[MAX_DIMENSIONS+1];
    int outPutDimensions;
    CountInt zSize, frameSize;
    BCInt inPutOffset, outPutOffset;	        //offset in bytes from begnning of handle to a wave to the actual data
    DataFolderHandle outPutDFHandle;	// Handle to the datafolder where we will put the output wave
    DFPATH outPutPath;	// string to hold data folder path of input wave
    WVNAME outPutWaveName;	        // C strings to hold names of input and output waves
    UInt8 overWriteOK = (UInt8)(p->overWrite);	    // 0 to not overwrite output wave if it already exists, 1 to overwrite old waves
    UInt16 kWidth = (UInt16)(p-> kWidth);           // width of the area over which the median will be calculated
    UInt8 isOverWritingInput;                       // non-zero if output is overwriting input wave
    char *inPutDataStartPtr = nullptr;
    char *outPutDataStartPtr = nullptr;
    char **bufferPtr = nullptr;
    // for threads
    UInt8 iThread, nThreads;
    CountInt threadFrames;
    MedianFramesThreadParamsPtr paramArrayPtr = nullptr;
#ifdef __GNUC__
    pthread_t* threadsPtr = nullptr;   // pointer to an array of pThread_t structures we will make
#endif
#ifdef  _WINDOWS_
    HANDLE *threadsPtr = nullptr;      // pointer to an array of Windows handles, which will point to thread structures we make
    DWORD *threadIDsPtr = nullptr;     // pointer to an array of DWORD (32-bit unsigned integer) that Windows uses for thread IDs
#endif
    
    try{
        // Check that kWidth is odd
        if ((kWidth % 2) == 0) throw result = BADKERNEL;
        // Get handle to input wave
        inPutWaveH = p->inPutWaveH;
        if (inPutWaveH == NIL) throw result = NON_EXISTENT_WAVE;
        // Get wave data type
        inPutWaveType = WaveType(inPutWaveH);
        if (inPutWaveType==TEXT_WAVE_TYPE) throw result = NOTEXTWAVES;
        // Get number of used dimensions in wave
        if (MDGetWaveDimensions(inPutWaveH, &inPutDimensions, inPutDimSizes)) throw result = WAVEERROR_NOS;
        // Check that inputwave is 2D or 3D
        if ((inPutDimensions == 1) || (inPutDimensions == 4)) throw result = INPUTNEEDS_2D3D_WAVE;
        frameSize = inPutDimSizes [ROWS] * inPutDimSizes [COLUMNS];
        // if z size is 0, make it 1 to calculate size
        if (inPutDimSizes [LAYERS] == 0)
            zSize = 1;
        else
            zSize=inPutDimSizes [LAYERS];
        // If outPutPath is empty string, we are overwriting existing wave
        if (WMGetHandleSize (p->outPutPath) == 0){
            if (!(overWriteOK)) throw result = OVERWRITEALERT;
            outPutWaveH = inPutWaveH;
            isOverWritingInput = 1;
        }else{ // Parse outPut path
            ParseWavePath (p->outPutPath, outPutPath, outPutWaveName);
            //No liberal wave names for output wave, maybe we don't really need to enforce non-liberal names?
            CleanupName (0, outPutWaveName, MAX_OBJ_NAME);
            // Get a handle to the output datafolder
            if (GetNamedDataFolder (NULL, outPutPath, &outPutDFHandle))throw result = WAVEERROR_NOS;
            // try to get a handle to the output wave specified in the string
            outPutWaveH = FetchWaveFromDataFolder(outPutDFHandle, outPutWaveName);
            if (outPutWaveH == nullptr){        // output wave specified in wave name does not exist, so make it
                if (MDMakeWave (&outPutWaveH, outPutWaveName, outPutDFHandle, inPutDimSizes, inPutWaveType, overWriteOK)) throw result = WAVEERROR_NOS;
            }else{  // wave specified in the input string already exists and outPutWaveH now points to it
                if (!(overWriteOK)) throw result = OVERWRITEALERT;  // because output wave already exists
                if (outPutWaveH == inPutWaveH){ // overwriting input wave
                    isOverWritingInput = 1;
                }else{ // outPut wave is an existing wave, and is not the input wave, redimension and change type if needed
                    outPutWaveType = WaveType(outPutWaveH);
                    if (outPutWaveType==TEXT_WAVE_TYPE) throw result = NOTEXTWAVES;
                    //Get number of dimensions and dimension sizes in output wave
                    if (MDGetWaveDimensions(outPutWaveH, &outPutDimensions, outPutDimSizes))throw result = WAVEERROR_NOS;
                    if ((outPutWaveType |= inPutWaveType) || (outPutDimSizes [0] |= inPutDimSizes [0]) || (outPutDimSizes [1] |= inPutDimSizes [1]) || (outPutDimSizes [2] |= inPutDimSizes [2])){
                        if (MDChangeWave(outPutWaveH, inPutWaveType, inPutDimSizes)) throw result = WAVEERROR_NOS;
                    }
                }
            }
        }
        
		//Get data offsets for the 2 waves (1 wave, if overwriting)
        if (MDAccessNumericWaveData(inPutWaveH, kMDWaveAccessMode0, &inPutOffset)) return WAVEERROR_NOS;
        if (isOverWritingInput){
            outPutOffset = inPutOffset;
        }else{
            if (MDAccessNumericWaveData(outPutWaveH, kMDWaveAccessMode0, &outPutOffset)) throw result = WAVEERROR_NOS;
        }
        inPutDataStartPtr = (char*)(*inPutWaveH) + inPutOffset;
        outPutDataStartPtr =  (char*)(*outPutWaveH) + outPutOffset; // true even if overwriting
        // multiprocessor init
        nThreads = gNumProcessors;
        if (zSize < nThreads) nThreads = zSize;
        threadFrames = zSize/nThreads;
       
        // make a buffer of frames for threads if overwriting src
        if (isOverWritingInput){
            // make frame-sized buffers for threads
            bufferPtr = (char**)WMNewPtr(sizeof(char*) * nThreads);
            if (bufferPtr == nullptr) throw result = MEMFAIL;
            for (iThread = 0 ; iThread < nThreads; iThread +=1){
                switch (inPutWaveType) {
                    case NT_I64 | NT_UNSIGNED:
                    case NT_I64:
                    case NT_FP64:
                        bufferPtr [iThread] = (char*)WMNewPtr (inPutDimSizes[ROWS] * inPutDimSizes[COLUMNS]  * 8);
                        if (bufferPtr [iThread] == nullptr) throw result = MEMFAIL;
                        break;
                    case NT_I32 | NT_UNSIGNED:
                    case NT_I32:
                    case NT_FP32:
                        bufferPtr[iThread] = (char*)WMNewPtr (inPutDimSizes[ROWS] * inPutDimSizes[COLUMNS] * 4);
                        if (bufferPtr [iThread] == nullptr) throw result = MEMFAIL;
                        break;
                    case NT_I16 | NT_UNSIGNED:
                    case NT_I16:
                        bufferPtr [iThread] = (char*)WMNewPtr (inPutDimSizes[ROWS] * inPutDimSizes[COLUMNS]   * 2);
                        if (bufferPtr [iThread] == nullptr) throw result = MEMFAIL;
                        break;
                    case NT_I8 | NT_UNSIGNED:
                    case NT_I8:
                        bufferPtr [iThread] = (char*)WMNewPtr (inPutDimSizes[ROWS] * inPutDimSizes[COLUMNS]  * 1);
                        if (bufferPtr [iThread] == nullptr) throw result = MEMFAIL;
                        break;
                    default:
                        throw result = NUMTYPE;
                        break;
                }
            }
        }
        // make array of parameter structures */
        paramArrayPtr = (MedianFramesThreadParamsPtr)WMNewPtr(nThreads * sizeof(MedianFramesThreadParams));
        if (paramArrayPtr == nullptr) throw result = NOMEM;
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
    }catch (int result){
        if (bufferPtr != nullptr){
            for (iThread = 0 ; iThread < nThreads; iThread +=1){
                if (bufferPtr [iThread] == nullptr)
                    break;
                else
                    WMDisposePtr ((Ptr)bufferPtr [iThread]);
            }
            WMDisposePtr ((Ptr)bufferPtr);
        }
        if (threadsPtr != nullptr) WMDisposePtr ((Ptr)threadsPtr);
#ifdef _WINDOWS_
        if (threadIDsPtr != nullptr) WMDisposePtr((Ptr)threadIDsPtr);
#endif
        if (paramArrayPtr != nullptr) WMDisposePtr ((Ptr)paramArrayPtr);
        WMDisposeHandle (p->outPutPath);    // free input string for output path
        p -> result = (double)(result - FIRST_XOP_ERR);
        #ifdef NO_IGOR_ERR
            return (0);
        #else
            return (result);
        #endif
    }
    // fill params array
    for (iThread = 0; iThread < nThreads; iThread++){
        paramArrayPtr[iThread].inPutWaveType = inPutWaveType;
        paramArrayPtr[iThread].inPutDataStartPtr = inPutDataStartPtr;
        paramArrayPtr[iThread].outPutDataStartPtr = outPutDataStartPtr;
        paramArrayPtr[iThread].bufferPtr = bufferPtr[iThread];
        paramArrayPtr[iThread].xSize = inPutDimSizes [0];
        paramArrayPtr[iThread].ySize = inPutDimSizes [1];
        paramArrayPtr[iThread].zSize = zSize;
        paramArrayPtr[iThread].startFrame = iThread * threadFrames;
        paramArrayPtr[iThread].threadFrames = threadFrames;
        paramArrayPtr[iThread].kWidth = (int)kWidth;
    }
    paramArrayPtr[nThreads -1].threadFrames += zSize % threadFrames;    // last thread gets any left-over frames
#ifdef __GNUC__
    // create the threads
    for (iThread = 0; iThread < nThreads; iThread++){
        pthread_create(&threadsPtr[iThread], NULL, MedianFramesThread, (void*)&paramArrayPtr[iThread]);
    }
    // wait for threads to finish
    for (iThread = 0; iThread < nThreads; iThread++) {
        pthread_join(threadsPtr[iThread], NULL);
    }
#endif
#ifdef _WINDOWS_
    // create the threads with Default security attributes, Default stack size, MedianFramesThread function, a pointer to a parameter structure, Default creation flags, and a pointer to a threadId
    for (iThread = 0; iThread < nThreads; iThread++){
        threadsPtr[iThread] = CreateThread( NULL, 0, (LPTHREAD_START_ROUTINE)MedianFramesThread, (LPVOID)&paramArrayPtr[iThread], 0, &threadIDsPtr[iThread]);
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
    if (isOverWritingInput){
           for (iThread = 0 ; iThread < nThreads; iThread +=1){
               WMDisposePtr ((Ptr)bufferPtr [iThread]);
           }
           WMDisposePtr ((Ptr)bufferPtr);
       }
    WMDisposePtr ((Ptr)threadsPtr);     // free memory for pThreads Array
    WMDisposePtr ((Ptr)paramArrayPtr);  // free memory for paramaters Array
    WMDisposeHandle (p->outPutPath);    // free input string for output path
    WaveHandleModified(outPutWaveH);	// Inform Igor that we have changed the output wave.
    p -> result = (0);
    return (0);
}
