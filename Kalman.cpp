#include "twoPhoton.h"
/* ---------------------------------------------Kalman----------------------------------------------------------
 Code for Kalman averaging of frames in a 3d wave, 2d wave, or a list of waves
 Last Modified 2025/06/23 by Jamie Boyd
 ---------------------------------------------------------------------------------------------------------------*/

/* The following template is used to handle any one of the 8 types of wave data, for any of the three Kalman averaging functions
 for 3D waves
 Last modified 2014/01/28 by Jamie Boyd */
template <typename T> int KalmanT(T *srcWaveStart, T *destWaveStart, CountInt pixPerThread, CountInt pixToNextFrame, CountInt numLayers, float multiplier) {
	T* destWaveEnd = destWaveStart + pixPerThread;	// End of the output layer we are putting the average into
	T* srcWave, *destWave;	// pointers used to iterate through source wave and destination wave
	/* If multiplier is < 1, use standard averaging across layers with a floating point temporary value */
	if (multiplier < 1){
		double tempVal;
		CountInt pixPerLayer = pixPerThread + pixToNextFrame;
		T* srcWaveEnd = srcWaveStart + (numLayers * pixPerLayer);
		CountInt toNextSrcPix = (numLayers * pixPerLayer) - 1;
		for (srcWave = srcWaveStart, destWave = destWaveStart ; destWave < destWaveEnd ; srcWave -= toNextSrcPix, destWave ++){			for (tempVal =0; srcWave < srcWaveEnd; srcWave += pixPerLayer)
			tempVal += *srcWave;
            *destWave = tempVal/numLayers;
		}
	}else{
		CountInt layer; // used to iterate through layers
		// Special stuff for first frame
		if (srcWaveStart == destWaveStart){	//this happens when collapsing a wave into the first frame
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


/* Structure to pass data to each KalmanThread
 Last Modified Feb 21 2010 by Jamie Boyd */
typedef struct KalmanThreadParams{
    int inPutWaveType;
    char* inPutDataStartPtr;
    char* outPutDataStartPtr;
    CountInt startLayer;
    CountInt outPutLayer;
    CountInt xSize;
    CountInt ySize;
    CountInt zSize;
    float multiplier;
    UInt8 ti; // number of this thread, starting from 0
    UInt8 tN; // total number of threads
} KalmanThreadParams, *KalmanThreadParamsPtr;


/* Each thread to do Kalman averaging starts with this function
 Last Modified 2013/07/16 by Jamie Boyd */
void* KalmanThread (void* threadarg){
	struct KalmanThreadParams* p;
	p = (struct KalmanThreadParams*) threadarg;
	CountInt xSize= p->xSize;
	CountInt ySize = p->ySize;
	CountInt zSize = p->zSize;
	CountInt startLayer = p->startLayer;
	CountInt outPutLayer = p->outPutLayer;
	float multiplier = p->multiplier;
	UInt8 ti = p->ti;
	UInt8 tN = p->tN;
	CountInt frameSize = ySize * xSize;
	CountInt pixPerThread = frameSize/tN;
	CountInt startPos = ti * pixPerThread; // which pixel to start this thread on depends on thread number * points per thread. ti is 0 based
	if (ti == (tN - 1)) pixPerThread += (frameSize % tN); // last thread gets any extra pixels
	CountInt pixToNextFrame = frameSize - pixPerThread;
	int result;
	switch (p->inPutWaveType) {
	case NT_I8:
		result = KalmanT ((char*)p->inPutDataStartPtr + (startLayer * frameSize) + startPos, (char*)p->outPutDataStartPtr + (outPutLayer * frameSize) + startPos, pixPerThread, pixToNextFrame, zSize, multiplier);
		break;
	case (NT_I8 | NT_UNSIGNED):
		result = KalmanT ((unsigned char*)p->inPutDataStartPtr + (startLayer * frameSize) + startPos, (unsigned char*)p->outPutDataStartPtr + (outPutLayer * frameSize) + startPos, pixPerThread, pixToNextFrame, zSize, multiplier);
		break;
	case NT_I16:
		result = KalmanT ((short*)p->inPutDataStartPtr + (startLayer * frameSize) + startPos, (short*) p->outPutDataStartPtr + (outPutLayer * frameSize) + startPos, pixPerThread, pixToNextFrame, zSize, multiplier);
		break;
	case (NT_I16 | NT_UNSIGNED):
		result = KalmanT ((unsigned short*)p->inPutDataStartPtr + (startLayer * frameSize) + startPos, (unsigned short*) p->outPutDataStartPtr + (outPutLayer * frameSize) + startPos, pixPerThread, pixToNextFrame, zSize, multiplier);
		break;
	case NT_I32:
		result = KalmanT ((long*)p->inPutDataStartPtr + (startLayer * frameSize) + startPos, (long*) p->outPutDataStartPtr + (outPutLayer * frameSize) + startPos, pixPerThread, pixToNextFrame, zSize, multiplier);
		break;
	case (NT_I32| NT_UNSIGNED):
		result = KalmanT ((unsigned long*)p->inPutDataStartPtr + (startLayer * frameSize)+ startPos, (unsigned long*) p->outPutDataStartPtr + (outPutLayer * frameSize) + startPos, pixPerThread, pixToNextFrame, zSize, multiplier);
		break;
	case NT_FP32:
		result = KalmanT ((float*)p->inPutDataStartPtr + (startLayer * frameSize) + startPos, (float*) p->outPutDataStartPtr + (outPutLayer * frameSize) + startPos, pixPerThread, pixToNextFrame, zSize, multiplier);
		break;
	case NT_FP64:
		result = KalmanT ((double*)p->inPutDataStartPtr + (startLayer * frameSize) + startPos, (double*) p->outPutDataStartPtr + (outPutLayer * frameSize) + startPos, pixPerThread, pixToNextFrame, zSize, multiplier);
		break;
	default:	// Unknown data type - possible in a future version of Igor.
		result= NT_FNOT_AVAIL;
		break;
	}
	if (result != 0) throw result;
    return 0;
}

/* KalmanAllFrames XOP entry function
 Does Kalman averaging across all layers in a 3D wave and places results in a new 2D wave
 KalmanAllFramesParams
 waveHndl inPutWaveH    handle to a 3D input wave
 Handle outPutPath      A handle to a string containing path to output wave we want to make
 double multiplier      Multiplier for,e.g., 16 bit waves containing less than 16 bits of data
 double overWrite       0 to give errors when wave already exists. non-zero to overwrite existing wave without warning.
 result                 0 for success, error code for failure
Last Modified 2015/06/18 by Jamie Boyd  */
extern "C" int KalmanAllFrames(KalmanAllFramesParamsPtr p) {
	int result = 0;	// The error returned from various Wavemetrics functions
	waveHndl inPutWaveH = NULL, outPutWaveH = NULL;	// Handles to the input and output waves
    int inPutWaveType;	// Wavemetrics numeric code for data type of wave
	int inPutDimensions;		// The number of dimensions used in the input wave
    CountInt inPutOffset, outPutOffset; //offset in bytes from begnning of handle to a wave to the actual data - size of headers, units, etc.
	CountInt inPutDimensionSizes[MAX_DIMENSIONS+1];	// An array used to hold the sizes of each dimension of the input wave
	UInt16 outPutPathLen; //Length of the path to the target folder (output path:wave name)
	DataFolderHandle inPutDFHandle, outPutDFHandle;	// Handle to the datafolder where we will put the output wave
	DFPATH inPutPath, outPutPath;	// string to hold data folder path of input wave
	WVNAME inPutWaveName, outPutWaveName;	// C string to hold name of input wave
	UInt8 overWrite = (int)p->overWrite;	// 0 to not overwrite output wave if it already exists, 1 to overwrite old waves
	UInt8 isOverWriting; // non-zero if output is overwriting input wave
	char *inPutDataStartPtr, *outPutDataStartPtr;
	float multiplier = (float)(p->multiplier);	// multiplier for integer waves containing less than their full range of data
	CountInt zSize;
	UInt8 iThread, nThreads;
    KalmanThreadParamsPtr paramArrayPtr = NULL;
    pthread_t* threadsPtr = NULL;
    try {
		// Get handle to input wave.
		inPutWaveH = p ->inPutWaveH;
		if(inPutWaveH == NIL)throw result = NON_EXISTENT_WAVE;
		// get wave data type and check that we don't have a text wave
		inPutWaveType = WaveType(inPutWaveH);
		if (inPutWaveType==TEXT_WAVE_TYPE) throw result = NOTEXTWAVES;
		//Get number of used dimensions in input wave.
		if (MDGetWaveDimensions(inPutWaveH, &inPutDimensions, inPutDimensionSizes))throw result = WAVEERROR_NOS;
		// Check that input wave is 3D
		if (inPutDimensions != 3)
			throw result = INPUTNEEDS_3D_WAVE;
		// Save z Size as we will be resizing dimensions array
		zSize = inPutDimensionSizes [LAYERS];
		// If outPutPath is empty string, we are overwriting existing wave
		outPutPathLen = WMGetHandleSize (p->outPutPath);
		if (outPutPathLen == 0){
			if (overWrite == NO_OVERWITE) throw result = OVERWRITEALERT;
			outPutWaveH = inPutWaveH;
			isOverWriting = 1;
		}else{ // Parse outPut path
			ParseWavePath (p->outPutPath, outPutPath, outPutWaveName);
			//check that data folder is valid and get a handle to the datafolder
			if (GetNamedDataFolder (NULL, outPutPath, &outPutDFHandle))throw result = WAVEERROR_NOS;
			// Test name and data folder for output wave against the input wave to prevent accidental overwriting, if src and dest are the same
			WaveName (inPutWaveH, inPutWaveName);
            if (GetWavesDataFolder (inPutWaveH, &inPutDFHandle)) throw result = WAVEERROR_NOS;
            if (GetDataFolderNameOrPath (inPutDFHandle, 1, inPutPath)) throw result = WAVEERROR_NOS;
			if ((!(CmpStr (inPutPath,outPutPath))) && (!(CmpStr (inPutWaveName,outPutWaveName)))){	// Then we would overwrite wave
				isOverWriting = 1;
				if (overWrite == NO_OVERWITE) throw result = OVERWRITEALERT;
				outPutWaveH = inPutWaveH;
			}else{
				isOverWriting = 0;
				// make the output wave
				//Chop off the layers dimension, as we will use dimensionSizes array to make 2D output wave. We pass numlayers in a separate variable
				inPutDimensionSizes [2] = 0;
				//No liberal wave names for output wave
				CleanupName (0, outPutWaveName, MAX_OBJ_NAME);
				if (MDMakeWave (&outPutWaveH, outPutWaveName, outPutDFHandle, inPutDimensionSizes, inPutWaveType, overWrite)) throw result = WAVEERROR_NOS;
			}
		}
        //Get data offsets for the 2 waves (1 wave, if overwriting)
		if (MDAccessNumericWaveData(inPutWaveH, kMDWaveAccessMode0, &inPutOffset) != 0) throw result = WAVEERROR_NOS;
		inPutDataStartPtr = (char*)(*inPutWaveH) + inPutOffset;
		if (isOverWriting){
			outPutOffset = inPutOffset;
			outPutDataStartPtr = inPutDataStartPtr;
		}else{
			if (MDAccessNumericWaveData(outPutWaveH, kMDWaveAccessMode0, &outPutOffset)) throw result = WAVEERROR_NOS;
			outPutDataStartPtr =  (char*)(*outPutWaveH) + outPutOffset;
		}
        // get ready for multiprocessor initialization
        if (zSize < gNumProcessors){
            nThreads =zSize;
        }else{
            nThreads = gNumProcessors;
        }
        // make array of parameter structures
        paramArrayPtr = (KalmanThreadParamsPtr)WMNewPtr (nThreads * sizeof(KalmanThreadParams));
        if (paramArrayPtr == NULL) throw result = MEMFAIL;
        // make an array of pthread_t
        threadsPtr =(pthread_t*)WMNewPtr(nThreads * sizeof(pthread_t));
        if (threadsPtr == NULL) throw result = MEMFAIL;
        // catch error before starting threads
	}catch (int result){
        if (paramArrayPtr != NULL) WMDisposePtr ((Ptr)paramArrayPtr);
        if (threadsPtr != NULL) WMDisposePtr ((Ptr)paramArrayPtr);
        WMDisposeHandle(p->outPutPath);  // dispose passed in string paramater
        p -> result = (double)(result - FIRST_XOP_ERR);
#ifdef NO_IGOR_ERR
		return (0);
#else
        return (result);
#endif
	}
    // fill paramater array
    for (iThread = 0; iThread < nThreads; iThread++){
        paramArrayPtr[iThread].inPutWaveType = inPutWaveType;
        paramArrayPtr[iThread].inPutDataStartPtr = inPutDataStartPtr;
        paramArrayPtr[iThread].outPutDataStartPtr = outPutDataStartPtr;
        paramArrayPtr[iThread].startLayer = 0;
        paramArrayPtr[iThread].outPutLayer =0;
        paramArrayPtr[iThread].xSize = inPutDimensionSizes [0];
        paramArrayPtr[iThread].ySize = inPutDimensionSizes [1];
        paramArrayPtr[iThread].zSize =zSize;
        paramArrayPtr[iThread].multiplier =multiplier;
        paramArrayPtr[iThread].ti=iThread; // number of this thread, starting from 0
        paramArrayPtr[iThread].tN =nThreads; // total number of threads
    }
    // create the threads
    for (iThread = 0; iThread < nThreads; iThread++){
        pthread_create (&threadsPtr[iThread], NULL, KalmanThread, (void *) &paramArrayPtr[iThread]);
    }
    // Wait till all the threads are finished
    for (iThread = 0; iThread < nThreads; iThread++){
        pthread_join (threadsPtr[iThread], NULL);
    }
    WMDisposePtr ((Ptr)threadsPtr);      // free memory for pThreads Array
    WMDisposePtr ((Ptr)paramArrayPtr);  // Free paramaterArray memory
    WMDisposeHandle(p->outPutPath);    // dispose passed in string paramater
    if (isOverWriting){	//then collapsing a 3D wave to 2 D
        inPutDimensionSizes [0] = -1;
        inPutDimensionSizes [1] = -1;
        inPutDimensionSizes [2] = 0;
        inPutDimensionSizes [3] = 0;
        MDChangeWave (outPutWaveH, -1, inPutDimensionSizes);
    }
    // Inform Igor that we have changed output wave
    WaveHandleModified(outPutWaveH);
    p -> result = (0);
    return (0);
}

/* KalmanSpecFrames XOP entry function
 Averages a specified range of layers of the input wave into a specified layer of the output wave
 KalmanSpecFramesParams
 inPutWaveH          handle to input wave
 startLayer          start of layers to average for input wave
 endLayer            end of layers to average
 outPutWaveH         handle to output wave
 outPutLayer         layer of output wave to receive results of averaging
 multiplier          Multiplier, as for 16 bit waves containing less than 16 bits of data
 result              0 or error code
 Last Modified 2014/02/13 by Jamie Boyd */
extern "C" int KalmanSpecFrames(KalmanSpecFramesParamsPtr p) {
	int result = 0;	// The error returned from various Wavemetrics functions
	waveHndl inPutWaveH = NIL, outPutWaveH = NIL;	// Handles to the input and output waves
	int inPutWaveType, outPutWaveType;	// Wavemetrics numeric code for data type of wave
    int inPutDimensions,outPutDimensions;	// number of dimensions in input and output waves
	CountInt inPutDimensionSizes[MAX_DIMENSIONS+1], outPutDimensionSizes[MAX_DIMENSIONS+1];	// an array used to hold the width, height, layers, and chunk sizes
	CountInt inPutOffset, outPutOffset;	//offset in bytes from begnning of handle to a wave to the actual data - size of headers, units, etc.
	CountInt startLayer, endLayer, layersToDo;  //vaiables for iterating through the data.
	CountInt outPutLayer;	//The layer of the output wave that gets the result
	CountInt pointsPerLayer;	// The number of points in a layer, needed information for iterating through a layer
	char *inPutDataStartPtr, *outPutDataStartPtr;
	float multiplier = (float)(p->multiplier);	// multiplier for integer waves containing less than their full range of data
    UInt8 iThread, nThreads;
    KalmanThreadParamsPtr paramArrayPtr = NULL;
    pthread_t* threadsPtr = NULL;
    // try/catch before starting threads
	try {
		// Get handles to input wave and kernel. Make sure both waves exist.
		inPutWaveH = p ->inPutWaveH;	// Get Handle to the input Wave and make sure it exists
		outPutWaveH =  p->outPutWaveH;
		if ((inPutWaveH == NIL) || (outPutWaveH == NIL))throw result = NON_EXISTENT_WAVE;
		// get wave data type and check that datatypes are the same and that neither is a text wave
		inPutWaveType = WaveType(inPutWaveH);
		outPutWaveType = WaveType(outPutWaveH);
		if (inPutWaveType != outPutWaveType)throw result = NOTSAMEWAVETYPE;
		if ((inPutWaveType==TEXT_WAVE_TYPE) || (outPutWaveType==TEXT_WAVE_TYPE))throw result = NOTEXTWAVES;
		// Get number of used dimensions in waves.
		if (MDGetWaveDimensions(inPutWaveH, &inPutDimensions, inPutDimensionSizes))throw result=WAVEERROR_NOS;
		if (MDGetWaveDimensions(outPutWaveH, &outPutDimensions, outPutDimensionSizes))throw result=WAVEERROR_NOS;
		// Check that input wave is 3D and output wave is 2D or 3D
		if (inPutDimensions != 3)throw result = INPUTNEEDS_3D_WAVE;
		if (!((outPutDimensions == 2) || (outPutDimensions == 3))) throw result = OUTPUTNEEDS_2D3D_WAVE;
		//	Check that X and Y dimensions of the 2 waves are the same size.
		if (!((inPutDimensionSizes[0] == outPutDimensionSizes [0]) && (inPutDimensionSizes[1] == outPutDimensionSizes [1]))) throw result = NOTSAMEDIMSIZE;
		// Load outputlayer and the input startlayer and endlayer into local variables and check that they are o.k. wrt number of frames
		// InPut layer must be 0 to use 2 D wave as output wave
		outPutLayer = (long)p ->outPutLayer;
		if (((outPutLayer != 0) && (outPutLayer > outPutDimensionSizes [2] - 1)) || (outPutLayer < 0)) throw result = INVALIDOUTPUTFRAME;
		startLayer = p -> startLayer;
		if (startLayer > inPutDimensionSizes [2] -1)throw result = INVALIDINPUTFRAME;
        // Clip start layer to first layer, if start layer < 0
        if (startLayer < 0) startLayer = 0;
		endLayer = p -> endLayer;
		if (endLayer < startLayer){
            CountInt temp;
            SWAP(startLayer, endLayer);
        }
		// Clip endlayer to the last layer of the input wave
		if (endLayer > inPutDimensionSizes [2] -1) endLayer = inPutDimensionSizes [2] -1;
		// Calculate number of layers to do
		layersToDo = endLayer - startLayer + 1;
		// As X and Y are same size for input and output waves, we need only look at the input wave to get points per layer
		pointsPerLayer = inPutDimensionSizes[0] * inPutDimensionSizes[1];
        //Get data offsets for the 2 waves
		if (MDAccessNumericWaveData(inPutWaveH, kMDWaveAccessMode0, &inPutOffset)) throw result = WAVEERROR_NOS;
		inPutDataStartPtr = (char*)(*inPutWaveH) + inPutOffset;
		if (MDAccessNumericWaveData(outPutWaveH, kMDWaveAccessMode0, &outPutOffset))throw result = WAVEERROR_NOS;
		outPutDataStartPtr =  (char*)(*outPutWaveH) + outPutOffset;
        // multiprocessor initialization
        nThreads =gNumProcessors;
        if (inPutDimensionSizes [LAYERS] < nThreads) nThreads = inPutDimensionSizes [LAYERS];
        // make an array of parameter structures
        paramArrayPtr = (KalmanThreadParamsPtr)WMNewPtr (nThreads * sizeof(KalmanThreadParams));
        if (paramArrayPtr == NULL) throw result = MEMFAIL;
        // make an array of pthread_t
        threadsPtr =(pthread_t*)WMNewPtr(nThreads * sizeof(pthread_t));
        if (threadsPtr == NULL) throw result = MEMFAIL;
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
        paramArrayPtr[iThread].inPutWaveType = inPutWaveType;
        paramArrayPtr[iThread].inPutDataStartPtr = inPutDataStartPtr;
        paramArrayPtr[iThread].outPutDataStartPtr = outPutDataStartPtr;
        paramArrayPtr[iThread].startLayer = startLayer;
        paramArrayPtr[iThread].outPutLayer = outPutLayer;
        paramArrayPtr[iThread].xSize = inPutDimensionSizes [0];
        paramArrayPtr[iThread].ySize = inPutDimensionSizes [1];
        paramArrayPtr[iThread].zSize = layersToDo;
        paramArrayPtr[iThread].multiplier =multiplier;
        paramArrayPtr[iThread].ti=iThread; // number of this thread, starting from 0
        paramArrayPtr[iThread].tN =nThreads; // total number of threads
    }
    // create the threads
    for (iThread = 0; iThread < nThreads; iThread++){
        pthread_create (&threadsPtr[iThread], NULL, KalmanThread, (void *) &paramArrayPtr[iThread]);
    }
    // Wait till all the threads are finished
    for (iThread = 0; iThread < nThreads; iThread++){
        pthread_join (threadsPtr[iThread], NULL);
    }
    WMDisposePtr ((Ptr)threadsPtr);         // free memory for pThreads Array
    WMDisposePtr ((Ptr)paramArrayPtr);      // Free paramaterArray memory
    WaveHandleModified(outPutWaveH);        // Inform Igor that we have changed the output wave.
    p -> result = (0);
    return (0);
}


/* KalmanWaveToFrame XOP entry function
 Collapses a 3D input wave into a single 2D frame. You can get the same result with KalmanAllFrames by
 using "" as outPut String
 KalmanWaveToFrameParams
 inPutWaveH     handle to input wave
 multiplier     Multiplier, as for  for 16 bit waves containing less than 16 bits of data
 result         0 or error code
 Last Modified: 2025/06/18 by Jamie Boyd */
int KalmanWaveToFrame (KalmanWaveToFrameParamsPtr p) {
	int result=0;	// The error returned from various Wavemetrics functions
	waveHndl inPutWaveH = NIL;		// handle to the input wave
    int inPutWaveType; //  Wavetypes numeric codes for things like 32 bit floating point, 16 bit int, etc
    int inPutDimensions;	// number of dimensions in input and output waves
	CountInt inPutDimensionSizes[MAX_DIMENSIONS+1];	// an array used to hold the width, height, layers, and chunk sizes
	CountInt inPutOffset;	//offset in bytes from begnning of handle to a wave to the actual data - size of headers, units, etc.
	char *inPutDataStartPtr;
	float multiplier = (float)(p->multiplier);
	CountInt zSize;
    // threading
	UInt8 iThread, nThreads;
    KalmanThreadParamsPtr paramArrayPtr = NULL;
    pthread_t* threadsPtr = NULL;
	try {
		// Get handle to input wave. Make sure input wave exists.
		inPutWaveH = p->inPutWaveH;
		if(inPutWaveH == NIL)throw result = NON_EXISTENT_WAVE;
		// get wave data type and check that we don't have a text wave
		inPutWaveType = WaveType(inPutWaveH);
		if (inPutWaveType==TEXT_WAVE_TYPE) throw result = NOTEXTWAVES;
		//Get number of used dimensions in input wave.
		if (MDGetWaveDimensions(inPutWaveH, &inPutDimensions, inPutDimensionSizes))throw result = WAVEERROR_NOS;
		// Check that input wave is 3D
		if (inPutDimensions != 3) throw result = INPUTNEEDS_3D_WAVE;
		// Save z Size as we will be resizing dimensions array
		zSize = inPutDimensionSizes [2];
        //Get data offset for the wave
		if (MDAccessNumericWaveData(inPutWaveH, kMDWaveAccessMode0, &inPutOffset)) throw result = WAVEERROR_NOS;
		inPutDataStartPtr = (char*)(*inPutWaveH) + inPutOffset;
        // for multiprocessor initialization
        nThreads = gNumProcessors;
        // make an array of parameter structures
        paramArrayPtr = (KalmanThreadParamsPtr)WMNewPtr (nThreads * sizeof(KalmanThreadParams));
        if (paramArrayPtr == NULL) throw result = MEMFAIL;
        // make an array of pthread_t
        threadsPtr = (pthread_t*)WMNewPtr(nThreads * sizeof(pthread_t));
        if (threadsPtr == NULL) throw result = MEMFAIL;
	}catch (int result){
        if (paramArrayPtr != NULL)  WMDisposePtr ((Ptr)paramArrayPtr);
        if (threadsPtr != NULL) WMDisposePtr ((Ptr)threadsPtr);
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
        paramArrayPtr[iThread].outPutDataStartPtr = inPutDataStartPtr;
        paramArrayPtr[iThread].startLayer = 0;
        paramArrayPtr[iThread].outPutLayer =0;
        paramArrayPtr[iThread].xSize = inPutDimensionSizes [0];
        paramArrayPtr[iThread].ySize = inPutDimensionSizes [1];
        paramArrayPtr[iThread].zSize =zSize;
        paramArrayPtr[iThread].multiplier =multiplier;
        paramArrayPtr[iThread].ti=iThread; // number of this thread, starting from 0
        paramArrayPtr[iThread].tN =nThreads; // total number of threads
    }
    // create the threads
    for (iThread = 0; iThread < nThreads; iThread++){
        pthread_create (&threadsPtr[iThread], NULL, KalmanThread, (void *) &paramArrayPtr[iThread]);
    }
    // Wait till all the threads are finished
    for (iThread = 0; iThread < nThreads; iThread++){
        pthread_join (threadsPtr[iThread], NULL);
    }
    WMDisposePtr ((Ptr)threadsPtr);     // free memory for pThreads Array
    WMDisposePtr ((Ptr)paramArrayPtr);  // Free paramaterArray memory
	// Redimension wave
	inPutDimensionSizes [0] = -1;
	inPutDimensionSizes [1] = -1;
	inPutDimensionSizes [2] = 0;
	inPutDimensionSizes [3] = 0;
	MDChangeWave (inPutWaveH, -1, inPutDimensionSizes);
	WaveHandleModified(inPutWaveH);     // Inform Igor that we have changed the input wave.
	p -> result = (0);
	return (0);
}


/* Template for handling all data types for KalmanList function
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
		if (*srcWaveStarts == destWaveStart){	//this happens when collapsing a wave into the first frame
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
    CountInt nPnts;
    float multiplier;
    UInt8 ti; // number of this thread, starting from 0
    UInt8 tN; // total number of threads
} KalmanListThreadParams, *KalmanListThreadParamsPtr;


/* Each thread to average a list of waves starts with this function
Last Modified 2013/07/16 by Jamie Boyd */
void* KalmanListThread (void* threadarg){
	struct KalmanListThreadParams* p;
	p = (struct KalmanListThreadParams*) threadarg;
	CountInt nPnts= p->nPnts;
	float multiplier = p->multiplier;
	UInt8 ti = p->ti;
	UInt8 tN = p->tN;
	CountInt pntsPerThread = nPnts/tN;
	CountInt startPos = ti * pntsPerThread; // which point to start this thread on depends on thread number * points per thread. ti is 0 based
	if (ti == (tN - 1)) pntsPerThread += (nPnts % tN); // last thread gets any extra points
	switch (p->inPutWaveType) {
	case NT_I8:
		KalmanListT ((char**)p->inPutDataStartsPtr, (char*)p->outPutDataStartPtr, p->nWaves, startPos, startPos + pntsPerThread, multiplier);
		break;
	case (NT_I8 | NT_UNSIGNED):
		KalmanListT ((unsigned char**)p->inPutDataStartsPtr, (unsigned char*)p->outPutDataStartPtr, p->nWaves, startPos, startPos + pntsPerThread, multiplier);				break;
	case NT_I16:
		KalmanListT ((short**)p->inPutDataStartsPtr, (short*)p->outPutDataStartPtr, p->nWaves, startPos, startPos + pntsPerThread, multiplier);
		break;
	case (NT_I16 | NT_UNSIGNED):
		KalmanListT ((unsigned short**)p->inPutDataStartsPtr, (unsigned short*)p->outPutDataStartPtr, p->nWaves, startPos, startPos + pntsPerThread, multiplier);
		break;
	case NT_I32:
		KalmanListT ((long**)p->inPutDataStartsPtr, (long*)p->outPutDataStartPtr, p->nWaves, startPos, startPos + pntsPerThread, multiplier);
		break;
	case (NT_I32| NT_UNSIGNED):
		KalmanListT ((unsigned long**)p->inPutDataStartsPtr, (unsigned long*)p->outPutDataStartPtr, p->nWaves, startPos, startPos + pntsPerThread, multiplier);
		break;
	case NT_FP32:
		KalmanListT ((float**)p->inPutDataStartsPtr, (float*)p->outPutDataStartPtr, p->nWaves, startPos, startPos + pntsPerThread, multiplier);
		break;
	case NT_FP64:
		KalmanListT ((double**)p->inPutDataStartsPtr, (double*)p->outPutDataStartPtr, p->nWaves, startPos, startPos + pntsPerThread, multiplier);
		break;
	default:	// Unknown data type - possible in a future version of Igor.
		throw NT_FNOT_AVAIL;
		break;
	}
	return 0;
}


/* KalmanList XOP entry function
 Averages a semicolon-separated list of waves. Each wave must have same data type and same dimensions. This is not used by the
 twoPhoton acquisition code so we print some more information in the error cases with XOPNotice
 KalmanListParams
 inPutList          semicolon separated list of input waves, with paths
 outPutPath         path and wavename of output wave
 multiplier         Multiplier for 16 bit waves containing less than 16 bits of data
 overwrite          0 to give errors when wave already exists. non-zero to overwrite existing wave.
 result             0 for success, else error code
 Last Modified 2014/02/13 by Jamie Boyd */
extern "C" int KalmanList (KalmanListParamsPtr p) {
	int result = 0;	// The error returned from various Wavemetrics functions
	waveHndl outPutWaveH = NULL; // handle to output wave
	waveHndl* handleList; // pointer to an array of handles for input waves
	DFPATH inPutPath, outPutPath ;	// string to hold data folder path of input wave
	WVNAME inPutWaveName, outPutWaveName;	// C strings to hold names of input and output waves
	DataFolderHandle inPutDFHandle, outPutDFHandle;	// Handles to datafolders of input and output waves
	int inPutWaveType; //  Wavetypes numeric codes for things like 32 bit floating point, 16 bit int, etc
    int inPutDimensions;	// number of numDimensions in input and output waves
	CountInt inPutDimensionSizes[MAX_DIMENSIONS+1];	// an array used to hold the width, height, layers, and chunk sizes
	char* outPutDataStartPtr; // Pointer to start of output wave
	Ptr* inPutDataStartsPtr = nullptr; // Pointer to an array of pointers for starts of input data
	UInt8 overWrite = p->overwrite; // if it is O.K. to overwrite an existing wave
	UInt8 isOverWriting = 0; // 0 if using a separate output wave, 1 for overwriting first wave in list with results
	float multiplier = p->multiplier;
	UInt16 numWaves;	//number of input waves in the input list
	CountInt waveOffset;	//offset in bytes from begnning of handle to a wave to the actual data - size of headers, units, etc.
	CountInt nPnts;
    UInt8 iThread, nThreads;
    KalmanListThreadParamsPtr paramArrayPtr = nullptr;
    pthread_t* threadsPtr = nullptr;
	try {
		// Check that input string exists
		if (WMGetHandleSize (p->inPutList) == 0) throw result = NON_EXISTENT_WAVE;
		// If outPutPath is empty string, we are overwriting first wave in list with results
		if (WMGetHandleSize (p->outPutPath) == 0){
			if (overWrite == NO_OVERWITE) throw result = OVERWRITEALERT;
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
		int tInPutDimensions;	// temp number of numDimensions for each wave in array
		CountInt tInPutDimensionSizes[MAX_DIMENSIONS+1];	// temp width, height, layers, and chunk sizes for each wave in array
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
		if (MDGetWaveDimensions(handleList[0], &inPutDimensions, inPutDimensionSizes))throw result = WAVEERROR_NOS;
		nPnts = inPutDimensionSizes [0];
		for (int id =1; id < inPutDimensions; id +=1){
			nPnts *= inPutDimensionSizes [id];
		}
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
			if (handleList [iw] == NULL){
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
			if (MDGetWaveDimensions(handleList[iw], &tInPutDimensions, tInPutDimensionSizes))throw result = WAVEERROR_NOS;
			if (tInPutDimensions != inPutDimensions){
#ifdef IPARSEWAVE_INF0
				sprintf(XOPbuffer, "The number of dimensions of wave %d in the input list did not match the number of dimenisons of the first wave in the list.\r", iw);
				XOPNotice (XOPbuffer);
#endif
				throw result = NOTSAMEDIMSIZE;
			}
			// check sizes of each dimension
			for (id=0; id < MAX_DIMENSIONS; id +=1){
				if (tInPutDimensionSizes [id] != inPutDimensionSizes [id]){
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
			if ( MDMakeWave (&outPutWaveH, outPutWaveName, outPutDFHandle, inPutDimensionSizes, inPutWaveType, overWrite)) throw result = WAVEERROR_NOS;
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
        nThreads = gNumProcessors;
        paramArrayPtr = (KalmanListThreadParamsPtr)WMNewPtr(nThreads * sizeof(KalmanListThreadParams));
        if (paramArrayPtr == nullptr) throw result = MEMFAIL;
        // make an array of pthread_t
        threadsPtr =(pthread_t*)WMNewPtr(nThreads * sizeof(pthread_t));
        if (threadsPtr == nullptr) throw result = MEMFAIL;
	}catch (int result){
        if (threadsPtr != nullptr) WMDisposePtr((Ptr)threadsPtr);
        if (paramArrayPtr != nullptr) WMDisposePtr ((Ptr)paramArrayPtr);
		if (inPutDataStartsPtr != nullptr) WMDisposePtr ((Ptr)inPutDataStartsPtr);
        WMDisposeHandle(p->inPutList);
        WMDisposeHandle(p->outPutPath);
		// set result
		p -> result = (result - FIRST_XOP_ERR);	// XFUNC error code
		return (result);
	}
    // fill threadArray
	for (iThread = 0; iThread < nThreads; iThread++){
		paramArrayPtr[iThread].inPutWaveType = inPutWaveType;
        paramArrayPtr[iThread].inPutDataStartsPtr = inPutDataStartsPtr;
		paramArrayPtr[iThread].outPutDataStartPtr = outPutDataStartPtr;
		paramArrayPtr[iThread].nWaves=numWaves;
		paramArrayPtr[iThread].nPnts = nPnts;
		paramArrayPtr[iThread].multiplier = multiplier;
		paramArrayPtr[iThread].ti=iThread; // number of this thread, starting from 0
		paramArrayPtr[iThread].tN =nThreads; // total number of threads
	}
    
	// create the threads
	for (iThread = 0; iThread < nThreads; iThread++){
		pthread_create (&threadsPtr[iThread], NULL, KalmanListThread, (void *) &paramArrayPtr[iThread]);
	}
	// Wait till all the threads are finished
	for (iThread = 0; iThread < nThreads; iThread++){
		pthread_join (threadsPtr[iThread], NULL);
	}
	WMDisposePtr ((Ptr)threadsPtr);         // free memory for pThreads Array
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


/* Template for doing sequential Kalman averaging, src wave is the new wave,
 dest wave is the old wave already averaged iKal times
 Last Modified 2014/01/28 by Jamie Boyd */
template <typename T> void KalmanNextT (T* srcWaveStart, T* destWaveStart, CountInt nPoints, UInt16 iKal) {
	T* srcWave;
	T* destWave;
	T* destWaveEnd = destWaveStart + nPoints;
	if (iKal == 0){
		for (srcWave = srcWaveStart, destWave = destWaveStart; destWave < destWaveEnd; srcWave++, destWave++)
			*destWave = *srcWave;
	}else{
		for (srcWave = srcWaveStart, destWave = destWaveStart; destWave < destWaveEnd; srcWave++, destWave++)
			*destWave= (*destWave  * iKal +  *srcWave)/(iKal + 1);
	}
}


/* Structure to pass data to each KalmanNext Thread
 Last Modified: 2014/01/28 by Jamie Boyd */
typedef struct KalmanNextThreadParams{
    int inPutWaveType;
    char* inPutDataStartPtr;
    char* outPutDataStartPtr;
    CountInt nPnts;
    UInt16 iKal;
    UInt8 ti; // number of this thread, starting from 0
    UInt8 tN; // total number of threads
} KalmanNextThreadParams, *KalmanNextThreadParamsPtr;


/* Each thread to do sequential Kalmaning starts with this function
 Last Modified 2014/01/29 by Jamie Boyd */
void* KalmanNextThread (void* threadarg){
	struct KalmanNextThreadParams* p;
	p = (struct KalmanNextThreadParams*) threadarg;
	CountInt nPnts= p->nPnts;
	UInt8 ti = p->ti;
	UInt8 tN = p->tN;
	CountInt pntsPerThread = nPnts/tN;
	CountInt startPos = ti * pntsPerThread; // which point to start this thread on depends on thread number * points per thread. ti is 0 based
	if (ti == (tN - 1)) pntsPerThread += (nPnts % tN); // last thread gets any extra points
	switch (p->inPutWaveType) {
        case NT_I8:
            KalmanNextT ((char*)p->inPutDataStartPtr + startPos, (char*)p->outPutDataStartPtr + startPos, pntsPerThread, p->iKal);
            break;
        case (NT_I8 | NT_UNSIGNED):
            KalmanNextT ((unsigned char*)p->inPutDataStartPtr + startPos, (unsigned char*)p->outPutDataStartPtr + startPos, pntsPerThread, p->iKal);
            break;
        case NT_I16:
            KalmanNextT ((short*)p->inPutDataStartPtr + startPos, (short*)p->outPutDataStartPtr+ startPos, pntsPerThread, p->iKal);
            break;
        case (NT_I16 | NT_UNSIGNED):
            KalmanNextT ((unsigned short*)p->inPutDataStartPtr + startPos, (unsigned short*)p->outPutDataStartPtr + startPos, pntsPerThread, p->iKal);
            break;
        case NT_I32:
            KalmanNextT ((long*)p->inPutDataStartPtr + startPos, (long*)p->outPutDataStartPtr + startPos, pntsPerThread, p->iKal);
            break;
        case (NT_I32| NT_UNSIGNED):
            KalmanNextT ((unsigned long*)p->inPutDataStartPtr + startPos, (unsigned long*)p->outPutDataStartPtr + startPos, pntsPerThread, p->iKal);
            break;
        case NT_FP32:
            KalmanNextT ((float*)p->inPutDataStartPtr + startPos, (float*)p->outPutDataStartPtr + startPos, pntsPerThread, p->iKal);
            break;
        case NT_FP64:
            KalmanNextT ((double*)p->inPutDataStartPtr + startPos, (double*)p->outPutDataStartPtr + startPos, pntsPerThread, p->iKal);
            break;
	}
	return nullptr;
}


/* KalmanNext XOP entry function
 Averages a wave into an already averaged wave. Both waves must have same data type and same dimensions.
 KalmanNextParams:
 inPutWaveH     handle to input wave
 outPutWaveH    handle to output wave
 iKal           index of wave are we adding
 result         0 or error code
 Last Modified 2014/01/28 by Jamie Boyd */
extern "C" int KalmanNext (KalmanNextParamsPtr p) {
	int result = 0;	// The error returned from various Wavemetrics functions
    waveHndl outPutWaveH = NULL; // handle to output wave
	waveHndl inPutWaveH; // handle to input wave
    int inPutWaveType, outPutWaveType; //  Wavetypes numeric codes for things like 32 bit floating point, 16 bit int, etc
	int inPutDimensions, outPutDimensions;	// number of numDimensions in input and output waves
	CountInt inPutDimensionSizes[MAX_DIMENSIONS+1];	// an array used to hold the width, height, layers, and chunk sizes
	CountInt outPutDimensionSizes[MAX_DIMENSIONS+1];
    char* outPutDataStartPtr; // Pointer to start of output wave
	char* inPutDataStartPtr;
	CountInt inPutOffset, outPutOffset;	//offset in bytes from begnning of handle to a wave to the actual data - size of headers, units, etc.
	CountInt nPnts = 1;
	UInt16 iKal;
    UInt8 iThread, nThreads;
    KalmanNextThreadParamsPtr paramArrayPtr = NULL;
    pthread_t* threadsPtr = NULL;
	try {
		// get iKal
        iKal = (UInt16)p->iKal;
        // Get handle to input wave. Make sure input wave exists.
		inPutWaveH = p->inPutWaveH;
		if(inPutWaveH == NIL)throw result = NON_EXISTENT_WAVE;
		// get wave data type and check that we don't have a text wave
		inPutWaveType = WaveType(inPutWaveH);
		if (inPutWaveType==TEXT_WAVE_TYPE) throw result = NOTEXTWAVES;
		//Get number of used dimensions in input wave.
		if (MDGetWaveDimensions(inPutWaveH, &inPutDimensions, inPutDimensionSizes))throw result = WAVEERROR_NOS;
        // Get handle to outPut wave. Make sure outPut wave exists.
		outPutWaveH = p->outPutWaveH;
		if(outPutWaveH == NIL)throw result = NON_EXISTENT_WAVE;
		// get wave data type and check that we don't have a text wave
		outPutWaveType = WaveType(outPutWaveH);
		if (outPutWaveType==TEXT_WAVE_TYPE) throw result = NOTEXTWAVES;
		//Get number of used dimensions in outPut wave.
		if (MDGetWaveDimensions(outPutWaveH, &outPutDimensions, outPutDimensionSizes))throw result = WAVEERROR_NOS;
        // check that waves are the same types and dimension sizes
        if (inPutWaveType != outPutWaveType) throw result = NOTSAMEDIMSIZE;
        if (inPutDimensions != outPutDimensions) throw result = NOTSAMEDIMSIZE;
        // check sizes of each dimension, and calculate total number of points as well
        UInt8 id;
        for (id=0; id < inPutDimensions; id +=1){
            if (inPutDimensionSizes [id] != outPutDimensionSizes [id]) throw result = NOTSAMEDIMSIZE;
            nPnts *= inPutDimensionSizes [id];
        }
		//Get data offset for the waves
		if (MDAccessNumericWaveData(inPutWaveH, kMDWaveAccessMode0, &inPutOffset)) throw result = WAVEERROR_NOS;
		inPutDataStartPtr = (char*)(*inPutWaveH) + inPutOffset;
        if (MDAccessNumericWaveData(outPutWaveH, kMDWaveAccessMode0, &outPutOffset)) throw result = WAVEERROR_NOS;
		outPutDataStartPtr = (char*)(*outPutWaveH) + outPutOffset;
        // multiprocessor initialization
        nThreads =gNumProcessors;
        // make an array of parameter structures
        paramArrayPtr = (KalmanNextThreadParamsPtr)WMNewPtr (nThreads * sizeof(KalmanNextThreadParams));
        // make an array of pthread_t
        threadsPtr =(pthread_t*)WMNewPtr(nThreads * sizeof(pthread_t));
    }catch (int result){
        if (threadsPtr != NULL) WMDisposePtr ((Ptr)threadsPtr);
        if (paramArrayPtr != NULL) WMDisposePtr ((Ptr)paramArrayPtr);
        p -> result = (double)(result - FIRST_XOP_ERR);
#ifdef NO_IGOR_ERR
            return (0);
 #else
            return (result);
#endif
    }
    // make an array of parameter structures
    for (int iThread = 0; iThread < nThreads; iThread++){
        paramArrayPtr[iThread].inPutWaveType = inPutWaveType;
        paramArrayPtr[iThread].inPutDataStartPtr = inPutDataStartPtr;
        paramArrayPtr[iThread].outPutDataStartPtr = outPutDataStartPtr;
        paramArrayPtr[iThread].nPnts =nPnts;
        paramArrayPtr[iThread].iKal =iKal;
        paramArrayPtr[iThread].ti=iThread; // number of this thread, starting from 0
        paramArrayPtr[iThread].tN =nThreads; // total number of threads
    }
    // create the threads
    for (iThread = 0; iThread < nThreads; iThread++){
        pthread_create (&threadsPtr[iThread], NULL, KalmanNextThread, (void *) &paramArrayPtr[iThread]);
    }
    // Wait till all the threads are finished
    for (iThread = 0; iThread < nThreads; iThread++){
        pthread_join (threadsPtr[iThread], NULL);
    }
    // free memory for pThreads Array
    WMDisposePtr ((Ptr)threadsPtr);
    // Free paramaterArray memory
    WMDisposePtr ((Ptr)paramArrayPtr);
	MDChangeWave (outPutWaveH, -1, outPutDimensionSizes);
    // Inform Igor that we have changed the input wave.
	WaveHandleModified(outPutWaveH);
	p -> result = (0);
	return (0);
}
