/*
    twoP.h -- equates for twoP XOP
*/
#ifndef TWOPHOTON_H_
#define TWOPHOTON_H_

#include "ParseWavePath.h"              // Utility to parse strings into data folder paths and wave names
#include "XOPResources.h"                // Contains definition of XOP_TOOLKIT_VERSION
#include "XOPStandardHeaders.h"            // Include ANSI headers, Mac headers, IgorXOP.h, XOP.h and XOPSupport.h


//#define NO_IGOR_ERR   // when defined, all functions return 0 to avoid modal dialogs
#undef NO_IGOR_ERR      // when not defined, functions return error codes that invoke modal dialogs

/* twoP custom error codes - also serve as return values from functions, after subtracting FIRST_XOP_ERR */
#define OLD_IGOR                1 + FIRST_XOP_ERR
#define NON_EXISTENT_WAVE       2 + FIRST_XOP_ERR
#define INPUTNEEDS_3D_WAVE      3 + FIRST_XOP_ERR
#define OUTPUTNEEDS_2D3D_WAVE   4 + FIRST_XOP_ERR
#define NOTSAMEWAVETYPE         5 + FIRST_XOP_ERR
#define NOTSAMEDIMSIZE          6 + FIRST_XOP_ERR
#define INVALIDOUTPUTFRAME      7 + FIRST_XOP_ERR
#define INVALIDINPUTFRAME       8 + FIRST_XOP_ERR
#define OUTPUTNEEDS_2D_WAVE     9 + FIRST_XOP_ERR
#define BADKERNEL               10 + FIRST_XOP_ERR
#define INPUTNEEDS_2D3D_WAVE    11 + FIRST_XOP_ERR
#define NO_INPUT_STRING         12 + FIRST_XOP_ERR
#define BADFACTOR               13 + FIRST_XOP_ERR
#define BADDSTYPE               14 + FIRST_XOP_ERR
#define WAVEERROR_NOS           15 + FIRST_XOP_ERR
#define OVERWRITEALERT          16 + FIRST_XOP_ERR
#define NOTEXTWAVES             17 + FIRST_XOP_ERR
#define BADDIMENSION            18 + FIRST_XOP_ERR
#define NOT16OR32               19 + FIRST_XOP_ERR
#define OUTPUTNEEDS_3D_WAVE     20 + FIRST_XOP_ERR
#define BADWAVEINLIST           21 + FIRST_XOP_ERR
#define BADSYMKERNEL            22 + FIRST_XOP_ERR
#define NOTUNSIGNED             23 + FIRST_XOP_ERR
#define MEMFAIL                 24 + FIRST_XOP_ERR
#define NUMTYPE                 25 + FIRST_XOP_ERR
#define INPUT_RANGE             26 + FIRST_XOP_ERR

// mnemonic defines
#define OVERWRITE 1
#define NO_OVERWITE 0
#define KILLOUTPUT 1
#define NO_KILLOUTPUT 0

//preprocessor macro to swap 2 values
#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp

// Threading platform-dependent globals,includes, and macros
extern UInt8 gNumProcessors;

// include pThreads library on Windows
#ifdef  _WINDOWS_
// include pThreads library on Windows
#include "pthread.h"
#include "sched.h"
#include "semaphore.h"
static inline int num_processors() {
    SYSTEM_INFO info;
    GetSystemInfo(&info);
    return info.dwNumberOfProcessors;
}
#endif

// include native pThreads library on MacOS
#ifdef __GNUC__
#include <pthread.h>
#include <sys/types.h>
#include <sys/sysctl.h>
static inline int num_processors() {
    int np = 1;
    size_t length = sizeof(np);
    sysctlbyname("hw.ncpu", &np, &length, NULL, 0);
    return np;
}
#endif

// Structure definitions. All structures passed to Igor are two-byte aligned
#pragma pack(2)

//Kalman Averaging
typedef struct KalmanAllFramesParams {
    double overWrite;    //0 to give errors when wave already exists. non-zero to overwrite existing wave without warning.
    double multiplier;    // Multiplier for,e.g., 16 bit waves containing less than 16 bits of data
    Handle outPutPath;    // A handle to a string containing path to output wave we want to make
    waveHndl inPutWaveH;    // handle to a 3D input wave
    UserFunctionThreadInfoPtr tp; // Pointer to Igor private data.
    double result;
}KalmanAllFramesParams, *KalmanAllFramesParamsPtr;

typedef struct KalmanSpecFramesParams {
    double multiplier;    // Multiplier for 16 bit waves containing less than 16 bits of data
    double outPutLayer;    // layer of output wave to modify
    waveHndl outPutWaveH;//handle to output wave
    double endLayer;    // end of lyaers to average
    double startLayer;    // start of layers to average for input wave
    waveHndl inPutWaveH;// handle to input wave
    UserFunctionThreadInfoPtr tp; // Pointer to Igor private data.
    double result;
}KalmanSpecFramesParams, * KalmanSpecFramesParamsPtr;

typedef struct KalmanWaveToFrameParams {
    double multiplier;    // Multiplier for 16 bit waves containing less than 16 bits of data
    waveHndl inPutWaveH;// handle to input wave
    UserFunctionThreadInfoPtr tp; // Pointer to Igor private data.
    double result;
} KalmanWaveToFrameParams, * KalmanWaveToFrameParamsPtr;

typedef struct KalmanListParams {
    double overwrite;    //0 to give errors when wave already exists. non-zero to overwrite existing wave.
    double multiplier; // Multiplier for 16 bit waves containing less than 16 bits of data
    Handle outPutPath;    // path and wavename of output wave
    Handle inPutList;    //semicolon separated list of input waves, with paths
    UserFunctionThreadInfoPtr tp; // Pointer to Igor private data.
    double result;
}KalmanListParams, * KalmanListParamsPtr;

typedef struct KalmanNextParams {
    double iKal; //which number of wave are we adding
    waveHndl outPutWaveH;//handle to output wave
    waveHndl inPutWaveH; // handle to iinput wave
    UserFunctionThreadInfoPtr tp; // Pointer to Igor private data.
    double result;
}KalmanNextParams, * KalmanNextParamsPtr;

// Project Image
typedef struct ProjectAllFramesParams {
    double projMode;        // variable for kind of projection, 0 is minimum intensity, 1 is maximum intensity
    double overwrite;    //0 to give errors when wave already exists. non-zero to cheerfully overwrite existing wave.
    double flatDimension;    //Which dimension we want to collapse on, 0 for x, 1 for y, 2 for z
    Handle outPutPath;    // A handle to a string containing path to output wave we want to make
    waveHndl inPutWaveH; //handle to the input wave
    UserFunctionThreadInfoPtr tp; // Pointer to Igor private data.
    double result;
} ProjectAllFramesParams, * ProjectAllFramesParamsPtr;

typedef struct ProjectSpecFramesParams {
    double projMode;             // 0 if minimum intensity projection, 1 for maximum intensity projection
    double flatDimension;    // the dimension that we are projecting along
    double outPutLayer;        //the layer in the output wave that receives the projection
    waveHndl outPutWaveH;    // A handle to output wave
    double inPutEndLayer;    //end of range of layers to project
    double inPutStartLayer;    //start of range of layers to project
    waveHndl inPutWaveH;    //handle to the input wave
    UserFunctionThreadInfoPtr tp; // Pointer to Igor private data.
    double result;
} ProjectSpecFramesParams, * ProjectSpecFramesParamsPtr;

typedef struct ProjectSliceParams {
    double slice;    // X, Y  or Z slice to get
    waveHndl outPutWaveH;    //handle to the output wave
    waveHndl inPutWaveH;//handle to the input wave
    UserFunctionThreadInfoPtr tp; // Pointer to Igor private data.
    double result;
}ProjectSliceParams, * ProjectSliceParamsPtr;

// LSM Utilities
typedef struct GetSetNumProcessorsParams{
    double result;
}GetSetNumProcessorsParams, *GetSetNumProcessorsParamsPtr;

typedef struct SwapEvenParams {
    waveHndl w1;
    UserFunctionThreadInfoPtr tp; // Pointer to Igor private data.
    double result;
}SwapEvenParams, * SwapEvenParamsPtr;

typedef struct DownSampleParams {
    double dsType;
    double boxFactor;
    waveHndl w1;
    UserFunctionThreadInfoPtr tp; // Pointer to Igor private data.
    double result;
} DownSampleParams, * DownSampleParamsPtr;

typedef struct DecumulateParams {
    double bitSize;   //bitsize of the counter. either 24 or 32
    waveHndl w1;
    UserFunctionThreadInfoPtr tp; // Pointer to Igor private data.
    double result;
} DecumulateParams, * DecumulateParamsPtr;

typedef struct TransposeFramesParams {
    waveHndl w1;
    UserFunctionThreadInfoPtr tp; // Pointer to Igor private data.
    double result;
}TransposeFramesParams, * TransposeFramesParamsPtr;

// Filter
typedef struct ConvolveFramesParams {
    double overWrite; // 1 if it is o.k. to overwrite existing waves, 0 to exit with error if overwriting will occur
    waveHndl kernelH; // convolution kernel, a 2D wave odd number of pixels high and wide, or 1D odd number wave for Sym
    double outPutType; // 0 for same type as input wave, non-zero for floating point wave
    Handle outPutPath;    // A handle to a string containing path to output wave we want to make, or empty string to overwrite existing wave
    waveHndl inPutWaveH; //input wave. needs to be 2D or 3D wave
    UserFunctionThreadInfoPtr tp; // Pointer to Igor private data.
    double result;
} ConvolveFramesParams, * ConvolveFramesParamsPtr;

typedef struct MedianFramesParams {
    double overWrite; // 1 if it is o.k. to overwrite existing waves, 0 to exit with error
    double kWidth; //width of the area over which to calculate the median. Must be an odd number
    Handle outPutPath;    // A handle to a string containing path to output wave we want to make, or empty string to overwrite existing wave
    waveHndl inPutWaveH;//input wave. needs to be 2D or 3D wave
    UserFunctionThreadInfoPtr tp; // Pointer to Igor private data.
    double result;
} MedianFramesParams, * MedianFramesParamsPtr;

// Return to default structure packing
#pragma pack()

/* Prototypes */
HOST_IMPORT int XOPMain(IORecHandle ioRecHandle);
//LSM Utilities
extern "C" int GetSetNumProcessors(GetSetNumProcessorsParamsPtr p);
extern "C" int SwapEven(SwapEvenParamsPtr);
extern "C" int DownSample(DownSampleParamsPtr p);
extern "C" int Decumulate(DecumulateParamsPtr p);
extern "C" int TransposeFrames(TransposeFramesParamsPtr p);
// Kalman Averaging
extern "C" int KalmanAllFrames(KalmanAllFramesParamsPtr);
extern "C" int KalmanSpecFrames(KalmanSpecFramesParamsPtr);
extern "C" int KalmanWaveToFrame(KalmanWaveToFrameParamsPtr);
extern "C" int KalmanList(KalmanListParamsPtr p);
extern "C" int KalmanNext(KalmanNextParamsPtr p);
// Project Image
extern "C" int  ProjectAllFrames(ProjectAllFramesParamsPtr p);
extern "C" int  ProjectSpecFrames(ProjectSpecFramesParamsPtr p);
extern "C" int  ProjectXSlice(ProjectSliceParamsPtr p);
extern "C" int  ProjectYSlice(ProjectSliceParamsPtr p);
extern "C" int  ProjectZSlice(ProjectSliceParamsPtr p);
//Filter frames
extern "C" int  ConvolveFrames(ConvolveFramesParamsPtr p);
extern "C" int  SymConvolveFrames(ConvolveFramesParamsPtr p);
extern "C" int  MedianFrames(MedianFramesParamsPtr p);
template <typename T> T medianT(UInt32 n, T* dataStrtPtr);
#endif
