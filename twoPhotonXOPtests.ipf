#pragma TextEncoding = "UTF-8"
#pragma DefaultTab={3,20,4}		// Set default tab width in Igor Pro 9 and later
#pragma rtGlobals=3		// Use modern global access method and strict wave access.


Menu "Macros"
	"test twoPhotonXOP", test_twoPhotonXOP ()
End


function closeTwoPGraphs()
	string graphs =  WinList("twoPxop_*", ";", "WIN:1" )
	variable iGraph, nGraphs = itemsinlist (graphs, ";")
	for (iGraph = 0; iGraph < nGraphs; iGraph += 1)
		doWindow/k $StringFromList(iGraph, graphs, ";")
	endfor
end


function/wave makeKernel (width)
	variable width
	
	if (mod (width, 2) ==0)
		width +=1
	endif
	make/FREE/o/s/n= (width, width) gWave
	setscale/I x -1,1, "m" gWave
	setscale/I y -1, 1, "m" gWave
	gWave = Gauss (x, 0, 0.5, y, 0, 0.5)
	variable sumVal =  sum (gwave)
	gwave /= sumVal
	return gwave
end


function/wave makeSymkernel (width)
	variable width
	
	if (mod (width, 2) ==0)
		width +=1
	endif
	make/FREE/o/s/n= (width) gWave
	setscale/I x -1,1, "m" gWave
	gWave = Gauss (x, 0, 0.5)
	variable sumVal =  sum (gwave)
	gwave /= sumVal
	return gwave
end

function test_twoPhotonXOP ()
	// make sure XOP is loaded
	string IgorKind= stringbykey ("IGORKIND", igorinfo (0),":", ";")
	string XOPlist = igorinfo (10)
	if (cmpstr (IgorKind, "pro64") == 0)
		if (WhichListItem("twoPhoton64", XOPlist, ";") == -1)
			print "TwoPhoton64 XOP is not loaded!"
			return 1
		endif
	else
		if (WhichListItem("twoPhotonx86", XOPlist, ";") == -1)
			print "TwoPhotonx86 XOP is not loaded!"
			return 1
		endif
		return 1
	endif
	//timer stuff
	variable testNum =0
	Variable timerRefNum, timerMicroSeconds
	make/o/d/n=100 root:testScores
	WAVE testScores = root:testScores
	setscale d 0, 0, "s" testScores 
	make/o/t/n=100 root:testType
	WAVE/T testType = root:testType
	// How many cores?
	printf "%d processor cores are available to twoPhotonXOP.\r", GetSetNumProcessors()
	// close all 2p graph windows
	closeTwoPGraphs()
	// make a large-ish stack 
	make/o/w/u/n =(1000,500,300) root:theStack
	WAVE theStack = root:theStack
	// mimic a noisy 12 bit A/D
	MultiThread /NT=(ThreadProcessorCount) theStack = 2^11 + (2^10) * (sin ((x-z)/60) * cos ((y+z)/40)) + enoise (2^10)
	
	// display stack
	execute "twoPxop_theStack()"; doupdate; sleep/S 1
	
	// Kalman All Frames
	testType [testNum]=" Kalman All Frames"
	timerRefNum = StartMSTimer
	KalmanAllFrames (theStack, "root:KalmanAllFrames_out", 16, 1)
	timerMicroSeconds = StopMSTimer(timerRefNum)
	testScores [testNum] = timerMicroSeconds/1E6
	testNum += 1
	wave KalManAllFrames_out = root:KalManAllFrames_out
	execute "twoPxop_KalmanAllFrames_out()"
	doupdate; sleep/S 1
	
	// Kalman Spec Frames
	testType [testNum]=" Kalman Specified Frames x 10"
	testScores [testNum] = 0
	make/w/u/o/n = (1000, 500, 10) root:KalmanSpecFrames_out
	execute "twoPxop_KalmanSpecFrames_out()" ;doupdate
	WAVE KalmanSpecFrames_out = root:KalmanSpecFrames_out
	variable ispec
	for (iSpec = 0; iSPec < 10; iSpec += 1)
		timerRefNum = StartMSTimer
		KalmanSpecFrames (theStack, (iSpec * 20), (20 + iSpec * 20) , KalmanSpecFrames_out, iSpec, 16)
		timerMicroSeconds = StopMSTimer(timerRefNum)
		testScores [testNum] += timerMicroSeconds/1E6
		ModifyImage/w= twoPxop_KalmanSpecFrames_out KalmanSpecFrames_out plane=iSPec; doupdate;sleep/S 0.25
	endfor
	testNum +=1
	
	// Kalman wave to Frame
	duplicate/o/r=[] [] [0,100] theStack, root:KalmanWaveToFrame_out
	wave KalmanWaveToFrame_out = root:KalmanWaveToFrame_out
	execute "twoPxop_KalmanWaveToFrame_out()"
	testType [testNum]="Kalman Wave to Frame"
	timerRefNum = StartMSTimer
	KalmanWaveToFrame (KalmanWaveToFrame_out, 16)
	timerMicroSeconds = StopMSTimer(timerRefNum)
	testScores [testNum] = timerMicroSeconds/1E6
	testNum +=1
	doupdate;sleep/S 1
	
	
	DoAlert /T="Testing twoPhotonXOP" 0, "Kalman Averaging"
	// close all 2p graph windows
	closeTwoPGraphs()
	
	// Project Frames X
	//X min
	testType [testNum]="Project All Frames dim =X, mode = Min"
	timerRefNum = StartMSTimer
	ProjectAllFrames (theStack, "root:ProjectFrames_output_x", 0, 1, 0)
	timerMicroSeconds = StopMSTimer(timerRefNum)
	testScores [testNum] = timerMicroSeconds/1E6
	testNum +=1
	execute "twoPxop_ProjFrames_out_X()"
	doupdate;sleep/S 1
	
	// x Max
	testType [testNum]="Project All Frames dim =X, mode = Max"
	timerRefNum = StartMSTimer
	ProjectAllFrames (theStack, "root:ProjectFrames_output_x", 0, 1, 1)
	timerMicroSeconds = StopMSTimer(timerRefNum)
	testScores [testNum] = timerMicroSeconds/1E6
	testNum +=1
	DoWindow/T twoPxop_ProjFrames_out_x "Project All Frames X Maximum"
	doupdate;sleep/S 1
	
	//x Average
	testType [testNum]="Project All Frames dim =X, mode = Mean"
	timerRefNum = StartMSTimer
	ProjectAllFrames (theStack, "root:ProjectFrames_output_x", 0, 1, 2)
	timerMicroSeconds = StopMSTimer(timerRefNum)
	testScores [testNum] = timerMicroSeconds/1E6
	testNum +=1
	DoWindow/T twoPxop_ProjFrames_out_x "Project All Frames X Mean"
	doupdate;sleep/S 1
	
	// X Median
	testType [testNum]="Project All Frames dim =X, mode = Median"
	timerRefNum = StartMSTimer
	ProjectAllFrames (theStack, "root:ProjectFrames_output_x", 0, 1, 3)
	timerMicroSeconds = StopMSTimer(timerRefNum)
	testScores [testNum] = timerMicroSeconds/1E6
	testNum +=1
	DoWindow/T twoPxop_ProjFrames_out_x "Project All Frames X Median"
	doupdate;sleep/S 1
	
	// Project X slice
	WAVE ProjectFrames_output = root:ProjectFrames_output_x
	ProjectFrames_output = 0
	testType [testNum]="Project X Slice x 10"
	testScores [testNum] = 0
	for (iSpec = 0; iSPec < 10; iSpec += 1)
		timerRefNum = StartMSTimer
		ProjectXSlice (theStack, ProjectFrames_output, 100 + iSpec)
		timerMicroSeconds = StopMSTimer(timerRefNum)
		testScores [testNum] += timerMicroSeconds/1E6
		DoWindow/T twoPxop_ProjFrames_out_x "Project X Slice " + num2str (100 + iSpec)
		doupdate;sleep/S 0.5
	endfor
	testNum +=1
	
	
	// Project Frames Y
	// Y min
	testType [testNum]="Project All Frames dim =Y, mode = Min"
	timerRefNum = StartMSTimer
	ProjectAllFrames (theStack, "root:ProjectFrames_output_y", 1, 1, 0)
	timerMicroSeconds = StopMSTimer(timerRefNum)
	testScores [testNum] = timerMicroSeconds/1E6
	testNum +=1
	execute "twoPxop_ProjFrames_out_y()"
	doupdate;sleep/S 1
	
	// Y Max
	testType [testNum]="Project All Frames dim =Y, mode = Max"
	timerRefNum = StartMSTimer
	ProjectAllFrames (theStack, "root:ProjectFrames_output_y", 1, 1, 1)
	timerMicroSeconds = StopMSTimer(timerRefNum)
	testScores [testNum] = timerMicroSeconds/1E6
	testNum +=1
	DoWindow/T twoPxop_ProjFrames_out_y "Project All Frames Y Maximum"
	doupdate;sleep/S 1
	
	//Y Average
	testType [testNum]="Project All Frames dim =Y, mode = Mean"
	timerRefNum = StartMSTimer
	ProjectAllFrames (theStack, "root:ProjectFrames_output_y", 1, 1, 2)
	timerMicroSeconds = StopMSTimer(timerRefNum)
	testScores [testNum] = timerMicroSeconds/1E6
	testNum +=1
	DoWindow/T twoPxop_ProjFrames_out_y "Project All Frames Y Mean"
	doupdate;sleep/S 1
	
	// Y Median
	testType [testNum]="Project All Frames dim =Y, mode = Median"
	timerRefNum = StartMSTimer
	ProjectAllFrames (theStack, "root:ProjectFrames_output_y", 1, 1, 3)
	timerMicroSeconds = StopMSTimer(timerRefNum)
	testScores [testNum] = timerMicroSeconds/1E6
	testNum +=1
	DoWindow/T twoPxop_ProjFrames_out_y "Project All Frames Y Median"
	doupdate;sleep/S 1
	
	// Project Y slice
	WAVE ProjectFrames_output = root:ProjectFrames_output_y
	ProjectFrames_output = 0
	testType [testNum]="Project Y Slice x 10"
	testScores [testNum] = 0
	for (iSpec = 0; iSPec < 10; iSpec += 1)
		timerRefNum = StartMSTimer
		ProjectYSlice (theStack, ProjectFrames_output, 100 + iSpec)
		timerMicroSeconds = StopMSTimer(timerRefNum)
		testScores [testNum] += timerMicroSeconds/1E6
		DoWindow/T twoPxop_ProjFrames_out_y "Project Y Slice " + num2str (100 + iSpec)
		doupdate;sleep/S 0.5
	endfor
	testNum +=1
	
	// Project Frames Z
	// Z min
	testType [testNum]="Project All Frames dim =Z, mode = Min"
	timerRefNum = StartMSTimer
	ProjectAllFrames (theStack, "root:ProjectFrames_output_z", 2, 1, 0)
	timerMicroSeconds = StopMSTimer(timerRefNum)
	testScores [testNum] = timerMicroSeconds/1E6
	testNum +=1
	execute "twoPxop_ProjFrames_out_z()"
	doupdate;sleep/S 1
	
	// Z Max
	testType [testNum]="Project All Frames dim =Z, mode = Max"
	timerRefNum = StartMSTimer
	ProjectAllFrames (theStack, "root:ProjectFrames_output_z", 2, 1, 1)
	timerMicroSeconds = StopMSTimer(timerRefNum)
	testScores [testNum] = timerMicroSeconds/1E6
	testNum +=1
	DoWindow/T twoPxop_ProjFrames_out_z "Project All Frames Z Maximum"
	doupdate;sleep/S 1
	
	//Z Average
	testType [testNum]="Project All Frames dim =Z, mode = Mean"
	timerRefNum = StartMSTimer
	ProjectAllFrames (theStack, "root:ProjectFrames_output_z", 2, 1, 2)
	timerMicroSeconds = StopMSTimer(timerRefNum)
	testScores [testNum] = timerMicroSeconds/1E6
	testNum +=1
	DoWindow/T twoPxop_ProjFrames_out_z "Project All Frames Z Mean"
	doupdate;sleep/S 1
	
	// Z Median
	testType [testNum]="Project All Frames dim =Z, mode = Median"
	timerRefNum = StartMSTimer
	ProjectAllFrames (theStack, "root:ProjectFrames_output_z", 2, 1, 3)
	timerMicroSeconds = StopMSTimer(timerRefNum)
	testScores [testNum] = timerMicroSeconds/1E6
	testNum +=1
	DoWindow/T twoPxop_ProjFrames_out_z "Project All Frames Z Median"
	doupdate;sleep/S 1
	
	// Project Z slice
	WAVE ProjectFrames_output = root:ProjectFrames_output_z
	ProjectFrames_output = 0
	testType [testNum]="Project Z Slice x 10"
	testScores [testNum] = 0
	for (iSpec = 0; iSPec < 10; iSpec += 1)
		timerRefNum = StartMSTimer
		ProjectZSlice (theStack, ProjectFrames_output, 100 + iSpec)
		timerMicroSeconds = StopMSTimer(timerRefNum)
		testScores [testNum] += timerMicroSeconds/1E6
		DoWindow/T twoPxop_ProjFrames_out_z "Project Z Slice " + num2str (100 + iSpec)
		doupdate;sleep/S 0.5
	endfor
	testNum +=1
			
	DoAlert /T="Testing twoPhotonXOP" 0, "Project Frames"
	// close all 2p graph windows
	closeTwoPGraphs()
	doUpdate
	
	//Convolve frames
	testType [testNum]="Convolve Frames w=5"
	WAVE gwave = makeKernel(5) // a free wave for kernel
	timerRefNum = StartMSTimer
	ConvolveFrames (theStack, "root:Convolve_Out", 0, gwave, 1)
	timerMicroSeconds = StopMSTimer(timerRefNum)
	testScores [testNum] = timerMicroSeconds/1E6
	testNum += 1
	execute "twoPxop_Convole_Out()"
	doupdate;sleep/S 1
	
	testType [testNum]="Convolve Frames w=11"
	wave gwave = makeKernel(11) // a free wave for kernel, bigger
	timerRefNum = StartMSTimer
	ConvolveFrames (theStack, "root:Convolve_Out", 0, gwave, 1)
	timerMicroSeconds = StopMSTimer(timerRefNum)
	testScores [testNum] = timerMicroSeconds/1E6
	testNum += 1
	DoWindow/T twoPxop_Convole_Out "Gaussian Convolve Width = 11"
	doupdate;sleep/S 1
	
	testType [testNum]="Symetrical Convolve Frames w=11"
	WAVE gwave = makeSymkernel (11)
	timerRefNum = StartMSTimer
	SymConvolveFrames (theStack, "root:Convolve_Out", 0, gwave, 1)
	timerMicroSeconds = StopMSTimer(timerRefNum)
	testScores [testNum] = timerMicroSeconds/1E6
	testNum += 1
	DoWindow/T twoPxop_Convole_Out "Symetrical Gaussian Convolve Width = 11"
	doupdate;sleep/S 1
	
	testType [testNum]="Symetrical Convolve Frames w=23"
	WAVE gwave = makeSymkernel (23)
	timerRefNum = StartMSTimer
	SymConvolveFrames (theStack, "root:Convolve_Out", 0, gwave, 1)
	timerMicroSeconds = StopMSTimer(timerRefNum)
	testScores [testNum] = timerMicroSeconds/1E6
	testNum += 1
	DoWindow/T twoPxop_Convole_Out "Symetrical Gaussian Convolve Width = 23"
	doupdate;sleep/S 1
	
	testType [testNum]="Median Frames w=5"
	timerRefNum = StartMSTimer
	MedianFrames (theStack,  "root:Convolve_Out", 5, 1)
	timerMicroSeconds = StopMSTimer(timerRefNum)
	testScores [testNum] = timerMicroSeconds/1E6
	testNum += 1
	DoWindow/T twoPxop_Convole_Out "Median Frames Width = 5"
	doupdate;sleep/S 1
	
	DoAlert /T="Testing twoPhotonXOP" 0, "Convolve Frames"
	
	// LSM utilities
	WAVE convolveOut = root:Convolve_Out
	testType [testNum]="SwapEven x 2"
	testScores [testNum] =0
	timerRefNum = StartMSTimer
	swapEven (convolveOut)
	timerMicroSeconds = StopMSTimer(timerRefNum)
	testScores [testNum] += timerMicroSeconds/1E6
	DoWindow/T twoPxop_Convole_Out "Swap Even Once"
	doupdate;sleep/S 2
	timerRefNum = StartMSTimer
	swapEven (convolveOut)
	timerMicroSeconds = StopMSTimer(timerRefNum)
	testScores [testNum] += timerMicroSeconds/1E6
	testNum += 1
	DoWindow/T twoPxop_Convole_Out "Swap Even Twice"
	doupdate;sleep/S 2
	
	testType [testNum]="Transpose Frames x 2"
	testScores [testNum] =0
	timerRefNum = StartMSTimer
	TransposeFrames(convolveOut)
	timerMicroSeconds = StopMSTimer(timerRefNum)
	testScores [testNum] += timerMicroSeconds/1E6
	DoWindow/T twoPxop_Convole_Out "Transpose Frames Once"
	doupdate;sleep/S 2
	timerRefNum = StartMSTimer
	TransposeFrames(convolveOut)
	timerMicroSeconds = StopMSTimer(timerRefNum)
	testScores [testNum] += timerMicroSeconds/1E6
	testNum += 1
	DoWindow/T twoPxop_Convole_Out "Transpose Frames Twice"
	doupdate;sleep/S 2
	DoAlert /T="Testing twoPhotonXOP" 0, "LSM utilities"
	
end

	
	
function tse ()
// LSM utilities
	WAVE convolveOut = root:Convolve_Out
	variable timerMicroSeconds
	//testType [testNum]="SwapEven x 2"
	//testScores [testNum] =0
	timerMicroSeconds =0
	variable timerRefNum = StartMSTimer
	TransposeFrames (convolveOut)
	timerMicroSeconds += StopMSTimer(timerRefNum)
	//testScores [testNum] += timerMicroSeconds/1E6
	DoWindow/T twoPxop_Convole_Out "Transpose Once"
	doupdate;sleep/S 2
	timerRefNum = StartMSTimer
	TransposeFrames (convolveOut)
	timerMicroSeconds += StopMSTimer(timerRefNum)
	//testScores [testNum] += timerMicroSeconds/1E6
	//testNum += 1
	DoWindow/T twoPxop_Convole_Out "Transpose Twice"
	doupdate;sleep/S 2
	print timerMicroSeconds/1E6
end
	
	
end




Window twoPxop_theStack() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(0,66,528,394) as "3D Image Stack"
	AppendImage/T theStack
	ModifyImage theStack ctab= {0,4096,Rainbow,1}
	ModifyGraph margin(left)=14,margin(bottom)=14,margin(top)=14,margin(right)=14,gFont="Helvetica"
	ModifyGraph width={Plan,1,top,left}
	ModifyGraph mirror=2
	ModifyGraph nticks(left)=7,nticks(top)=10
	ModifyGraph font="Helvetica"
	ModifyGraph minor=1
	ModifyGraph fSize=9
	ModifyGraph standoff=0
	ModifyGraph tkLblRot(left)=90
	ModifyGraph btLen=3
	ModifyGraph tlOffset=-2
EndMacro


Window twoPxop_KalmanAllFrames_out() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(0,422,528,750) as "Kalman All Frames Output"
	AppendImage/T KalManAllFrames_out
	ModifyImage KalManAllFrames_out ctab= {0,4096,Rainbow,1}
	ModifyGraph margin(left)=14,margin(bottom)=14,margin(top)=14,margin(right)=14,gFont="Helvetica"
	ModifyGraph mirror=2
	ModifyGraph nticks(left)=7,nticks(top)=10
	ModifyGraph minor=1
	ModifyGraph fSize=9
	ModifyGraph standoff=0
	ModifyGraph tkLblRot(left)=90
	ModifyGraph btLen=3
	ModifyGraph tlOffset=-2
EndMacro

Window twoPxop_KalmanSpecFrames_out() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(529,66,1057,394) as "Kalman Spec Frames Output"
	AppendImage/T KalmanSpecFrames_out
	ModifyImage KalmanSpecFrames_out ctab= {0,4096,Rainbow,1}
	ModifyGraph margin(left)=14,margin(bottom)=14,margin(top)=14,margin(right)=14,gFont="Helvetica"
	ModifyGraph mirror=2
	ModifyGraph nticks(left)=7,nticks(top)=10
	ModifyGraph minor=1
	ModifyGraph fSize=9
	ModifyGraph standoff=0
	ModifyGraph tkLblRot(left)=90
	ModifyGraph btLen=3
	ModifyGraph tlOffset=-2
EndMacro

Window twoPxop_KalmanWaveToFrame_out() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(529,423,1057,751) as "Kalman Wave To Frame Output"
	AppendImage/T KalmanWaveToFrame_out
	ModifyImage KalmanWaveToFrame_out ctab= {0,4096,Rainbow,1}
	ModifyGraph margin(left)=14,margin(bottom)=14,margin(top)=14,margin(right)=14,gFont="Helvetica"
	ModifyGraph mirror=2
	ModifyGraph nticks(left)=7,nticks(top)=10
	ModifyGraph minor=1
	ModifyGraph fSize=9
	ModifyGraph standoff=0
	ModifyGraph tkLblRot(left)=90
	ModifyGraph btLen=3
	ModifyGraph tlOffset=-2
EndMacro


Window twoPxop_ProjFrames_out_X() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(4.5,40.25,336,250.25) as "Project All Frames X Minimum"
	AppendImage/G=1/T ProjectFrames_output_x
	ModifyImage ProjectFrames_output_x ctab= {0,4096,Rainbow,1}
	ModifyGraph margin(left)=14,margin(bottom)=14,margin(top)=14,margin(right)=14,gFont="Arial"
	ModifyGraph width={Plan,1,top,left}
	ModifyGraph mirror=2
	ModifyGraph nticks(top)=7
	ModifyGraph minor=1
	ModifyGraph fSize=9
	ModifyGraph standoff=0
	ModifyGraph tkLblRot(left)=90
	ModifyGraph btLen=3
	ModifyGraph tlOffset=-2
EndMacro

Window twoPxop_ProjFrames_out_y() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(346.5,40.25,983.25,250.25) as "Project All Frames Y Minimum"
	AppendImage/G=1/T ProjectFrames_output_y
	ModifyImage ProjectFrames_output_y ctab= {0,4096,Rainbow,1}
	ModifyGraph margin(left)=14,margin(bottom)=14,margin(top)=14,margin(right)=14,gFont="Arial"
	ModifyGraph width={Plan,1,top,left}
	ModifyGraph mirror=2
	ModifyGraph nticks(top)=7
	ModifyGraph minor=1
	ModifyGraph fSize=9
	ModifyGraph standoff=0
	ModifyGraph tkLblRot(left)=90
	ModifyGraph btLen=3
	ModifyGraph tlOffset=-2
EndMacro

Window twoPxop_ProjFrames_out_z() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(346.5,273.5,984,605.75) as "Project All Frames Z Minimum"
	AppendImage/G=1/T ProjectFrames_output_z
	ModifyImage ProjectFrames_output_z ctab= {0,4096,Rainbow,1}
	ModifyGraph margin(left)=14,margin(bottom)=14,margin(top)=14,margin(right)=14,gFont="Arial"
	ModifyGraph width={Plan,1,top,left}
	ModifyGraph mirror=2
	ModifyGraph nticks(top)=7
	ModifyGraph minor=1
	ModifyGraph fSize=9
	ModifyGraph standoff=0
	ModifyGraph tkLblRot(left)=90
	ModifyGraph btLen=3
	ModifyGraph tlOffset=-2
EndMacro


Window twoPxop_ProjFrames_out() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(2,66,841,635) as "Project All Frames X Minimum"
	AppendImage/G=1/T ProjectFrames_output
	ModifyImage ProjectFrames_output ctab= {0,4096,Rainbow,1}
	ModifyGraph margin(left)=14,margin(bottom)=14,margin(top)=14,margin(right)=14,gFont="Helvetica"
	ModifyGraph width={Plan,1,top,left}
	ModifyGraph mirror=2
	ModifyGraph nticks(top)=7
	ModifyGraph minor=1
	ModifyGraph fSize=9
	ModifyGraph standoff=0
	ModifyGraph tkLblRot(left)=90
	ModifyGraph btLen=3
	ModifyGraph tlOffset=-2
EndMacro

Window twoPxop_Convole_Out() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(6,41,643.5,373.25) as "Convolve Gaussian Width = 5"
	AppendImage/T Convolve_Out
	ModifyImage Convolve_Out ctab= {0,4096,Rainbow,1}
	ModifyImage Convolve_Out plane= 30
	ModifyGraph margin(left)=14,margin(bottom)=14,margin(top)=14,margin(right)=14,gFont="Arial"
	ModifyGraph width={Plan,1,top,left}
	ModifyGraph mirror=2
	ModifyGraph nticks(left)=7,nticks(top)=10
	ModifyGraph minor=1
	ModifyGraph fSize=9
	ModifyGraph standoff=0
	ModifyGraph tkLblRot(left)=90
	ModifyGraph btLen=3
	ModifyGraph tlOffset=-2
	SetAxis/A/R left
EndMacro

