# CoherentDiffractiveImaging
This repository holds the MATLAB code for a project in Coherent Diffractive Imaging. CDI is a lenless approach to imaging nanoscale objects (much smaller than what a conventional microscope is capable of). This technique scatters coherent beams on the target object. The diffraction pattern is then collected by a detector. The beauty of this process is that the diffraction pattern collected is a Fourier Transform of the original image. The issue with this is that the fourier transform of an image is completely unrecognizable. This process spawns the well documented Phase Problem, which is that when the Fourier Transform is taken the magnitude of the image remains but ALL phase information is lost. This makes the process of reconstructing the original image much more difficult although possible through iterative algorithms. The code in this respository utilizes the **Hybrid Input Output** technique to reconstruct the phase information and in turn reconstruct the original image. 
The First file is reconstructing an artificial image with no issues. 
The Second file is reconstrcuting the same artificial image but with missing information at the center of the diffraction pattern. (Using a beamstop to block the center of a CDI image is common in real CDI experiments because the coherent beam is intense an can damage the detectors)
The third file is reconstructing a real image of a yeast cell. This diffraction pattern also has missing information at the center as well as blurry spots and other issues that come up in CDI experiments. 
