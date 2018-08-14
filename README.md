# HiLo Structured Illumination Microscopy

This is a very basic MATLAB implementation of the HiLo algorithm for structured illumination. 

<img src="https://github.com/raacampbell/HiLo_SIM/blob/master/example_images/before-after-zeiss.jpg" />


## Usage
```
H=hilospeckle(ImageUniform, ImageSpeckle); %Run and show result image

clf, imagesc(H.HiloFinal), axis equal off %re-plot result image

%Change a setting and re-run
H.targetThickness=4; 
H.runHiLo
H.showResults

%Swap the speckle image with the unform image and re-run to compare with above
H.speckleIm = H.uniformIm;
H.runHiLo
figure
H.showResults
```


###
References
* [Optically sectioned in vivo imaging with speckle illumination HiLo microscopy
Daryl Lim, Tim N. Ford, Kengyeh K. Chu, Jerome Mertz, Journal of Biomedical Optics
2012](https://www.spiedigitallibrary.org/journals/Journal-of-Biomedical-Optics/volume-16/issue-1/016014/Optically-sectioned-in-vivo-imaging-with-speckle-illumination-HiLo-microscopy/10.1117/1.3528656.full)
* [ImageJ plugin](http://biomicroscopy.bu.edu/resources/4)
