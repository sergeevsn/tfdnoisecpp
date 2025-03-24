## TFDNOISE.CPP

C++ STFT Median Thresholding Filter, based on the FreeUSP TFDNoise utility. The filtering function itself is located in ```tfd.cpp```. The utility ``tfdnoise`` takes a SEG-Y file, filters all gathers within it, and stores the result in a new SEG-Y file, while also creating a difference SEG-Y file. Parameters for the filtering process are stored in ```tfdparams.py```. The utility ```makepics.py``` generates a comparison plot of these three SEG-Y files for specified trace numbers.

First, download, build and install:

```FFTW3``` for spectrum calcuation: https://github.com/FFTW/fftw3.git

```SegyIO``` SEG-Y reading: https://github.com/equinor/segyio.git

```Indicators``` library headers are already in ```include``` folder

FreeUSP Toolkit https://stuartschmitt.com/FreeUSP/

### Parameter file explanation
```n_fft=64```  - number of samples in STFT window

```trace_aperture=25```  - number of adjacent traces from which the median amplitude is taken for replacing those larger than threshold

```input_file=../test_data/01_test_shots.sgy```  - input SEG-Y file

```output_file=../test_data/02_test_shots_denoised.sgy```   - output SEG-Y file
```threshold_multipliers=0 - 1000, 0.5 - 1, 1.5 - 1```   - you can specify time values versus threshold multipliers, for example if you want to leave the top of gather untouched

```gather_byte=9```  - starting byte of gather identifier in SEG-Y trace header, 9 is usually for Field Number, 21 for CDP and so on

### Install
```bash
git clone https://github.com/sergeevsn/tfdnoisepy.git
cd tfdnoisepy
mkdir build && cd build
cmake ..
make
```

### Test run
```bash
mpirun -np 4 ./tfdnoise ../params.tfd
```
