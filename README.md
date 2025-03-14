## tfdnoise.py

C++ STFT Median Thresholding Filter, based on the FreeUSP TFDNoise utility. The filtering function itself is located in ```tfd.cpp```. The utility ``tfdnoise`` takes a SEG-Y file, filters all gathers within it, and stores the result in a new SEG-Y file, while also creating a difference SEG-Y file. Parameters for the filtering process are stored in ```tfdparams.py```. The utility ```makepics.py``` generates a comparison plot of these three SEG-Y files for specified trace numbers.

First, download, build and install:

```FFTW3``` for spectrum calcuation: https://github.com/FFTW/fftw3.git

```SegyIO``` SEG-Y reading: https://github.com/equinor/segyio.git

```Indicators``` library headers are already in ```include``` folder

FreeUSP Toolkit https://stuartschmitt.com/FreeUSP/

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
./tfdnoise ../params.tfd
```
