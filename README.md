# image-Source-Reverb

This Matlab code creates an impulse response based on image source modelling within a simulated space and between simulated listener and source positions.  An audio file can be loaded into the code and convolved with the impulse response in order to simulate the audio as the source in the model.

## Getting Started

Open and run the code in Matlab, a GUI will appear and audio can be loaded as the dry signal.  Room and model parameters can be changed in the GUI for different effects.

### Prerequisites

Make sure all files in this repository are located in the chosen Matlab directory.  Any future audio files must be kept in this same directory as well.

## Running the tests

The impulse response is automatically calculated but use "Calculate IR" just in case this does not occur.  You can play and save the IR of the model simulation using the buttons below.  Play wet audio will output the audio convolved with the IR at the modelled positions and save wet audio exports an audio file with the reverb effect active.
  
## Versioning

Version 1.0

### Current Concerns/Issues

* GUI version can be tightened to provide easier use.

## Authors

* **Dan Roth** - *Initial work* 

## Acknowledgments

* Thanks to Dr. Adam Hill at University of Derby for his help and support with this project and his code for calculating wall hits and getFile.m.

