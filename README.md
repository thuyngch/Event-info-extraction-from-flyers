# Event-info-extraction-from-flyers

Extract key information (i.e., time and location) of an event from a poster/flyer using Digital image processing technique. In this project, we adopt and develop the assignment project 'Event Info Extraction from Mobile Camera Images' from [Standford students](https://web.stanford.edu/class/ee368/Project_Winter_1314/)

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

What things you need to install the software and how to install them

```
Ubuntu 16.04, the OS that we use to run this system and we also recommend you to run the system in this OS.

MATLAB R2016b or higher versions, the main IDE tool to run the algorithm.

Tesseract-OCR, the text-extraction tool used in the system.
```

### Installing

In Ubuntu 16.04, you can install Tesseract by the following command:

```
sudo apt-get install tesseract
```

## Project organization

There are two main folders **Document** and **MATLAB**.

Folder **Document** contains documents, e.g. *Report.pdf* (of Standford students) and *Vietnamese translation.pdf* (our translation). 

Folder **MATLAB** is the whole of system. Within it, *fnc* is old algorithm, *refine* is improved algorithm, *im* is the image database, *result* is the buffer to store files generated during running the system, *gui* contains a GUI of the system, *quality.xlsx* is the metric of two algorithms.

## Deployment

At the beginning, you must change the current folder to *MATLAB* folder and run the command ```startup``` to prepare for needed components for the system.

Subsequently, there are two way to run the system. First, you can use the GUI by typing the command ```gui``` to run the system with the improved algorithm. Second, you may implement the old algorithm by the command ```run_im```, however, keep in mind that the you must change the path of image that you want to run in the file *run_im.m*.

[See demo video](https://youtu.be/TrWnVXtxauQ).

## Team members

* **Do Tieu Thien** - *Leader*, *Pre-processing*
* **Nguyen Chinh Thuy** - *Text Detection*, *Algorithm improvement*
* **Le Van Hoang Phuong** - *Post-processing*, *Tesseract-OCR*
