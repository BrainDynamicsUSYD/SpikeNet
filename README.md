# SpikeNet
SpikeNet is a software that has three standard alone components.
1. User interface for configuring spiking neuronal networks
2. A c++ simulator 
3. User interface for parsing and post-analyzing the simulation results.
The design of SpikeNet provides the following four main features.

* **Configurability** SpikeNet supports any user-defined structure of synaptic connectivity topologies, coupling strengths and conduction delays. It can be easily extended by developers to support any variations of integrate-and-fire neuron and synapse models.

* **Performance**  Simulation of spiking neuronal network quickly becomes computationally intensive if the number of neurons in the network exceeds a few thousand. To achieve superior performance, various measures have been taken at both algorithmic and implementation level. 

* **User-friendly interface** In SpikeNet, although c++ is used for heavy-duty computation, its user-interface is written in  high-level programming language (Matlab) for user-friendliness and fast prototyping. This means SpikeNet does not require non-developer users to be familiar with c++.

* **Scalability** The design of the SpikeNet c++ simulator readily supports parallel computing using Message Passing Interface (MPI). Additionally, the HDF5-based I/O file format provides big data handling capability. Also Portable Batch System (PBS) scripts for array jobs are provided if the you have access to a cluster.


## Getting Started

### Prerequisites
* Autoconf, a tandard tool on OSX and linux distributions
* A c++ compiler that supports c++11 standard (GCC 4.2.1 or later; Intel C++ 12.0 or later)
* HDF5 c/c++ API (open source)
* Matlab (2013a or later) is optional but highly recommended
* Portable Batch System (PBS) is optional but highly recommended

#### FAQ
Q: What if I am using Windows?
A: Sorry, you are on your own. Well, you can always run a Linux virtual machine on Windows.

Q: What if I do not have Matlab or simply hate it?
A: You can either request I/O interface in Python from us or contribute to the project by translating the existing Matlab I/O interface into Python or other langangues.

### Installing
Following are the steps to set up the SpikeNet c++ simulator.
1. Ask for read permission from one of the contributors with admin rights.
\item Make a new directory: \mylstinline{mkdir tmp; cd tmp}
\item Clone SpikeNet: \mylstinline{git clone git@github.com:BrainDynamicsUSYD/SpikeNet.git}
\newline (Method 2: \mylstinline{git clone https://github.com/BrainDynamicsUSYD/SpikeNet})

\item Go into the directory: \mylstinline{cd SpikeNet}
\item Make a copy of the makefile: \mylstinline{cp other_scripts/make*bak makefile}
\item Build the C\texttt{++} simulator:  \mylstinline{make; make clean}
\end{itemize}

Now you should see the ``simulator'' in the current directory, with which you can run simulations by creating input files according to \nameref{sec:IO protocols}.
However, it will be much easier to work with the Matlab user interface.
If you do not have access to Matlab, try to contact the contributors to request interfaces with other high-level programming languages (Python for example).
Following are the steps to use the Matlab user interface.
\begin{itemize}
\item Make a new directory for storing data: \mylstinline{cd ..; mkdir tmp_data; ls}
\item Start Matlab: \mylstinline{matlab -nodisplay}
\item Set up the environment for Matlab: \mylstinline{cd SpikeNet; addpath(genpath(cd));}
\item Generate the example input files: \mylstinline{cd ../tmp_data; main_demo;}
\item Quit Matlab: \mylstinline{quit}
\item Run the simulator with the input files: \mylstinline{cd tmp_data; ../simulator *ygin;}
\item Start Matlab: \mylstinline{cd ..; matlab -nodisplay}
\item Set up the environment for Matlab: \mylstinline{cd SpikeNet; addpath(genpath(cd));}
\item Parse the output files, run some basic post-processing and visualization: 
\newline \mylstinline{cd ../tmp_data; PostProcessYG()}
\item Load the simulation result: \mylstinline{d = dir(`*RYG.mat'); R = load(d(1).name)} (You may need to correct the single quotes in Matlab if you are directly copying the code from here.)
\end{itemize}

For HDF5, please make relevant changes in the makefile to specify the paths where you have installed your HDF5 API. Each of the matlab interfaces has a corresponding HDF5 version. For a complete list of both types of matlab interfaces,  
\begin{itemize}
\item \mylstinline{ls SpikeNet/matlab_interfaces/write2ygin}
\item \mylstinline{ls SpikeNet/matlab_interfaces/write2HDF5}
\end{itemize}


For those who have access to a high-performance computing cluster with PBS, SpikeNet also provides bash script that fully automates the above Matlab $\rightarrow$  C\texttt{++}  $\rightarrow$ Matlab workflow for PBS job array submission. 
The script all\_in\_one.sh has the following features:
\begin{itemize}
\item It automatically detects which stage each array job (with a unique 4-digit integer array ID) has reached: pre-processing done, simulation done or post-simulation data parsing done. The script will start each array job from the last unfinished stage instead of the first stage. This feature comes in handy when hundreds of array jobs end prematurely at different stages, say, due to the HPC being shut down unexpectedly, in which case simply a re-submission of the script will clean up the mess.
\item It passes the array ID as an input argument to the matlab pre-processing script.
\item It automatically saves a copy of the pre-processing Matlab script to the data directory when starting the array job with ID 0001.
\end{itemize}

Following are the steps to use the PBS script.
\begin{itemize}
\item Make sure you have set up your PBS environment correctly (e.g., modelue load HDF5-1.10.0). This can be done by editing shell config file.
\item Go to the tmp directory
\item Make a copy of the script: \mylstinline{cp SpikeNet/other_scripts/all*bak all_in_one.sh}
\item Change it to executable: \mylstinline{chmod +x all_in_one.sh}
\item Edit the following variables in the bash script accordingly: 
\begin{itemize}
\item {\footnotesize MATLAB\_SOURCE\_PATH\_2=`your\_path'}
\item {\footnotesize MATLAB\_PRE\_PROCESS\_FUNC=`your\_functions'}
\end{itemize}
\item Make a directory for PBS output: \mylstinline{mkdir PBSout}
\item Submit the job: \mylstinline{qsub -t 1-X -q queue_name all_in_one.sh}
\end{itemize}


For MPI jobs with SpikeNet, please contact Yifan Gu for more technical details.


```
Give the example
```

And repeat

```
until finished
```

End with an example of getting some data out of the system or using it for a little demo


## Authors

* **Yifan Gu** - *Initial work* - [yigu8115](https://github.com/yigu8115)
* **James A Henderson** - *HDF5-based I/O and learning schemes* - [JamesAHenderson](https://github.com/JamesAHenderson)

See also the list of [contributors](https://github.com/BrainDynamicsUSYD/SpikeNet/graphs/contributors) who participated in this project.

## License

This project is licensed under the Apache License 2.0 - see the [LICENSE.md](LICENSE.md) file for details.


