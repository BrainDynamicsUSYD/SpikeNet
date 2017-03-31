# SpikeNet
SpikeNet is a software for simulating and analyzing spiking neuronal networks, of which the design provides the following four main features.

* **Configurability** SpikeNet supports any user-defined structure of synaptic connectivity topologies, coupling strengths and conduction delays. It can be easily extended by developers to support any variations of integrate-and-fire neuron and synapse models.

* **Performance**  Simulation of spiking neuronal network quickly becomes computationally intensive if the number of neurons $N$ in the network exceeds a few thousand. To achieve superior performance, various measures have been taken at both algorithmic and implementation level. Notably, at the algorithmic level, kinetic models are exclusively chosen for their algorithmic efficiency. Furthermore, for such kinetic models, computational cost of the simulation can be dramatically reduced  with the commonly used mathematical abstraction that synapses can be simplified into a few groups, within which their dynamics (or more specifically, the time constants) are identical. At the implementation level, c++, a programming language renowned for its high-performance computing, is used.

* **User-friendly interface** In SpikeNet, although C++ is used for heavy-duty computation, its user-interface is written in a high-level programming language (Matlab) for user-friendliness and fast prototyping. 

* **Scalability** The design of the SpikeNet C++ simulator readily supports parallel computing using Message Passing Interface (MPI). Additionally, the HDF5-based I/O file format provides big data handling capability. Also Portable Batch System (PBS) scripts for array jobs are provided if the you have access to a cluster.


## Getting Started

### Prerequisites
* Autoconf (a tandard tool on OSX and linux distributions)
* A c++ compiler that supports c++11 standard (GCC 4.2.1 or later; Intel C++ 12.0 or later).
* HDF5 c/c++ API (open source)
* Matlab (2013a or later) is optional but highly recommended.
* Portable Batch System (PBS) is optional but highly recommended.

#### FAQ
Q: What if I am using Windows?
A: Sorry you are on your own. 

Q: What if I do not have or hate Matlab?
A: You can either request I/O interface in Python from us or contribute to the project by translating the existing Matlab I/O interface into Python or other langangues.


### Installing

A step by step series of examples that tell you have to get a development env running

Say what the step will be

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


