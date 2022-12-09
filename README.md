# Modular Modeling of Analog Organic Neuromorphic Circuits: Towards Prototyping of Hardware-Level Spiking Neural Networks

![](https://img.shields.io/badge/license-GPL-blue.svg)
![DOI](https://img.shields.io/badge/DOI-10.1109%2FTCSI.2022.3226163-brightgreen.svg)

GitHub Repository Detailing the Software Code for *'Modular Modeling of Analog Organic Neuromorphic Circuits: Towards Prototyping of Hardware-Level Spiking Neural Networks'*, published in IEEE Transactions on Circuits and Systems I: Regular Papers (IEEE TCAS-I).

## Abstract
This work proposes a novel modeling approach for analog organic circuits using very simple to customize circuit topology and parameters of individual p- and n-type organic field effect transistors (OFETs). Aided with the combination of primitive elements (OFETs, capacitors, resistors), the convoluted behavior of analog organic neuromorphic circuits (ONCs) and even other general analog organic circuits, can be predicted. The organic log-domain integrator (oLDI) synaptic circuit, the organic differential-pair integrator (oDPI) synaptic circuit, and the organic Axon-Hillock (oAH) somatic circuit are designed and serve as the modular circuit primitives of more complicated ONCs. We first validate our modeling approach by comparing the simulated oDPI and oAH circuit responses to their experimental measurements. Thereafter, the summation effects of the excitatory and inhibitory oDPI circuits in prototyped ONCs are investigated. We also predict the dynamic power dissipation of modular ONCs and show an average power consumption of 2.1 μJ per spike for the oAH soma at a ~1 Hz spiking frequency. Furthermore, we compare our modeling approach with other two representative organic circuit models and prove that our approach outperforms the other two in terms of accuracy and convergence speed.

## Citation

To cite our work, kindly use the following BibTex entry:

```
@ARTICLE{9976306,
author={Y. Yang and M. J. Mirshojaeian Hosseini and W. Kruger and R. A. Nawrocki},
journal={IEEE Transactions on Circuits and Systems I: Regular Papers},
title={Modular Modeling of Analog Organic Neuromorphic Circuits: Toward Prototyping of Hardware-Level Spiking Neural Networks},
year={2022},
doi={10.1109/TCSI.2022.3226163},
ISSN={1549-8328},
}
```

## Get Started

### Prerequisites
- Matlab R2021a or higher version
- Simulink with Simscape toolbox installed
- Cadence Virtuoso IC618

### Start Simulations

- Fit device parameters of p- and n-OFETs using `./p-OFET/pofet_main.m` and `./n-OFET/nofet_main.m`. Note: Manually set the inflexion point that separates the double-exponential region and the above-threshold region in `./p-OFET/ID_Marinov_m2.m` and `./n-OFET/ID_Marinov_m2.m`.
- Set device parameter values in `./oAH/oAH_main1.m`, `./oAH/oAH_main2.m`, `./oAH/oAH_main3.m` and `./oAH/oAH_main4.m` and start the oAH somatic circuit simulation. Note: Manually set the inflexion point in `./oAH/P_Marinov_OFET.ssc` and `./oAH/N_Marinov_OFET.ssc`. The other device parameters are set in the main scripts.
- Set device parameter values in `./oLDI/oLDI_main_Vdd.m`, and `./oLDI/oLDI_main_Vtau.m` and start the oLDI synaptic circuit simulation.
- Set device parameter values in `./oDPI/oDPI_main_Vw1.m`, `./oDPI/oDPI_main_Vw2.m`, `./oDPI/oDPI_main_Vtau.m`, and `./oDPI/oDPI_main_Vth.m` and start the oDPI synaptic circuit simulation.
- Set device parameter values in `./ONCs/SNN_run_main11.m`, `./ONCs/SNN_run_main22.m` and `./ONCs/SNN_run_main33.m` and start the simulation of the ONCs in three configurations.
- Import the Verilog-A folder as an OFET_benchmark library in Cadence Virtuoso and start the DC simulation of single p-OFET or n-OFET device. Open the schematic of the resistive load inverter circuit in Virtuoso and start RL circuit simulation.





## License
All code is licensed under the GNU General Public License v3.0. Details pertaining to this are available at: https://www.gnu.org/licenses/gpl-3.0.en.html

