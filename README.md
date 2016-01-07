KappaNEURON
===========

KappaNEURON integrates the [SpatialKappa][SpatialKappa] simulator with
[NEURON][NEURON] to allow rule-based simulations of molecular systems
embedded in neurons.

For example, the package facilitates simulation of dynamic models of
the postsynaptic proteome in the context of the spine head.

![KappNEURON demonstration simulation of postsynaptic proteome in context of spine head - first 6 seconds](doc/figs/neuron_kappa_Very_short_6000.png)
![KappNEURON demonstration simulation of postsynaptic proteome in context of spine head - first 65 seconds](doc/figs/neuron_kappa_Very_short_65000.png)

Further details
---------------

* Sterratt, D. C., Sorokina, O. and Armstrong, J. D. (in
  press). ‘Integration of rule-based models and compartmental models
  of neurons’. Lecture Notes in Bioinformatics 7699.  Presented to the
  Third International Workshop on Hybrid Systems Biology Vienna,
  Austria, July 23-24, 2014 at the International Conference on
  Computer-Aided Verification 2014.  Preprint at <a title="Abstract"
  href="http://arxiv.org/abs/1411.4980">arXiv:1411.4980</a>
* [Presentation given to Hybrid Systems Biology 2014](doc/2014-07-24-rb-compartmental-method.pdf)

Installation & testing
------------------------

* First install NEURON from source according to the commands in
  [doc/INSTALL-neuron][INSTALL-neuron]
* Then follow the instructions in [doc/INSTALL.md][INSTALL]

Authorship & License
--------------------

KappaNEURON is Copyright © 2013-2014 David C. Sterratt
<<david.c.sterratt@ed.ac.uk>> and is released under the
[GPL Version 3](http://www.gnu.org/copyleft/gpl.html).

We thank Anatoly Sorokin for his help with SpatialKappa and Vincent
Danos for helpful discussions.

Acknowledgements
----------------

The research leading to these results has received funding from the
European Union Seventh Framework Programme (FP7/2007-2013) under grant
agreement nos. 241498 (EUROSPIN project), 242167 (SynSys-project) and
604102 (Human Brain Project). 

[SpatialKappa]: https://github.com/davidcsterratt/SpatialKappa "SpatialKappa"

[NEURON]: http://neuron.yale.edu/neuron/ "NEURON"

[INSTALL-neuron]: doc/INSTALL-neuron "NEURON installation instructions"

[INSTALL]: doc/INSTALL.md "KappaNEURON installation instructions"

<!--  LocalWords:  KappaNEURON SpatialKappa KappNEURON Sterratt Danos
 -->
<!--  LocalWords:  Anatoly Sorokin FP EUROSPIN SynSys
 -->
