# Notes on NEURON rxd interface and KappaNEURON implementation

## Table of Contents

1.  [Overview of rxd module in NEURON 7.4](#org56aaac9)
2.  [Overview of KappaNEURON in NEURON 7.4](#org78be12a)
    1.  [`Kappa` class](#orgee9ea18)
    2.  [`KappaFlux` class](#orgdb7386c)
3.  [`generalisedReaction` methods not overridden by KappaNEURON classes](#orgec03626)
    1.  [`_update_indices()`](#org8a7dc56)
    2.  [`_mult`](#org3ba9ea3)
    3.  [`_setup_membrane_fluxes(node_indices, cur_map)`](#org54cc755)
4.  [`generalisedReaction` methods overridden by `Kappa` and `KappaFlux` classes that return null](#org87c4f1e)
    1.  [`_evaluate(states)`](#orgf2da5f3)
    2.  [`_jacobian_entries()`](#orgde405c0)
5.  [`MultiCompartmentalReaction` variables and methods to be overridden in KappaFlux class](#orge008524)
    1.  [`_get_memb_flux()`](#orga0eb36a)
    2.  [`_memb_scales`](#org32d5cf0)
    3.  [`_do_memb_scales()`](#orga95f3db)
    4.  [`_sources`, `_dests`](#org8a2622f)
    5.  [`_cur_ptrs`](#org0d79a76)
    6.  [`_cur_mapped`](#org8ae7352)
6.  [`rxd.py` functions](#orga7c5631)
    1.  [`_rxd_reaction(states)`](#orga3faa09)
7.  [Variables](#org4f2495a)
8.  [NEURON integration procedure](#orgb9504f2)



<a id="org56aaac9"></a>

## Overview of rxd module in NEURON 7.4

The `rxd` module sets up biochemical reactions using a
`generalizedReaction` abstract class from which classes such as
`Reaction`, `Rate` and `MultiCompartmentReaction` are derived. These
classes associate `Species` objects with `Region` objects and set up
the equations for Species concentrations that are passed to NEURON's
solver. 

Communication with the NEURON solver is achieved via the
`neuron.nonvint_block_supervisor` module, which is used, in
`rxd._do_nbs_register()`, to register a vector callback functions
defined in the `rxd` module.

The `generalizedReaction._get_memb_flux()` method is called
from `rxd._currents()`, which passes information about induced
currents back to NEURON's solver.


<a id="org78be12a"></a>

## Overview of KappaNEURON in NEURON 7.4

Essentially, KappaNEURON sets up Kappa simulations and maps the
variables in these simulations to variables that exist in NEURON. At
the start of a deterministic timestep, transmembrane fluxes that
originate in NEURON, e.g. calcium transmembrane currents, are passed
to the the appropriate variable in the appropriate Kappa
simuation. The Kappa simulation is then advanced in time, and the
change in variables is sent back to NEURON at the end of the time
step, and mapped back to the relevant variable in NEURON. To achieve
all this, there is quite a lot of fiddly mapping of data structures
from Kappa to NEURON and vice-versa, and various other NEURON
functions have to be overridden.

The KappaNEURON module overwrites the fourth callback function
(`rxd._fixed_step_solve()`) defined by `rxd._callbacks` by replacing it
with a funtion `_kn_fixed_step_solve()`. At each
timestep of the deterministic solver `_kn_fixed_step_solve()` is called,
and carries out the following:

1. Use `states = rxd._node_get_states()` to get a reference to the state variables
   (concentrations and voltages) in the NEURON solver. The length of `states` is
   the product of the number of segments and the number of state
   variables in the model.
2. Use `b = rxd._rxd_reaction(states)` to get fluxes across the membrane
3. Call `states = _run_kappa_continuous(states, b, dt)`. This function:
   1. Converts the units of `b` into absolute  transition rates, which 
      are passed to the creation rules of each Kappa simulation (stored in
      a `Kappa` object)
   2. Passes the voltage to Kappa
   3. Runs Kappa for one timestep `dt`
   4. To each `KappaFlux` appends the differences between the flux specified
      by `b` and the actual average flux generated over the time step by the
      random Kappa process. The transmembrane fluxes induced by any Kappa rules that
      control ionic flow through the membrane are stored in the `_memb_flux`
      variable of `KappaFlux` objects which are members of the associated
      `Kappa` object.
   5. Updates the continuous variables using `rxd._reaction_matrix_solve()`
   6. Ensures that concentration states are the same as the simulated ones.
4. Call `rxd._section1d_transfer_to_legacy()` to ensure states in rxd module are
   returned to NEURON legacy solver.
   
To achieve the model specification, and mapping of variables the
module defines two classes `Kappa`, derived from `GeneralizedReaction`
and `KappaFlux`, derived from `MultiCompartmentalReaction`.
<a id="orgee9ea18"></a>

### `Kappa` class

This class reads a Kappa file, and starts the Kappa solver (currently
SpatialKappa). It also creates `KappaFlux` objects for each species
that crosses the membrane. Typically this will be calcium. 

When a `Kappa` object is instantiated, it is registered with
`neuron.rxd.rxd._register_reaction()`, which means that various of its
member functions are called at each timestep. As these functions are
used to set up the derivatives and Jacobian in deterministic
equations, they will mostly return zeros.


<a id="orgdb7386c"></a>

### `KappaFlux` class

The `KappaFlux` object references the `Kappa` object, and implements
`_get_memb_flux()` to communicate fluxes at each NEURON time
step. This returns `_memb_flux` mentioned above.


<a id="orgec03626"></a>

## `generalisedReaction` methods not overridden by KappaNEURON classes


<a id="org8a7dc56"></a>

#### `_update_indices()`

-   Sets up:
    -   **`self._indices_dict`:** Mapping from species onto indices in the
        state vector.
    -   **`self._mult`        :** Multiplier used in `rxd._rxd_reaction()`
        to convert output of `self._evaluate()` into units of mM/ms
        when `self._trans_membrane` is true.


<a id="org3ba9ea3"></a>

#### `_mult`

Set in `generalizedReaction._update_indices()`, its value depends on the value of
  `self._trans_membrane` and `self._scale_by_area`:

-   **if** `self._trans_membrane` is `True` and `self.scale_by_area` is `True`
    
    `-areas / volumes / molecules_per_mM_um3` for `source_indices`
    
    `areas / volumes / molecules_per_mM_um3` for `dest_indices`
    
    Units are mol/mM/µm<sup>3</sup>, which are dimensionless

-   **if** `self._trans_membrane` is `True` and `self.scale_by_area` is
    `False`
    
    `-1 / volumes / molecules_per_mM_um3` for `source_indices`
    
    `1 / volumes / molecules_per_mM_um3` for `dest_indices`

-   **if** `self._trans_membrane` is `False`
    
    `-1` for `source_indices`
    
    `1` for `dest_indices`


<a id="org54cc755"></a>

#### `_setup_membrane_fluxes(node_indices, cur_map)`

-   Set up `node_indices`, indices mapping source species
    `self._sources` and destination sepcies `self._dests` to segment
    indices


<a id="org87c4f1e"></a>

## `generalisedReaction` methods overridden by `Kappa` and `KappaFlux` classes that return null


<a id="orgf2da5f3"></a>

#### `_evaluate(states)`

-   Evaluates states to give rate of change of states

-   **Returns**
    -   **`self._indices`    :** Indicies of state variables in the reaction
        in the global state vector
    -   **`self._mult`       :** Set in `_update_indices()`
    -   **`self._rate(*args)`:** Rate of change of states in mM/ms (by
        default) or molecules um<sup>-2</sup> ms<sup>-1</sup> (in
        `multiCompartmentReaction` class, where `self._trans_membrane`
        is `True`)
-   **KappaNEURON implementation**: As algorithm does not utilise the rate
    of change of states, returns three empty lists.


<a id="orgde405c0"></a>

#### `_jacobian_entries()`

-   **Returns**
    -   **`_jac_rows`:** Indicies of rows of entries
    -   **`_jac_cols`:** Indicies of columns of entries
    -   **`data`     :** Values of Jacobian entries in ms<sup>\(-1\)</sup>
-   **KappaNEURON implementation**: As algorithm does not add to Jacobian
    entries, returns three empty lists.


<a id="orge008524"></a>

## `MultiCompartmentalReaction` variables and methods to be overridden in KappaFlux class

-   This needs KappaNEURON to provide methods for any code in loops in
    which `_all_reactions` is iterated over.


<a id="orga0eb36a"></a>

#### `_get_memb_flux()`

-   Gets flux across the membrane due to univalent ion in mA/cm<sup>2</sup>
-   **Returns**: `self._memb_scales*rates` where `rates` comes from
    `_evaluate()` and is in molecules/µm<sup>2</sup>/ms and `self._memb_scales`
    gives the charge per molecule.
-   Thus the units returned by `get_memb_flux()` are 
    
    10<sup>-14</sup> C molecules<sup>-1</sup> molecules um<sup>-2</sup> ms<sup>-1</sup>
    
    = 10<sup>-14</sup> C um<sup>-2</sup> 10<sup>3</sup> s<sup>-1</sup>
    
    = 10<sup>-11</sup> A um<sup>-2</sup>
    
    = 10<sup>-11</sup> A 10<sup>12</sup> 10<sup>-4</sup> cm<sup>-2</sup>
    
    = 10<sup>-3</sup> A cm<sup>-2</sup>
    
    = mA cm<sup>-2</sup>
-   **KappaNEURON implementation:** Picks up flux contained in member
    variable `_memb_flux` which is set in
    `KappaNEURON._kn_fixed_step_solve()` and which results from the net
    change in ions during the preceding time step.


<a id="org32d5cf0"></a>

#### `_memb_scales`

-   Charge per molecule in units of 10<sup>-14</sup> C/molecule, scaled
    according to areas of membranes.
-   Set in `KappaNEURON._do_memb_scales()`
-   This is  `-area_ratios * FARADAY / (10000 * molecules_per_mM_um3)` where
    `area_ratios` is normally 1.
-   Thus scaling factor is `FARADAY / (10000 * molecules_per_mM_um3)`,
     which has units 
    10<sup>4</sup> C mol<sup>-1</sup>/(10<sup>3</sup> molecules mol<sup>-1</sup> dm<sup>3</sup> um<sup>-3</sup>)
    
    = 10<sup>4</sup> C /(10<sup>3</sup> molecules 10<sup>15</sup>) 
    
    = 10<sup>-14</sup> C/molecules


<a id="orga95f3db"></a>

#### `_do_memb_scales()`

-   Set up `self._memb_scales` and sets a current map `rxd._cur_map`
    used by `rxd._update_node_data()` in which is it is passed to
    `Species._setup_currents()`.

-   **KappaNEURON implementation:** Small changes from
    `multicompartmentReaction._do_memb_scales()`: 
    -   Ignore the '`o`' variables in the `cur_map`, otherwise there is a
        crash
    -   Use `numpy.concatenate()` instead of `itertools.chain` command,
        which doesn't seem to work


<a id="org8a2622f"></a>

#### `_sources`, `_dests`

-   List of `weakref` to `SpeciesOnRegion` used in
    `_setup_membrane_fluxes()` and `_kn_fixed_step_solve()`
-   **KappaNEURON implementation:** Set up in `self.__init__()`.


<a id="org0d79a76"></a>

#### `_cur_ptrs`

Pointers to species currents (`ica`, `ik` etc) in `nrn.Segment`.
Set up in `_do_memb_scales()`


<a id="org8ae7352"></a>

#### `_cur_mapped`

Mapping from species conc (`cai`, `cao`) and segment to index in
  `states()`.  Set up in `_do_memb_scales()`


<a id="orga7c5631"></a>

## `rxd.py` functions


<a id="orga3faa09"></a>

#### `_rxd_reaction(states)`

-   Return reaction rate in mM/ms
-   **Returns**
    -   **`self._mult*rate`:** as returned by `_evaluate()`.
-   **Units** If `self._trans_membrane` is `True` (as in
    `multiCompartmentReaction`) then the units of `rate` are molecules
    um<sup>-2</sup> ms<sup>-1</sup>. If `self._scale_by_area` is `True` then the units
    of `self._mult` are um<sup>2</sup>/(um<sup>3</sup> molecules mM<sup>-1</sup> um<sup>-3</sup>). Thus the
    units returned by `_rxd_reaction()` are mM ms<sup>-1</sup>.


<a id="org4f2495a"></a>

## Variables

All in `neuron.rxd` namespace

-   **`rxd._rxd_induced_currents` :** Transmembrane currents induced by
    reactions in `_current()` callback
-   **`rxd._cur_map`                   :** Map from species conc (`cai`,
    `cao`) and `nrn.Segment` to index of `rxd._curr_ptrs`
-   **`rxd._curr_ptrs`                 :** Pointer to species currents
    (`ica`, `ik` etc) in `nrn.Segment`
-   **`rxd._curr_ptr_vector`           :** Pointer from state indicies to
    species currents (`ica`, `ik` etc) in `nrn.Segment`
-   **`rxd._curr_scales` :** Set up in `section1d._setup_currents()`
    called from `species._setup_currents()`.  Converts from current
    density to change in
    concentration. `sign*surface_area[self.indices]*10000./`  
    `(self.species.charge*FARADAY*volumes[self.indices])`
-   **`rxd._curr_indices` :** Mapping from `_curr_ptr_storage` and `_rxd_induced_currents`
    onto state vector b in `_rxd_reaction()`


<a id="orgb9504f2"></a>

## NEURON integration procedure

    // Simplified from nrnoc/fadvance.c:
    void* nrn_fixed_step_thread(NrnThread* nth) {
      double wt;
      deliver_net_events(nth);
      wt = nrnmpi_wtime();
      nrn_random_play(nth);
      nth->_t += .5 * nth->_dt;
      fixed_play_continuous(nth);
      /* Calls nrn_nonvint_block_current() and nrn_nonvint_block_conductance()*/
      setup_tree_matrix(nth); 
      nrn_solve(nth);  /* Solve voltage */
      second_order_cur(nth); n
      update(nth);
      /* Updates t by 0.5dt and calls nrn_nonvint_block_fixed_step_solve*/
      nrn_fixed_step_lastpart(nth); 
      return (void*)0;
    }
    
    
    /* Simplified from nrnoc/treeset.c: */
    /* for the fixed step method */
    void* setup_tree_matrix(NrnThread* _nt){
      nrn_rhs(_nt);
      nrn_lhs(_nt);
      nrn_nonvint_block_current(_nt->end, _nt->_actual_rhs, _nt->id);
      nrn_nonvint_block_conductance(_nt->end, _nt->_actual_d, _nt->id);
      return (void*)0;
    }
    
    /* Simplified from nrnoc/fadvance.c: */
    void* nrn_fixed_step_lastpart(NrnThread* nth) {
      CTBEGIN
      nth->_t += .5 * nth->_dt;
      fixed_play_continuous(nth);
      nrn_extra_scatter_gather(0, nth->id);
      nonvint(nth); 	/* Calls nrn_nonvint_block_fixed_step_solve(_nt->id);*/
      nrn_ba(nth, AFTER_SOLVE);
      fixed_record_continuous(nth);
      CTADD
      nrn_deliver_events(nth) ; /* up to but not past texit */
      return (void*)0;
    }
    
    /* Simplified from nrnoc/fadvance.c: */
    void nonvint(NrnThread* _nt)
    {
      int i;
      double w;
      int measure = 0;
      NrnThreadMembList* tml;
      if (_nt->id == 0 && nrn_mech_wtime_) { measure = 1; }
      errno = 0;
      for (tml = _nt->tml; tml; tml = tml->next) if (memb_func[tml->index].state) {
        Pvmi s = memb_func[tml->index].state;
        if (measure) { w = nrnmpi_wtime(); }
        (*s)(_nt, tml->ml, tml->index);
        if (measure) { nrn_mech_wtime_[tml->index] += nrnmpi_wtime() - w; }
        if (errno) {
          if (nrn_errno_check(i)) {
    hoc_warning("errno set during calculation of states", (char*)0);
          }
        }
      }
      long_difus_solve(0, _nt); /* if any longitudinal diffusion */
      nrn_nonvint_block_fixed_step_solve(_nt->id);
    }

