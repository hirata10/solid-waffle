# Flat simulator for validating solid-waffle

This folder contains the scripts and functions that can be used to generate a simulation with certain known detector effects.  Everything is very much under development, and some features like nonlinear interpixel capacitance is not simulated with the code as it stands.

To run the generate a simulated flat field, you run the command
```python simulate_flat.py <config_file>```
where ```ex_sim_config``` can be used as an example, although it is currently not fully updated with the latest features.
