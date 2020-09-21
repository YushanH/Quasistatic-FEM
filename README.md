# Quasistatic-FEM
This project implements quasistatic FEM solver for hyperelastic deformation modeling.

Three different hyperelastic models are included: Linear elasticity, St-Venant Kirchoff and Corotated elasticity.

Conjugate gradient method for solving linear equations is included.

The results of rest-state grid positions are written to Houdini-compatible obj files.

The example result below is the rest state of a size 1*/1 2D mesh(50*/50 resolution) streched to 1.15 times horizontally, with gravity in effect 
<p align="center">
  <img src="https://github.com/YushanH/Quasistatic-FEM/blob/master/result.png" width="50%" >
</p>
