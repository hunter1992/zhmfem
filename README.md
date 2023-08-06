# zhmfem
A simple FEM (finite element method) crate in Rust

## Introduction
***zhmfem*** is in the early stages of programming. The main function of the existing code is to construct linear algebraic equations of finite element model and solve them. 

The front end meshing component has not yet been developed (future plans are to use existing meshing tools such as gmesh to provide grid data); The backend component to display the results of the calculation is also not part of the development plan, only the interface for generating VTK data to display the results of the calculation with ParaView is planned.
