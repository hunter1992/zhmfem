<!-- https://github.com/hunter1992/zhmfem -->
<a id="readme-top"></a>

# zhmfem

A finite element calculation command line software based on Rust

<!-- PROJECT SHIELDS -->
<!--
*** I'm using markdown "reference style" links for readability.
*** Reference links are enclosed in brackets [ ] instead of parentheses ( ).
*** See the bottom of this document for the declaration of the reference variables
*** for contributors-url, forks-url, etc. This is an optional, concise syntax you may use.
*** https://www.markdownguide.org/basic-syntax/#reference-style-links
-->
[![Contributors][contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![Issues][issues-shield]][issues-url]
[![Unlicense License][license-shield]][license-url]
[![Author][linkedin-shield]][linkedin-url]



<!-- PROJECT LOGO -->
<br />
<div align="center">
  <a href="https://github.com/hunter1992/zhmfem">
    <img src="imgs/zfem.png" alt="Logo" width="1000" height="391">

  </a>

<h3 align="center">ZFEM</h3>

  <p align="center">
    A finite element calculation command line software based on Rust
    <br />
    <a href="https://github.com/hunter1992/zhmfem"><strong>Explore the docs »</strong></a>
    <br />
    <br />
    <a href="https://github.com/hunter1992/zhmfem?tab=readme-ov-file#example">View Demo</a>
    ·
    <a href="https://github.com/hunter1992/zhmfem/issues/new?labels=bug&template=bug-report---.md">Report Bug</a>
    ·
    <a href="https://github.com/hunter1992/zhmfem/issues/new?labels=enhancement&template=feature-request---.md">Request Feature</a>
  </p>
</div>



<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
      <ul>
        <li><a href="#built-with">Built With</a></li>
      </ul>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation-and-use">Installation and use</a></li>
      </ul>
    </li>
    <li><a href="#tutorial">Tutorial</a></li>
    <li><a href="#roadmap">Roadmap</a></li>
    <li><a href="#contributing">Contributing</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgments">Acknowledgments</a></li>
  </ol>
</details>



<!-- ABOUT THE PROJECT -->
## About The Project

![FEM](imgs/FEM.jpg "FEM application")

After several years of research on the theoretical analysis and application of finite element methods, ZHM initiated this open source project out of interest upon his PhD graduation.

The main goals of this open source project are: 
* to efficiently implement classical finite element algorithms using the Rust language, and to be able to perform engineering calculations with the required efficiency and precision; 
* Explore the application of new finite element algorithms (such as virtual element method, peridynamics method, etc.) in practical problems.

The answers to some basic questions about this project are as follows:
* Why Rust? -- Performance, Reliability and Productivity. [(About Rust Language)](https://www.rust-lang.org/)
* pre/post processing module? -- ZHMFEM is in the early stages of core function development and there are currently no plans for pre/post processing modules. The nodes, cells, loads and boundary conditions required for pre-processing can be set by code; Post-processing is implemented in ParaView after the calculation result is output to vtk file.


ZHMFEM is currently in the early stages of development, and the authors look forward to receiving suggestions from various user groups, and will accept good suggestions after full consideration and discussion.

<p align="right">
(<a href="#readme-top">back to top</a>)</p>



### Built With

* [![Rust][Rustc]][Rust-url]
<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- GETTING STARTED -->
## Getting Started

### Prerequisites

#### 1. Operating system
The author developed the original version of ZHMFEM using the Manjaro Linux system. Currently in the early stages of core functionality development, the compiled version of ZHMFEM is __not__ compatible with MacOS and Windows (you can fix this by compiling the source code yourself on your system). The operating system information used by the author is as follows:
```
 ██████████████████  ████████     zhm@zhm
 ██████████████████  ████████     OS: Manjaro 24.2.1 Yonada
 ██████████████████  ████████     Kernel: x86_64 Linux 6.1.119-1-MANJARO
 ██████████████████  ████████     Uptime: 2h 26m
 ████████            ████████     Packages: 1600
 ████████  ████████  ████████     Shell: zsh 5.9
 ████████  ████████  ████████     Resolution: 1920x1080
 ████████  ████████  ████████     DE: KDE
 ████████  ████████  ████████     WM: KWin
 ████████  ████████  ████████     GTK Theme: Breath [GTK2/3]
 ████████  ████████  ████████     Icon Theme: breeze
 ████████  ████████  ████████     Disk: 199G / 326G (65%)
 ████████  ████████  ████████     CPU: Intel Core i5-8265U @ 8x 3.9GHz [35.0°C]
 ████████  ████████  ████████     GPU: NVIDIA GeForce MX250
                                  RAM: 4701MiB / 7699MiB
```

#### 2. Rust version

ZHMFEM was developed using the neightly version of the Rust language, and the latest version of Rust that currently makes ZHMFEM run is:
```
rustc 1.85.0-nightly (7c002ff9a 2024-12-25)
binary: rustc
commit-date: 2024-12-25
host: x86_64-unknown-linux-gnu
release: 1.85.0-nightly
LLVM version: 19.1.6
```

### Installation and use

1. Install the Rust language environment on your operating system. [install Rust now](https://www.rust-lang.org/tools/install)
2. Clone or download the ZHMFEM repository to local
   ```
   git clone https://github.com/hunter1992/zhmfem.git
   ```
3. Check out the FEM calculation examples located in the 'zhmfem/examples/' path. 
   Run an example with the following command:
   ```
   cargo run --examples <example-name>
   ```
   Or you can compile and run the examples more quickly using the following command:
   ```
   cargo run -j 8 --release --example <example-name>
   ```

<p align="right">(<a href="#readme-top">back to top</a>)</p>


<!-- USAGE EXAMPLES -->
## Tutorial

The using of ZHMFEM is introduced through modeling, calculation and output results of a practical example question.
This example question comes from 
__[Fundamentals of Finite Element Analysis](http://www.caemotion.com/pdf/Fundamentals%20of%20Finite%20Element%20Analysis.pdf)__
(Page105) written by Zeng Pan. 

(in Chinese: [曾攀](https://zh.wikipedia.org/wiki/%E6%9B%BE%E6%94%80) (清华大学)《有限元基础教程》北京：高等教育出版社，2009.7（2012.11重印）).

If you are interested in learning FEM, here is a great 
[lessen](https://www.bilibili.com/video/BV1iP4y1y7qh/?spm_id_from=333.337.search-card.all.click) 
by Prof. Zeng Pan.

The sample question is:
For the plane stress problem shown in the figure below, the given material parameters are:
E = 1 (Young's modulus), $\nu$ = 0.25 (Poisson's ratio). The thickness of the sheet is $1$. 

+ Known displacement boundary conditions: 

  $$u_A=0,\quad v_A=0,\quad u_D=0$$

+ Known external loads: 
  
  $$p_{Bx}=-1,\quad p_{By}=0,\quad p_{Cx}=1,\quad p_{Cy}=0,\quad p_{Dy}=0$$


Problems to be solved:

The linear triangular element (CST) and rectangular element are used to solve this problem respectively. The displacement, strain, stress, support reaction force at the nodes; strain energy, external force work and total potential energy of the system need to be calculated. Through the calculation results of the two schemes, the difference in calculation accuracy of the two elements for the same problem is compared.

Triangular elements and rectangular elements are used to mesh square thin plates, as shown in Figure 4-8(a) and 4-8(b).

![Plane stress problem 平面应力问题](imgs/CST.png "plane stress problem")

### Solution

1. **Solution of linear triangular element**

You can find the tri2d3n.rs file under 'zhmfem/examples/', 
which shows all the steps of setting material and size parameters, 
constructing nodes and cells, setting boundary conditions and external loads, calculating problems, output result files, etc.

```
#![allow(incomplete_features)]
#![feature(generic_const_exprs)]

use std::time::Instant;

use zhmfem::*;

fn main() {
    // set timing start
    let time_start = Instant::now();

    // -------- Part 0: Set initial parameters --------
    const E: Dtype = 0.0; // Exponent in scientific notation to base 10
    const CPU_CORES: usize = 2;

    let calc_method: &str = "lu"; // "lu" for LU decomposition algorithm or "gs" for gauss-seidel iteration method
    let calc_accuracy: Dtype = 0.001; // Calculation accuracy of iterative algorithm

    let parallel_or_singllel: &str = "s"; // "s" or "singllel" or "p" or "parallel"

    let thick: Dtype = 1.0; //Thickness of the plate
    let material: (Dtype, Dtype) = (1.0, 0.25); //Young's modulud & Poisson's ratio

    // -------- Part 1:  Meshing and applying boundary conditions --------
    // Set mesh and freedom parameters
    const R: usize = 2; // rows of nodes
    const C: usize = 2; // columns of nodes
    const M: usize = 3; // num of nodes in single element
    const F: usize = 2; // num of degree freedom at single node

    // Manually set coords and grouped nodes index
    /*
    let points: Vec<Vec<Dtype>> = vec![
        vec![0.0, 0.0],
        vec![1.0, 0.0],
        vec![1.0, 1.0],
        vec![0.0, 1.0],
    ];
    let grpdnidx: Vec<Vec<usize>> = vec![vec![0, 1, 3], vec![2, 3, 1]];

    // Set boundary conditions and external loads manually
    let zero_disp_index: Vec<usize> = vec![0, 1, 6];
    let force_index: Vec<usize> = vec![2, 4];
    let force_value: Vec<Dtype> = vec![-1.0, 1.0];
    let force_data: HashMap<usize, Dtype> = force_index
        .into_iter()
        .zip(force_value.into_iter())
        .collect();
    */

    // Automatically set coords and grouped nodes index
    // Auto-mesh generate coords and grouped nodes index
    const W: Dtype = 1.0; // width
    const H: Dtype = 1.0; // height
    let solid1 = Rectangle::new([0.0 as Dtype, 0.0 as Dtype], [W, H]);
    let (points, grpdnidx) = solid1.mesh_with_tri2d3n(R, C);

    // Set boundary conditions and external loads automatically
    let zero_disp_index: Vec<usize> = vec![0, 1, C * (R - 1) * F];
    let force_index: Vec<usize> = vec![(C - 1) * F, (C * R - 1) * F];
    let force_value: Vec<Dtype> = vec![-1.0, 1.0];

    let force_data: HashMap<usize, Dtype> = force_index
        .into_iter()
        .zip(force_value.into_iter())
        .collect();

    // -------- Part 2:  Construct nodes, elements and parts --------
    // Construct 2D nodes vector
    let nodes = nodes2d_vec(&points, &force_data);

    // Construct Tri2D3N elements vector
    let mut triangles = tri2d3n_vec(thick, &nodes, &grpdnidx, &material);

    // Construct 2D part & assembly global stiffness matrix
    let mut part: Part2D<'_, Tri2D3N<'_>, { R * C }, F, M> =
        Part2D::new(1, &nodes, &mut triangles, &grpdnidx);
    part.k_printer(parallel_or_singllel, CPU_CORES, E);

    // -------- Part 3:  Solve the problem --------
    // construct solver and solve the case
    let mut eqs: LinearEqs<{ R * C * F }> = LinearEqs::new(
        part.nodes_displacement(),
        part.nodes_force(),
        zero_disp_index,
        *part.k(parallel_or_singllel, CPU_CORES),
    );

    // 1) solve the linear equations of static system using direct method.
    // eqs.lu_direct_solver(); //LU decomposition method
    // let output_file = "LU.txt";

    // 2) solve the linear equations of static system using iter method.
    // eqs.gauss_seidel_iter_solver(0.001);
    // let output_file = "G-S.txt";

    // 3) or you can solve the problem with a more concise call:
    eqs.solve(calc_method, calc_accuracy);

    let calc_time: std::time::Duration = eqs.solver_time_consuming.unwrap();

    // write the displacement and force result into Node2D's field
    part.write_result(&eqs);

    // -------- Part 4:  Print all kinds of result --------
    print_1darr("qe", &part.nodes_displacement(), E, "v");
    print_1darr("fe", &part.nodes_force(), E, "v");

    println!("\n>>> System energy:");
    let strain_energy: Dtype = strain_energy(
        *part.k(parallel_or_singllel, CPU_CORES),
        part.nodes_displacement(),
    );
    let external_force_work: Dtype =
        external_force_work(part.nodes_force(), part.nodes_displacement());
    let potential_energy: Dtype = potential_energy(
        *part.k(parallel_or_singllel, CPU_CORES),
        part.nodes_force(),
        part.nodes_displacement(),
    );
    println!("\tE_d: {:-9.6} (deform energy)", strain_energy);
    println!("\tW_f: {:-9.6} (exforce works)", external_force_work);
    println!("\tE_p: {:-9.6} (potential energy)", potential_energy);

    part.elems
        .iter()
        .map(|elem| {
            println!("{}", elem.info(0.0));
        })
        .count();

    // -------- Part 5:  Write clac result into txt file --------
    let problem_type = "stress2D";
    let element_type = "Tri2D3N";
    let output_path = "/path/to/output/result/";
    let output_txt = format!(
        "{output_path}{problem_type}_{element_type}_{calc_method}_{parallel_or_singllel}.txt"
    );
    let output_vtk = format!(
        "{output_path}{problem_type}_{element_type}_{calc_method}_{parallel_or_singllel}.vtk"
    );

    // Output Calculation result into txt file
    part.txt_writer(
        &output_txt,
        calc_time,
        E,
        (strain_energy, external_force_work, potential_energy),
    )
    .expect(">>> !!! Failed to output text result file !!!");

    // Output Calculation result into vtk file
    part.vtk_writer(&output_vtk, element_type)
        .expect(">>> !!! Failed to output vtk file!");

    let total_time = time_start.elapsed();
    println!("\n>>> Total time consuming: {:?}", total_time);
}
```
    
  The rect2d4n.rs file under 'zhmfem/examples/' shows the steps to solve the example problem with linear rectangular elements. 
  Open the file to explore for more details.

2. **Build & Run**

Open Shell in the zhmfem root directory and use the following command
to compile and run the tri2d3n.rs file and rect2d4n.rs file:
```
cargo run -j 8 --release --example tri2d3n
```

```
cargo run -j 8 --release --example rect2d4n
```

3. **Check the results**

After compiling and running the tri2d3n.rs file, the following results are displayed in the Shell:

```
>>> Assembling Part2D(#1)'s stiffness matrix K1 in single thread ......

>>> Calculating Tri2D3N(#0)'s local stiffness matrix k0 ......

>>> Calculating Tri2D3N(#1)'s local stiffness matrix k1 ......

>>> Assembly Down!
        time consuming: 6.693µs

Part #1  K =  (* 10^0)
[[      0.733333       0.333333      -0.533333      -0.200000      -0.200000      -0.133333       0.000000       0.000000 ]
 [      0.333333       0.733333      -0.133333      -0.200000      -0.200000      -0.533333       0.000000       0.000000 ]
 [     -0.533333      -0.133333       0.733333       0.000000       0.000000       0.333333      -0.200000      -0.200000 ]
 [     -0.200000      -0.200000       0.000000       0.733333       0.333333       0.000000      -0.133333      -0.533333 ]
 [     -0.200000      -0.200000       0.000000       0.333333       0.733333       0.000000      -0.533333      -0.133333 ]
 [     -0.133333      -0.533333       0.333333       0.000000       0.000000       0.733333      -0.200000      -0.200000 ]
 [      0.000000       0.000000      -0.200000      -0.133333      -0.533333      -0.200000       0.733333       0.333333 ]
 [      0.000000       0.000000      -0.200000      -0.533333      -0.133333      -0.200000       0.333333       0.733333 ]]


>>> LU decomposition method down!
        time consuming = 18.356µs

qe = (10^0 *)
[[
           0.000000 
           0.000000 
          -1.718750 
          -0.937500 
           0.000000 
           0.781250 
           1.718750 
          -1.718750 
]]


fe = (10^0 *)
[[
           1.000000 
           0.000000 
          -1.000000 
           0.000000 
          -1.000000 
           0.000000 
           1.000000 
           0.000000 
]]


>>> System energy:
        E_d:  1.718750 (deform energy)
        W_f:  3.437500 (exforce works)
        E_p: -1.718750 (potential energy)

-----------------------------------------------------------------------------
Elem_Tri2D3N:
        Id:     0
        Area:     0.500000
        Mats:     1.000000 (Young's modulus)
                  0.250000 (Poisson's ratio)
        Nodes:
                Node_2D:
                        Id: 0
                        Coords: [    0.0000,     0.0000]
                        Displs: [    0.0000,     0.0000]
                        Forces: [    1.0000,     0.0000]

                Node_2D:
                        Id: 1
                        Coords: [    1.0000,     0.0000]
                        Displs: [   -1.7188,    -0.9375]
                        Forces: [   -1.0000,     0.0000]

                Node_2D:
                        Id: 2
                        Coords: [    0.0000,     1.0000]
                        Displs: [    0.0000,     0.7812]
                        Forces: [   -1.0000,     0.0000]

        Strain:
                [   -1.718750,     0.781250,    -0.937500]
        Stress:
                [   -1.625000,     0.375000,    -0.375000]

        Stiffness Matrix K0 =  (*10^0)
[[     0.733333      0.333333     -0.533333     -0.200000     -0.200000     -0.133333 ]
 [     0.333333      0.733333     -0.133333     -0.200000     -0.200000     -0.533333 ]
 [    -0.533333     -0.133333      0.533333      0.000000      0.000000      0.133333 ]
 [    -0.200000     -0.200000      0.000000      0.200000      0.200000      0.000000 ]
 [    -0.200000     -0.200000      0.000000      0.200000      0.200000      0.000000 ]
 [    -0.133333     -0.533333      0.133333      0.000000      0.000000      0.533333 ]]

-----------------------------------------------------------------------------
Elem_Tri2D3N:
        Id:     1
        Area:     0.500000
        Mats:     1.000000 (Young's modulus)
                  0.250000 (Poisson's ratio)
        Nodes:
                Node_2D:
                        Id: 3
                        Coords: [    1.0000,     1.0000]
                        Displs: [    1.7188,    -1.7188]
                        Forces: [    1.0000,     0.0000]

                Node_2D:
                        Id: 2
                        Coords: [    0.0000,     1.0000]
                        Displs: [    0.0000,     0.7812]
                        Forces: [   -1.0000,     0.0000]

                Node_2D:
                        Id: 1
                        Coords: [    1.0000,     0.0000]
                        Displs: [   -1.7188,    -0.9375]
                        Forces: [   -1.0000,     0.0000]

        Strain:
                [    1.718750,    -0.781250,     0.937500]
        Stress:
                [    1.625000,    -0.375000,     0.375000]

        Stiffness Matrix K1 =  (*10^0)
[[     0.733333      0.333333     -0.533333     -0.200000     -0.200000     -0.133333 ]
 [     0.333333      0.733333     -0.133333     -0.200000     -0.200000     -0.533333 ]
 [    -0.533333     -0.133333      0.533333      0.000000      0.000000      0.133333 ]
 [    -0.200000     -0.200000      0.000000      0.200000      0.200000      0.000000 ]
 [    -0.200000     -0.200000      0.000000      0.200000      0.200000      0.000000 ]
 [    -0.133333     -0.533333      0.133333      0.000000      0.000000      0.533333 ]]

>>> Writing calc results into txt file ......
    Down!

>>> Writing calc results into vtk file ......
    Down!

>>> Total time consuming: 771.936µs

```

After compiling and running the rect2d4n.rs file, the following results are displayed in the Shell:
```
>>> Assembling Part2D(#1)'s global stiffness matrix K1 ......

>>> Calculating Quad2D4N(#0)'s stiffness matrix k0 ......

Part #1  K =  (* 10^0)
[[     0.488889      0.166667     -0.288889     -0.033333      0.044444      0.033333     -0.244444     -0.166667 ]
[     0.166667      0.488889      0.033333      0.044444     -0.033333     -0.288889     -0.166667     -0.244444 ]
[    -0.288889      0.033333      0.488889     -0.166667     -0.244444      0.166667      0.044444     -0.033333 ]
[    -0.033333      0.044444     -0.166667      0.488889      0.166667     -0.244444      0.033333     -0.288889 ]
[     0.044444     -0.033333     -0.244444      0.166667      0.488889     -0.166667     -0.288889      0.033333 ]
[     0.033333     -0.288889      0.166667     -0.244444     -0.166667      0.488889     -0.033333      0.044444 ]
[    -0.244444     -0.166667      0.044444      0.033333     -0.288889     -0.033333      0.488889      0.166667 ]
[    -0.166667     -0.244444     -0.033333     -0.288889      0.033333      0.044444      0.166667      0.488889 ]]


>>> LU decomposition method down!
        time consuming = 2.172µs

qe = (10^0 *)
[[      0.000000       0.000000      -4.090909      -4.090909       0.000000      -0.000001       4.090909      -4.090909 ]]


fe = (10^0 *)
[[      1.000000       0.000000      -1.000000       0.000000      -1.000000      -0.000000       1.000000       0.000000 ]]


>>> System energy:
        E_d:  4.090909 (deform energy)
        W_f:  8.181817 (exforce works)
        E_p: -4.090909 (potential energy)
Quad2D4N k0 =  (* 10^0)
[[     0.488889     0.166667    -0.288889    -0.033333    -0.244444    -0.166667     0.044444     0.033333]
[     0.166667     0.488889     0.033333     0.044444    -0.166667    -0.244444    -0.033333    -0.288889]
[    -0.288889     0.033333     0.488889    -0.166667     0.044444    -0.033333    -0.244444     0.166667]
[    -0.033333     0.044444    -0.166667     0.488889     0.033333    -0.288889     0.166667    -0.244444]
[    -0.244444    -0.166667     0.044444     0.033333     0.488889     0.166667    -0.288889    -0.033333]
[    -0.166667    -0.244444    -0.033333    -0.288889     0.166667     0.488889     0.033333     0.044444]
[     0.044444    -0.033333    -0.244444     0.166667    -0.288889     0.033333     0.488889    -0.166667]
[     0.033333    -0.288889     0.166667    -0.244444    -0.033333     0.044444    -0.166667     0.488889]]


elem[0] strain:
        E_xx =        -4.090909
        E_yy =        -0.000001
        E_xy =        -4.090909

elem[0] stress:
        S_xx =        -4.363636
        S_yy =        -1.090910
        S_xy =        -1.636364

>>> Writing calc results into txt file ......
    Down!

>>> Total time consuming: 446.086µs
```
   
The calculation results of two kinds of elements for the same problem shows that:

1) The calculation result of linear quadrilateral element is higher than that of linear triangular element;

2) The calculation results of ZHMFEM are accurate and reliable.

There are several [examples](https://github.com/hunter1992/zhmfem/tree/main/examples) under 
zhmfem/exampls path.

These calculation examples show how to solve plane stress problems with various elements of ZHMFEM. 
The code for modeling, imposing boundary conditions, calculating systems of equations, 
and outputting the resulting files are clearly marked with comments in each sample file.


<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- ROADMAP -->
## Roadmap

- [✔] Add output .vtk file
- [ ] Add Tri2D6N element
- [ ] Add non-linear rod element
- [ ] Add Tri2D6N elements
- [ ] Add document
- [ ] Add non-linear FEM
    - [ ] plastic
    - [ ] geo-nonlinear

See the [open issues](https://github.com/hunter1992/zhmfem/issues) for a full list of proposed features (and known issues).

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- CONTRIBUTING -->
## Contributing

If you have a suggestion that would make this better, please fork the repo and create a pull request. 
You can also simply open an issue with the tag "enhancement".
Don't forget to give the project a star! Thanks again!

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/zhmfem`)
3. Commit your Changes (`git commit -m 'Add some zhmfem'`)
4. Push to the Branch (`git push origin feature/zhmfem`)
5. Open a Pull Request

### Top contributors:

<a href="https://github.com/othneildrew/Best-README-Template/graphs/contributors">
  <img src="https://contrib.rocks/image?repo=hunter1992/zhmfem" alt="contrib.rocks image" />
</a>

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- LICENSE -->
## License

Distributed under the GPL License. See [LICENSE](https://github.com/hunter1992/zhmfem/blob/main/LICENSE) for more information.

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- CONTACT -->
## Contact

[张洪铭](http://zhanghongming.com/) - [Zhang Hongming](https://github.com/hunter1992) - 1101510340zhm@gmail.com


<p align="right">(<a href="#readme-top">back to top</a>)</p>


<!-- ACKNOWLEDGMENTS -->
## Acknowledgments

* [追忆曾攀教授](https://www.tsinghua.org.cn/info/1951/37828.htm)
* [Choose an Open Source License](https://choosealicense.com)
* [GitHub Emoji Cheat Sheet](https://www.webpagefx.com/tools/emoji-cheat-sheet)
* [Malven's Flexbox Cheatsheet](https://flexbox.malven.co/)
* [Malven's Grid Cheatsheet](https://grid.malven.co/)
* [Img Shields](https://shields.io)
* [GitHub Pages](https://pages.github.com)
* [Font Awesome](https://fontawesome.com)
* [React Icons](https://react-icons.github.io/react-icons/search)

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://github.com/hunter1992/zhmfem -->
[contributors-shield]: https://img.shields.io/github/contributors/hunter1992/zhmfem.svg?style=for-the-badge
[contributors-url]: https://githubf.com/hunter1992/zhmfem/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/hunter1992/zhmfem.svg?style=for-the-badge
[forks-url]: https://github.com/hunter1992/zhmfem/network/members
[stars-shield]: https://img.shields.io/github/stars/hunter1992/zhmfem.svg?style=for-the-badge
[stars-url]: https://github.com/hunter1992/zhmfem/stargazers
[issues-shield]: https://img.shields.io/github/issues/hunter1992/zhmfem.svg?style=for-the-badge
[issues-url]: https://github.com/hunter1992/zhmfem/issues
[license-shield]: https://img.shields.io/github/license/hunter1992/zhmfem.svg?style=for-the-badge
[license-url]: https://github.com/hunter1992/zhmfem/blob/master/LICENSE.txt
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: http://zhanghongming.com
[product-screenshot]: images/screenshot.png
[Rustc]: https://img.shields.io/badge/Rust-000000?style=flat&logo=rust&logoColor=white
[Rust-url]: https://www.rust-lang.org/
