# zhmfem
A finite element calculation crate based on Rust.

<!-- PROJECT SHIELDS -->

[![Contributors][contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![Issues][issues-shield]][issues-url]
[![MIT License][license-shield]][license-url]
[![LinkedIn][linkedin-shield]][linkedin-url]

<!-- links -->
[your-project-path]:https://github.com/hunter1992/zhmfem
[contributors-shield]: https://img.shields.io/github/contributors/hunter1992/zhmfem.svg?style=flat-square
[contributors-url]: https://github.com/hunter1992/zhmfem/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/hunter1992/zhmfem.svg?style=flat-square
[forks-url]: https://github.com/hunter1992/zhmfem/network/members
[stars-shield]: https://img.shields.io/github/stars/hunter1992/zhmfem.svg?style=flat-square
[stars-url]: https://github.com/hunter1992/zhmfem/stargazers
[issues-shield]: https://img.shields.io/github/issues/hunter1992/zhmfem.svg?style=flat-square
[issues-url]: https://img.shields.io/github/issues/hunter1992/zhmfem.svg
[license-shield]: https://img.shields.io/github/license/hunter1992/zhmfem.svg?style=flat-square
[license-url]: https://github.com/hunter1992/zhmfem/blob/master/LICENSE.txt
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=flat-square&logo=linkedin&colorB=555
[linkedin-url]: https://github.com/hunter1992

## 目录
- [上手指南](#上手指南)
    - [配置要求](#配置要求)
    - [安装步骤](#安装与使用)
    - [计算示例](#计算示例)
- [文件目录](#文件目录说明)
- [开发的架构](#开发的架构)
- [部署](#部署)
- [使用到的框架](#使用到的框架)
- [贡献者](#贡献者)
    - [如何参与开源项目](#如何参与开源项目)
- [版本控制](#版本控制)
- [作者](#作者)
- [鸣谢](#鸣谢)

## 上手指南



### **安装与使用**

1.Install Rust on your operating system, [install Rust now](https://www.rust-lang.org/tools/install)

2.Clone the repo:
```
git clone https://github.com/hunter1992/zhmfem.git
```

3.Check out the examples under 'zhmfem/examples/' path using:
```
cargo run --examples <example-name>
```
Or you can compile and run the examples more quickly using the following command:
```
cargo run -j 8 --release --example <example-name>
```

4.Use the example file as a template
to write a .rs file for the problem 
you need to solve (still placed in 
the examples path), compile and run it, 
and view the output results.

### 计算示例

zhmfem目录的example路径下有许多示例，这些例子使用各种单元求解不同类型的有限元问题，并对计算的结果进行了展示（包括将结果
展示在命令行或输出成供人类阅读的.txt文件或供ParaView显示的.vtk文件）。

借助以CST（线性三角形单元）求解平面应力问题的例子简单介绍zhmfem的使用流程。

##### Step0   从crate根中引入构造问题/计算/输出结果所需的各种类型、trait、函数

### 配置要求

#### 操作系统 

作者开发zhmfem使用的操作系统是Manjaro Linux，由于zhmfem目前处于早期核心功能开发阶段，
因此未对MacOS及Windows系统进行适配。截至2024年12月，作者所使用的操作系统信息如下：
```
 ██████████████████  ████████     zhm@zhm
 ██████████████████  ████████     OS: Manjaro 24.2.0 Yonada
 ██████████████████  ████████     Kernel: x86_64 Linux 6.1.119-1-MANJARO
 ██████████████████  ████████     Uptime: 12h 4m
 ████████            ████████     Packages: 1599
 ████████  ████████  ████████     Shell: zsh 5.9
 ████████  ████████  ████████     Resolution: 2560x1440
 ████████  ████████  ████████     DE: KDE
 ████████  ████████  ████████     WM: KWin
 ████████  ████████  ████████     GTK Theme: Breath [GTK2/3]
 ████████  ████████  ████████     Icon Theme: breeze
 ████████  ████████  ████████     Disk: 185G / 325G (61%)
 ████████  ████████  ████████     CPU: Intel Core i5-8265U @ 8x 3.9GHz [44.0°C]
 ████████  ████████  ████████     GPU: NVIDIA GeForce MX250
                                  RAM: 3456MiB / 7699MiB
```

#### Rust版本

zhmfem的开发在Rust的nightly版本下展开，作者使用的Rust版本信息如下：
```
rustc 1.85.0-nightly (7db7489f9 2024-11-25)
binary: rustc
commit-date: 2024-11-25
host: x86_64-unknown-linux-gnu
release: 1.85.0-nightly
LLVM version: 19.1.4
```

### 文件目录说明
eg:

```
filetree 
├── Cargo.lock
├── Cargo.toml
├── examples
│   ├── beam1d2n.rs
│   ├── rect2d4n.rs
│   ├── rod1d2n.rs
│   ├── rod2d_10bar.rs
│   ├── rod2d2n.rs
│   └── tri2d3n.rs
├── README.md
├── src
│   ├── calc.rs
│   ├── elem
│   │   ├── dim1
│   │   │   ├── beam.rs
│   │   │   ├── mod.rs
│   │   │   └── rod.rs
│   │   ├── dim2
│   │   │   ├── mod.rs
│   │   │   ├── quadrila.rs
│   │   │   ├── rod.rs
│   │   │   └── triangle.rs
│   │   └── mod.rs
│   ├── lib.rs
│   ├── main.rs
│   ├── mesh
│   │   ├── mod.rs
│   │   └── plane.rs
│   ├── node.rs
│   └── part
│       ├── mod.rs
│       ├── part1d.rs
│       └── part2d.rs
└── tree.txt
```

### 使用到的框架

- [nalgebra](https://github.com/dimforge/nalgebra)
- [vtkio](https://github.com/elrnv/vtkio)

#### 如何参与开源项目

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request


### 版本控制

该项目使用Git进行版本管理。您可以在repository参看当前可用版本。

### 作者

- [Zhang Hongming](https://github.com/hunter1992) e-mail:1101510340zhm@gmail.com

### 版权说明

该项目签署了MIT 授权许可，详情请参阅 [LICENSE.txt](https://github.com/shaojintian/Best_README_template/blob/master/LICENSE.txt)

### 鸣谢

- [GitHub Emoji Cheat Sheet](https://www.webpagefx.com/tools/emoji-cheat-sheet)
- [Img Shields](https://shields.io)
- [Choose an Open Source License](https://choosealicense.com)
- [GitHub Pages](https://pages.github.com)
- [Animate.css](https://daneden.github.io/animate.css)
- [xxxxxxxxxxxxxx](https://connoratherton.com/loaders)
