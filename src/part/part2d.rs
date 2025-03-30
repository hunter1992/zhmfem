use crate::{ADtype, Dtype, Export, LinearEqs, Node2D, K};
use std::fmt::Write as _;
use std::io::{BufWriter, Write};
use std::sync::atomic::Ordering;
use std::sync::Arc;
use std::thread;
use std::time::Instant;

/// Three generic const: N for N_NODE, F for N_FREEDOM, M for N_NODE in single element
pub struct Part2D<'part2d, Elem: K, const N: usize, const F: usize, const M: usize>
where
    [[Dtype; N * F]; N * F]: Sized,
{
    pub id: usize,
    pub nodes: &'part2d [Node2D],
    pub elems: &'part2d mut [Elem],
    pub cplds: &'part2d [Vec<usize>], //cplds: coupled nodes index
    pub k_matrix: Option<[[Dtype; N * F]; N * F]>,
    pub assembly_time_consuming: Option<std::time::Duration>,
}

impl<'part2d, Elem: K, const N: usize, const F: usize, const M: usize>
    Part2D<'part2d, Elem, N, F, M>
where
    [[Dtype; N * F]; N * F]: Sized,
{
    /// Construct a Part2D
    /// Args: nodes --- Vec ptr to all 2D nodes
    ///       elems --- Slice of all element stiffness matrix pointers
    ///       cplds --- Slice of the ID combinations of Node2D within a single cell
    pub fn new(
        id: usize,
        nodes: &'part2d [Node2D],
        elems: &'part2d mut [Elem],
        cplds: &'part2d [Vec<usize>],
    ) -> Self {
        Part2D {
            id,
            nodes,
            elems,
            cplds,
            k_matrix: None,
            assembly_time_consuming: None,
        }
    }

    /// Get displacement of all nodes
    pub fn nodes_displacement(&self) -> [Dtype; N * F] {
        let mut data: [Dtype; N * F] = [0.0; N * F];
        for (idx, node) in self.nodes.iter().enumerate() {
            data[idx * 2] = node.displs.borrow()[0];
            data[idx * 2 + 1] = node.displs.borrow()[1];
        }
        data
    }

    /// Get force of all nodes
    pub fn nodes_force(&self) -> [Dtype; N * F] {
        let mut data: [Dtype; N * F] = [0.0; N * F];
        for (idx, node) in self.nodes.iter().enumerate() {
            data[idx * 2] = node.forces.borrow()[0];
            data[idx * 2 + 1] = node.forces.borrow()[1];
        }
        data
    }

    /// Calculate part's global stiffness matrix
    /// The value of arg parallel_or_singllel can be:
    /// "parallel" or "p" or "singllel" or "s"
    pub fn k(&mut self, parallel_or_singllel: &str, cpu_cores: usize) -> &[[Dtype; N * F]; N * F]
    where
        <Elem as K>::Kmatrix: std::ops::Index<usize, Output = [Dtype; M * F]> + Sync,
    {
        if self.k_matrix.is_none() {
            if self.cplds.len() != self.elems.len() {
                println!("\n---> Error! From Part2D.k_singllel func.");
                println!("     The count of elements not eq to K mat size.");
                panic!("---> Assembly global K failed!");
            }

            match parallel_or_singllel{
            "s"=> self.k_singllel(cpu_cores),
            "p"=> self.k_parallel(cpu_cores),
            "singllel"=> self.k_singllel(cpu_cores),
            "parallel"=> self.k_parallel(cpu_cores),
            _ => panic!("!!! Wrong arg for Part's stiffness matrix calc fun k(), from/src/part/part2d.rs fn k()"),
            }
        } else {
            self.k_matrix.as_ref().unwrap()
        }
    }

    fn k_singllel(&mut self, _cpu_cores: usize) -> &[[Dtype; N * F]; N * F]
    where
        <Elem as K>::Kmatrix: std::ops::Index<usize, Output = [Dtype; M * F]>,
    {
        println!(
            "\n>>> Assembling Part2D(#{})'s stiffness matrix K{} in single thread ......",
            self.id, self.id
        );

        let elems: Vec<&<Elem as K>::Kmatrix> =
            self.elems.iter_mut().map(|elem| elem.k()).collect();

        let mut part_stiffness_mat: [[Dtype; N * F]; N * F] = [[0.0; N * F]; N * F];

        let timing_start = Instant::now();

        self.cplds
            .iter()
            .zip(elems.iter())
            .map(|(cpld, &elem)| {
                for idx in 0..M {
                    for idy in 0..M {
                        for row in 0..F {
                            for col in 0..F {
                                part_stiffness_mat[cpld[idx] * F + row][cpld[idy] * F + col] +=
                                    elem[idx * F + row][idy * F + col];
                            }
                        }
                    }
                }
            })
            .count();

        let time_consuming = timing_start.elapsed();
        println!(
            "\n>>> Assembly Down!\n        time consuming: {:?}",
            time_consuming
        );
        self.assembly_time_consuming = Some(time_consuming);
        self.k_matrix.get_or_insert(part_stiffness_mat)
    }

    fn k_parallel(&mut self, cpu_cores: usize) -> &[[Dtype; N * F]; N * F]
    where
        <Elem as K>::Kmatrix: std::ops::Index<usize, Output = [Dtype; M * F]> + Sync,
    {
        println!(
            "\n>>> Assembling Part2D(#{})'s stiffness matrix K{} in {} threads ......",
            self.id, self.id, cpu_cores
        );

        let elems: Vec<&<Elem as K>::Kmatrix> =
            self.elems.iter_mut().map(|elem| elem.k()).collect();
        let mut part_stiffness_mat: [[Dtype; N * F]; N * F] = [[0.0; N * F]; N * F];
        let mat: Arc<Vec<Vec<ADtype>>> = Arc::new(
            (0..F * N)
                .map(|_| (0..F * N).map(|_| ADtype::new(0.0)).collect())
                .collect(),
        );

        let density: usize = (((self.cplds.len() as f64) / (cpu_cores as f64)).ceil()) as usize;
        let elemss: Vec<_> = elems.chunks(density).collect();
        let cpldss: Vec<_> = self.cplds.chunks(density).collect();

        let timing_start = Instant::now();

        for n in 0..cpu_cores {
            let elems = elemss[n];
            let cplds = cpldss[n];
            let mat_clone = Arc::clone(&mat);
            thread::scope(|s| {
                s.spawn(|| {
                    cplds
                        .iter()
                        .zip(elems.iter())
                        .map(|(cpld, &elem)| {
                            for idx in 0..M {
                                for idy in 0..M {
                                    for row in 0..F {
                                        for col in 0..F {
                                            mat_clone[cpld[idx] * F + row][cpld[idy] * F + col]
                                                .fetch_add(
                                                    elem[idx * F + row][idy * F + col],
                                                    Ordering::Relaxed,
                                                );
                                        }
                                    }
                                }
                            }
                        })
                        .count();
                });
            });
        }

        let time_consuming = timing_start.elapsed();
        println!(
            "\n>>> Assembly Down!\n        time consuming: {:?}",
            time_consuming
        );
        self.assembly_time_consuming = Some(time_consuming);
        for x in 0..F * N {
            for y in 0..F * N {
                part_stiffness_mat[x][y] = mat[x][y].load(Ordering::Relaxed);
            }
        }
        self.k_matrix.get_or_insert(part_stiffness_mat)
    }

    /// Print part's global stiffness matrix
    pub fn k_printer(&mut self, parallel_or_singllel: &str, cpu_cores: usize, n_exp: Dtype)
    where
        <Elem as K>::Kmatrix: std::ops::Index<usize, Output = [Dtype; M * F]> + Sync,
    {
        if self.k_matrix.is_none() {
            self.k(parallel_or_singllel, cpu_cores);
        }

        print!("\nPart #{}  K =  (* 10^{})\n[", self.id, n_exp as u8);
        for row in 0..(N * F) {
            if row == 0 {
                print!("[");
            } else {
                print!(" [");
            }
            for col in 0..(N * F) {
                print!(
                    " {:>-13.6} ",
                    self.k_matrix.unwrap()[row][col] / (10.0_f64.powf(n_exp as f64)) as Dtype
                );
            }
            if row == (N * F - 1) {
                println!("]]");
            } else {
                println!("]");
            }
        }
        println!("");
    }

    /// Return Part2D's stiffness matrix's formeted string
    pub fn k_string(&self, n_exp: Dtype) -> String {
        let mut k_matrix = String::new();
        write!(
            k_matrix,
            "\nPart #{}  K =  (* 10^{})\n",
            self.id, n_exp as u8
        )
        .expect("Write Part2D Stiffness Matrix Failed!");
        for row in 0..(N * F) {
            if row == 0 {
                write!(k_matrix, "[[").expect("Write Part2D Stiffness Matrix Failed!");
            } else {
                write!(k_matrix, " [").expect("Write Part2D Stiffness Matrix Failed!");
            }
            for col in 0..(N * F) {
                write!(
                    k_matrix,
                    " {:>-13.6} ",
                    self.k_matrix.unwrap()[row][col] / (10.0_f64.powf(n_exp as f64)) as Dtype
                )
                .expect("Write Part2D Stiffness Matrix Failed!");
            }
            if row == { N * F - 1 } {
                write!(k_matrix, "]]").expect("Write Part2D Stiffness Matrix Failed!");
            } else {
                write!(k_matrix, " ]\n").expect("Write Part2D Stiffness Matrix Failed!");
            }
        }
        write!(k_matrix, "\n").expect("Write Part2D Stiffness Matrix Failed!");
        k_matrix
    }

    /// Write the disp and force result into nodes
    pub fn write_result(&mut self, slv: &LinearEqs<{ N * F }>) {
        let disp = slv.disps;
        let force = slv.external_force.unwrap();
        for (idx, node) in self.nodes.iter().enumerate() {
            node.displs.borrow_mut()[0] = disp[idx * 2];
            node.displs.borrow_mut()[1] = disp[idx * 2 + 1];
            node.forces.borrow_mut()[0] = force[idx * 2];
            node.forces.borrow_mut()[1] = force[idx * 2 + 1];
        }
    }
}

impl<'a, Elem: K, const N: usize, const F: usize, const M: usize> Export
    for Part2D<'a, Elem, N, F, M>
where
    [[Dtype; N * F]; N * F]: Sized,
{
    fn txt_writer(
        &self,
        target_file: &str,
        calc_time: std::time::Duration,
        n_exp: Dtype,
        energy: (Dtype, Dtype, Dtype),
    ) -> std::io::Result<bool> {
        let txt_file = std::fs::File::create(target_file).unwrap();
        let mut text_writer = BufWriter::new(txt_file);

        println!("\n>>> Writing calc results into txt file ......");
        write!(text_writer, ">>> ZHMFEM calculating results:").expect("Write txt file error!");

        let (deform, exwork, potential) = energy;
        write!(text_writer, "\n").expect("Write parts' result into yxy file failed!!!");
        write!(text_writer, "\n>>> System energy:\n")
            .expect("Write parts' result into yxy file failed!!!");
        write!(text_writer, "\t\t\tE_d: {:-13.6} (deform energy)\n", deform)
            .expect("Write parts' result into yxy file failed!!!");
        write!(text_writer, "\t\t\tW_f: {:-13.6} (exforce works)\n", exwork)
            .expect("Write parts' result into yxy file failed!!!");
        write!(
            text_writer,
            "\t\t\tE_p: {:-13.6} (potential energy)\n",
            potential
        )
        .expect("Write parts' result into yxy file failed!!!");

        write!(
            text_writer,
            "\n>>> Solver time consuming: {:?}\n",
            calc_time
        )
        .expect("Write txt file error!");

        write!(text_writer, "\n>>> Part2D stiffness matrix:\n").expect("Write info failed!");
        write!(text_writer, "{}", self.k_string(n_exp)).expect("Write info failed!");

        write!(text_writer, "\n>>> Details of each element:")
            .expect("Write parts' result into txt file failed!!!");
        for elem in self.elems.iter() {
            write!(text_writer, "{}\n", elem.info(n_exp)).expect("Write info failed!");
        }
        text_writer.flush().expect("!!! Flush txt file failed!");
        println!("    Down!");
        Ok(true)
    }

    fn vtk_writer(&self, _target_file: &str) -> std::io::Result<bool> {
        Ok(true)
    }
}
