use zhmfem::*;

fn main() {
    let matrix: [[Dtype; 8]; 8] = [
        [7., 0., 1., 0., 0., 2., 7., 0.],
        [0., -4., 8., 0., 2., 0., 0., 0.],
        [1., 8., 1., 0., 0., 0., 0., 5.],
        [0., 0., 0., 7., 0., 0., 9., 0.],
        [0., 2., 0., 0., 5., -1., 5., 0.],
        [2., 0., 0., 0., -1., 0., 0., 5.],
        [7., 0., 0., 9., 5., 0., 11., 0.],
        [0., 0., 5., 0., 0., 5., 0., 5.],
    ];

    println!(">>> The matrix:");
    print_2darr("A", 0, &matrix, 0.0);

    let sks = compress_symmetry_matrix_sks(&matrix);
    let recover_sks: [[Dtype; 8]; 8] = sks.recover_square_arr();
    print!(">>> SKS format:\n{:-5.2}", sks);
    print_2darr("Recover by SKS", 0, &recover_sks, 0.0);
    println!("=================================================================================");
    let mut csr = compress_symmetry_matrix_csr(&matrix);
    let recover_csr: [[Dtype; 8]; 8] = *csr.recover_square_arr();
    print!(">>> CSR format:\n{:-5.2}", csr);
    print_2darr("Recover by CSR", 0, &recover_csr, 0.0);
    println!("=================================================================================");
    csr.convert_1base();
    print!(">>> CSR_1based format:\n{:-5.2}", csr);
    let recover_csr_1base: [[Dtype; 8]; 8] = *csr.recover_square_arr();
    print_2darr("Recover by CSR_1based", 0, &recover_csr_1base, 0.0);

    let index = vec![2, 4, 5];
    let sub_matrix = sks.get_sub_matrix(&index);
    println!("=================================================================================");
    print_2darr("A", 0, &matrix, 0.0);
    println!("extract index: {:?}", index);
    print_2dvec("SubMat", &sub_matrix, 0.0);
}
