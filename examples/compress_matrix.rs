use zhmfem::*;

fn main() {
    let a: [[Dtype; 8]; 8] = [
        [1., 0., 0., 0., 0., 0., 0., 0.7471593047883867],
        [
            0.,
            1.,
            0.7051461125222229,
            0.,
            0.6295785540137455,
            0.,
            0.,
            0.,
        ],
        [
            0.,
            0.7051461125222229,
            1.4972310400052034,
            0.,
            0.443944869890155,
            0.,
            0.,
            0.,
        ],
        [0., 0., 0., 1., 0., 0., 0., 0.],
        [
            0.,
            0.6295785540137455,
            0.443944869890155,
            0.,
            1.3963691556740387,
            0.,
            0.,
            0.,
        ],
        [0., 0., 0., 0., 0., 1., 0., 0.],
        [0., 0., 0., 0., 0., 0., 1., 0.],
        [
            0.7471593047883867,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            1.5582470267318653,
        ],
    ];

    print_2darr("A", 0, &a, 0.0);

    let comp_sks_a = compress_matrix_sks(&a);
    let comp_csr0_a = compress_matrix_csr_0based(&a);
    let comp_csr1_a = compress_matrix_csr_1based(&a);

    let recover_sks: [[Dtype; 8]; 8] = *comp_sks_a.recover();
    let recover_csr_0: [[Dtype; 8]; 8] = *comp_csr0_a.recover();
    let recover_csr_1: [[Dtype; 8]; 8] = *comp_csr1_a.recover();

    print!(">>> SKS:\n{:-10.4}", comp_sks_a);
    print_2darr("Recover by SKS", 0, &recover_sks, 0.0);
    println!("-------------------------------");
    print!(">>> CSR_0based:\n{:-10.4}", comp_csr0_a);
    print_2darr("Recover by csr_0", 0, &recover_csr_0, 0.0);
    println!("-------------------------------");
    print!(">>> CSR_1based:\n{:-10.4}", comp_csr1_a);
    print_2darr("Recover by SKS", 0, &recover_csr_1, 0.0);
}
