use extendr_api::prelude::*;
use std::f64;

/// Compute Jensen-Shannon divergence between two probability distributions.
/// @export
#[extendr]
fn jsd_r(p: RMatrix<f64>, q: RMatrix<f64>) -> f64 {
    // Work directly with slices - no copying
    let p_slice = p.as_real_slice().unwrap();
    let q_slice = q.as_real_slice().unwrap();

    let mut kl_pm = 0.0;
    let mut kl_qm = 0.0;

    // Single pass calculation to avoid multiple iterations
    for (pi, qi) in p_slice.iter().zip(q_slice.iter()) {
        let mi = 0.5 * (pi + qi);
        
        if *pi > 0.0 && mi > 0.0 {
            kl_pm += pi * (pi / mi).log2();
        }
        
        if *qi > 0.0 && mi > 0.0 {
            kl_qm += qi * (qi / mi).log2();
        }
    }

    0.5 * (kl_pm + kl_qm)
}

/// Compute Jensen-Shannon divergence matrix for a given input matrix.
/// @export
#[extendr]
fn dist_jsd_r(in_matrix: RMatrix<f64>) -> RMatrix<f64> {
    let rows = in_matrix.nrows();
    let cols = in_matrix.ncols();
    let mut result = RMatrix::new_matrix(rows, rows, |_, _| 0.0);
    
    // Get the underlying data slice (column-major order)
    let data = in_matrix.as_real_slice().unwrap();
    
    for i in 0..rows {
        for j in 0..i {
            let jsd_value = jsd_rows_from_slice(data, i, j, rows, cols);
            result[[i, j]] = jsd_value;
            result[[j, i]] = jsd_value;
        }
    }
    result
}


// Helper function for Method 2: JSD directly from slice data
fn jsd_rows_from_slice(data: &[f64], row_i: usize, row_j: usize, rows: usize, cols: usize) -> f64 {
    let mut kl_pm = 0.0;
    let mut kl_qm = 0.0;

    // In column-major order: element at (row, col) is at index: row + col * rows
    for col in 0..cols {
        let pi = data[row_i + col * rows];
        let qi = data[row_j + col * rows];
        let mi = 0.5 * (pi + qi);
        
        if pi > 0.0 && mi > 0.0 {
            kl_pm += pi * (pi / mi).log2();
        }
        
        if qi > 0.0 && mi > 0.0 {
            kl_qm += qi * (qi / mi).log2();
        }
    }

    0.5 * (kl_pm + kl_qm)
}

/// Compute Jensen-Shannon divergence matrix for two input matrices.
/// @export
#[extendr]
fn dist_jsd2_r(in_matrix_tr: RMatrix<f64>, in_matrix_co: RMatrix<f64>) -> RMatrix<f64> {
    let rows_tr = in_matrix_tr.nrows();
    let rows_co = in_matrix_co.nrows();
    let cols_tr = in_matrix_tr.ncols();
    let cols_co = in_matrix_co.ncols();
    
    // Initialize result matrix
    let mut result = RMatrix::new_matrix(rows_tr, rows_co, |_, _| 0.0);
    
    for i in 0..rows_tr {
        for j in 0..rows_co {
            // Extract row i from tr matrix and row j from co matrix
            let row_tr: Vec<f64> = (0..cols_tr).map(|c| in_matrix_tr[[i, c]]).collect();
            let row_co: Vec<f64> = (0..cols_co).map(|c| in_matrix_co[[j, c]]).collect();
            
            // Calculate JSD between the two rows
            let jsd_value = jsd_vectors(&row_tr, &row_co);
            result[[i, j]] = jsd_value;  // Note: [[i, j]] not [i][j]
        }
    }
    result
}

// Helper function for vector-based JSD calculation
fn jsd_vectors(p: &[f64], q: &[f64]) -> f64 {
    let mut kl_pm = 0.0;
    let mut kl_qm = 0.0;
    let len = p.len().min(q.len());  // Handle different lengths

    for i in 0..len {
        let pi = p[i];
        let qi = q[i];
        let mi = 0.5 * (pi + qi);
        
        if pi > 0.0 && mi > 0.0 {
            kl_pm += pi * (pi / mi).log2();
        }
        
        if qi > 0.0 && mi > 0.0 {
            kl_qm += qi * (qi / mi).log2();
        }
    }

    0.5 * (kl_pm + kl_qm)
}



// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod speed;
    // fn hello_world;
    fn jsd_r;
    fn dist_jsd_r;
    fn dist_jsd2_r;
}
