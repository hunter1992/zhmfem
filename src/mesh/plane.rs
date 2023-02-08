type Points = Vec<Vec<f64>>;
type Coupls = Vec<Vec<usize>>;

pub struct Rectangle {
    pub point1: [f64; 2],
    pub point2: [f64; 2],
}

impl Rectangle {
    pub fn new(point1: [f64; 2], point2: [f64; 2]) -> Self {
        Rectangle { point1, point2 }
    }

    pub fn mesh_with_tri(&self, rows: usize, cols: usize) -> (Points, Coupls) {
        let mut points: Vec<Vec<f64>> = vec![];
        let mut coupls: Vec<Vec<usize>> = vec![];

        let x0 = if self.point2[0] > self.point1[0] {
            self.point1[0]
        } else {
            self.point2[0]
        };
        let y0 = if self.point2[1] > self.point1[1] {
            self.point1[1]
        } else {
            self.point2[1]
        };
        let length = (self.point2[0] - self.point1[0]).abs();
        let height = (self.point2[1] - self.point1[1]).abs();
        let l_step = length / ((cols - 1) as f64);
        let h_step = height / ((rows - 1) as f64);

        // fill the nodes coords vector
        for r in 0..rows {
            for c in 0..cols {
                points.push(vec![x0 + c as f64 * l_step, y0 + r as f64 * h_step]);
            }
        }

        // generate the coupled nodes in single element
        /* the result of the mesh is similar to
         *    3—————4
         *    |\    |
         *    | \   |
         *    |  \  |
         *    |   \ |
         *    |    \|
         *    1—————2
         */
        for r in 0..(cols - 1) {
            for c in 0..(rows - 1) {
                coupls.push(vec![
                    r * cols + c + 1,
                    r * cols + c + 2,
                    (r + 1) * cols + c + 1,
                ]);
                coupls.push(vec![
                    (rows - r) * cols - c,
                    (rows - r) * cols - c - 1,
                    (rows - r - 1) * cols - c,
                ]);
            }
        }
        (points, coupls)
    }

    pub fn mesh_with_rect(&self, rows: usize, cols: usize) -> (Points, Coupls) {
        let mut points: Vec<Vec<f64>> = vec![];
        let mut coupls: Vec<Vec<usize>> = vec![];

        let x0 = if self.point2[0] > self.point1[0] {
            self.point1[0]
        } else {
            self.point2[0]
        };
        let y0 = if self.point2[1] > self.point1[1] {
            self.point1[1]
        } else {
            self.point2[1]
        };
        let length = (self.point2[0] - self.point1[0]).abs();
        let height = (self.point2[1] - self.point1[1]).abs();
        let l_step = length / ((cols - 1) as f64);
        let h_step = height / ((rows - 1) as f64);

        // generate the coupled nodes in single element
        /* the result of the mesh is similar to
         *    4—————3
         *    |     |
         *    |     |
         *    |     |
         *    1—————2
         */
        for r in 0..rows {
            for c in 0..cols {
                points.push(vec![x0 + c as f64 * l_step, y0 + r as f64 * h_step]);
            }
        }

        for r in 0..(cols - 1) {
            for c in 0..(rows - 1) {
                coupls.push(vec![
                    r * cols + c + 1,
                    r * cols + c + 2,
                    (r + 1) * cols + c + 2,
                    (r + 1) * cols + c + 1,
                ]);
                if rows < 3 || cols < 3 {
                    break;
                }
            }
        }
        (points, coupls)
    }
}
