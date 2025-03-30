use crate::Dtype;

type Points = Vec<Vec<Dtype>>;
type Coupls = Vec<Vec<usize>>;

pub struct Rectangle {
    pub point1: [Dtype; 2],
    pub point2: [Dtype; 2],
}

impl Rectangle {
    pub fn new(point1: [Dtype; 2], point2: [Dtype; 2]) -> Self {
        Rectangle { point1, point2 }
    }

    pub fn mesh_with_tri2d6n(&self, rows: usize, cols: usize) -> (Points, Coupls) {
        let mut coords: Vec<Vec<Dtype>> = vec![];
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
        let half_l_step = length / (2.0 * ((cols - 1) as Dtype));
        let half_h_step = height / (2.0 * ((rows - 1) as Dtype));

        // fill the nodes coords vector
        for r in 0..(2 * rows - 1) {
            for c in 0..(2 * cols - 1) {
                coords.push(vec![
                    x0 + c as Dtype * half_l_step,
                    y0 + r as Dtype * half_h_step,
                ]);
            }
        }

        // generate the coupled nodes in single element
        /* the result of the mesh is similar to
         *    6——7——8
         *    |\    |
         *    | \   |
         *    3  4  5
         *    |   \ |
         *    |    \|
         *    0——1——2
         *    tri2d6n cplds: [[0, 2, 6, 1, 4, 3], [8, 6, 2, 7, 4, 5]]
         */
        for r in 0..(cols - 1) {
            for c in 0..(rows - 1) {
                // 当cols=2, rows=2时,r和c都只等于一次0
                coupls.push(vec![
                    r * (4 * cols - 2) + c * 2,
                    r * (4 * cols - 2) + c * 2 + 2,
                    (r + 1) * (4 * cols - 2) + c * 2,
                    r * (4 * cols - 2) + c * 2 + 1,
                    r * (4 * cols - 2) + c * 2 + 1 + (2 * cols - 1),
                    r * (4 * cols - 2) + c * 2 + (2 * cols - 1),
                ]);
                coupls.push(vec![
                    (r + 1) * (4 * cols - 2) + c * 2 + 2,
                    (r + 1) * (4 * cols - 2) + c * 2,
                    r * (4 * cols - 2) + c * 2 + 2,
                    (r + 1) * (4 * cols - 2) + c * 2 + 1,
                    r * (4 * cols - 2) + c * 2 + 1 + (2 * cols - 1),
                    r * (4 * cols - 2) + c * 2 + 2 + (2 * cols - 1),
                ]);
            }
        }
        (coords, coupls)
    }

    pub fn mesh_with_tri2d3n(&self, rows: usize, cols: usize) -> (Points, Coupls) {
        let mut nodes: Vec<Vec<Dtype>> = vec![];
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
        let l_step = length / ((cols - 1) as Dtype);
        let h_step = height / ((rows - 1) as Dtype);

        // fill the nodes coords vector
        for r in 0..rows {
            for c in 0..cols {
                nodes.push(vec![x0 + c as Dtype * l_step, y0 + r as Dtype * h_step]);
            }
        }

        // generate the coupled nodes in single element
        /* the result of the mesh is similar to
         *    2—————3
         *    |\    |
         *    | \   |
         *    |  \  |
         *    |   \ |
         *    |    \|
         *    0—————1
         *    tri2d3n cplds: [[0, 1, 2], [3, 2, 1]]
         */
        for r in 0..(cols - 1) {
            for c in 0..(rows - 1) {
                // 当cols=2, rows=2时,r和c都只等于一次0
                coupls.push(vec![r * cols + c, r * cols + c + 1, (r + 1) * cols + c]);
                coupls.push(vec![
                    (rows - r) * cols - c - 1,
                    (rows - r) * cols - c - 2,
                    (rows - r - 1) * cols - c - 1,
                ]);
            }
        }
        (nodes, coupls)
    }

    pub fn mesh_with_rect(&self, rows: usize, cols: usize) -> (Points, Coupls) {
        let mut points: Vec<Vec<Dtype>> = vec![];
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
        let l_step = length / ((cols - 1) as Dtype);
        let h_step = height / ((rows - 1) as Dtype);

        // generate the coupled nodes in single element
        /* the result of the mesh is similar to
         *    3—————4
         *    |     |
         *    |     |
         *    |     |
         *    1—————2
         */
        for r in 0..rows {
            for c in 0..cols {
                points.push(vec![x0 + c as Dtype * l_step, y0 + r as Dtype * h_step]);
            }
        }

        for r in 0..(rows - 1) {
            for c in 0..(cols - 1) {
                coupls.push(vec![
                    r * cols + c,
                    r * cols + c + 1,
                    (r + 1) * cols + c + 1,
                    (r + 1) * cols + c,
                ]);
                if rows < 2 || cols < 2 {
                    break;
                }
            }
        }
        (points, coupls)
    }
}
