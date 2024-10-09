use crate::Dtype;
use std::fmt;

pub struct Node1D {
    pub id: usize,
    pub coords: [Dtype; 1],
    pub displs: [Dtype; 1],
    pub forces: [Dtype; 1],
}

impl Node1D {
    pub fn new(id: usize, coords: [Dtype; 1]) -> Node1D {
        Node1D {
            id,
            coords,
            displs: [0.0],
            forces: [0.0],
        }
    }
}

impl fmt::Display for Node1D {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "\n\t\tNode{}:  1D\n\t\t\tCoord: [{:-7.4}]\n\t\t\tDisps: [{:-7.4}]\n\t\t\tForce: [{:-7.4}]",
            self.id,
            self.coords[0],
            self.displs[0],
            self.forces[0],
        )
    }
}

#[derive(Copy, Clone)]
pub struct Node2D {
    pub id: usize,
    pub coords: [Dtype; 2],
    pub displs: [Dtype; 2],
    pub forces: [Dtype; 2],
}

impl Node2D {
    pub fn new(id: usize, coords: [Dtype; 2]) -> Node2D {
        Node2D {
            id,
            coords,
            displs: [0.0, 0.0],
            forces: [0.0, 0.0],
        }
    }
}

impl fmt::Display for Node2D {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "\n\t\tNode{}:  2D\n\t\t\tCoord: [{:-7.4}, {:-7.4}]\n\t\t\tDisps: [{:-7.4}, {:-7.4}]\n\t\t\tForce: [{:-7.4}, {:-7.4}]",
            self.id,
            self.coords[0],
            self.coords[1],
            self.displs[0],
            self.displs[1],
            self.forces[0],
            self.forces[1],
        )
    }
}

pub struct Node3D {
    pub id: usize,
    pub coords: [Dtype; 3],
    pub displs: [Dtype; 3],
    pub forces: [Dtype; 3],
}

impl Node3D {
    pub fn new(id: usize, coords: [Dtype; 3]) -> Node3D {
        Node3D {
            id,
            coords,
            displs: [0.0; 3],
            forces: [0.0; 3],
        }
    }
}

impl fmt::Display for Node3D {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "\nNode_3D Info:\n\tId:    {}\n\tCoord: [{}, {}, {}]",
            self.id, self.coords[0], self.coords[1], self.coords[2]
        )
    }
}
