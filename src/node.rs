use crate::{Data, Dtype};
use std::cell::RefCell;
use std::fmt;

pub struct Node1D {
    pub id: usize,
    pub coords: [Dtype; 1],
    pub displs: RefCell<[Dtype; 1]>,
    pub forces: RefCell<[Dtype; 1]>,
}

impl Node1D {
    pub fn new(id: usize, coords: [Dtype; 1]) -> Node1D {
        Node1D {
            id,
            coords,
            displs: RefCell::new([0.0]),
            forces: RefCell::new([0.0]),
        }
    }
}

impl fmt::Display for Node1D {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f,"\n\t\tNode_1D:\n\t\t\tId: {}\n\t\t\tCoords: [{:-10.4}]\n\t\t\tDispls: [{:-10.4}]\n\t\t\tForces: [{:-10.4}]\n",self.id, self.coords[0], self.displs.borrow()[0],self.forces.borrow()[0])
    }
}

pub struct Node2D {
    pub id: usize,
    pub coords: [Dtype; 2],
    pub displs: RefCell<[Dtype; 2]>,
    pub forces: RefCell<[Dtype; 2]>,
    pub stress: Vec<[Dtype; 3]>,
}

impl Node2D {
    pub fn new(id: usize, coords: [Dtype; 2]) -> Node2D {
        Node2D {
            id,
            coords,
            displs: RefCell::new([0.0, 0.0]),
            forces: RefCell::new([0.0, 0.0]),
            stress: Vec::with_capacity(6),
        }
    }
}

impl fmt::Display for Node2D {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f,"\n\t\tNode_2D:\n\t\t\tId: {}\n\t\t\tCoords: [{:-10.4}, {:-10.4}]\n\t\t\tDispls: [{:-10.4}, {:-10.4}]\n\t\t\tForces: [{:-10.4}, {:-10.4}]\n", self.id, self.coords[0], self.coords[1], self.displs.borrow()[0], self.displs.borrow()[1], self.forces.borrow()[0], self.forces.borrow()[1])
    }
}

pub struct Node3D {
    pub id: usize,
    pub coords: [Dtype; 3],
    pub displs: RefCell<[Dtype; 3]>,
    pub forces: RefCell<[Dtype; 3]>,
}

impl Node3D {
    pub fn new(id: usize, coords: [Dtype; 3]) -> Node3D {
        Node3D {
            id,
            coords,
            displs: RefCell::new([0.0, 0.0, 0.0]),
            forces: RefCell::new([0.0, 0.0, 0.0]),
        }
    }
}

impl fmt::Display for Node3D {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f,"\n\t\tNode_2D:\n\t\t\tId: {}\n\t\t\tCoords: [{:-10.4}, {:-10.4}, {:-10.4}]\n\t\t\tDispls: [{:-10.4}, {:-10.4}, {:-10.4}]\n\t\t\tForces: [{:-10.4}, {:-10.4}, {:-10.4}]\n", self.id, self.coords[0], self.coords[1], self.coords[2], self.displs.borrow()[0], self.displs.borrow()[1], self.displs.borrow()[2], self.forces.borrow()[0], self.forces.borrow()[1], self.forces.borrow()[2])
    }
}
