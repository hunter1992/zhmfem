use crate::Dtype;
use std::cell::RefCell;
use std::fmt;

pub struct Node1D {
    pub id: usize,
    pub coord: [Dtype; 1],
    pub disps: [RefCell<Dtype>; 1],
    pub forces: [RefCell<Dtype>; 1],
}

impl Node1D {
    pub fn new(id: usize, coord: [Dtype; 1]) -> Node1D {
        Node1D {
            id,
            coord,
            disps: [RefCell::new(0.0)],
            forces: [RefCell::new(0.0)],
        }
    }
}

impl fmt::Display for Node1D {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "\n\t\tNode{}:  1D\n\t\t\tCoord: [{:-7.4}]\n\t\t\tDisps: [{:-7.4}]\n\t\t\tForce: [{:-7.4}]",
            self.id,
            self.coord[0],
            self.disps[0].borrow(),
            self.forces[0].borrow(),
        )
    }
}

pub struct Node2D {
    pub id: usize,
    pub coord: [Dtype; 2],
    pub disps: [RefCell<Dtype>; 2],
    pub forces: [RefCell<Dtype>; 2],
}

impl Node2D {
    pub fn new(id: usize, coord: [Dtype; 2]) -> Node2D {
        Node2D {
            id,
            coord,
            disps: [RefCell::new(0.0), RefCell::new(0.0)],
            forces: [RefCell::new(0.0), RefCell::new(0.0)],
        }
    }
}

impl fmt::Display for Node2D {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "\n\t\tNode{}:  2D\n\t\t\tCoord: [{:-7.4}, {:-7.4}]\n\t\t\tDisps: [{:-7.4}, {:-7.4}]\n\t\t\tForce: [{:-7.4}, {:-7.4}]",
            self.id,
            self.coord[0],
            self.coord[1],
            self.disps[0].borrow(),
            self.disps[1].borrow(),
            self.forces[0].borrow(),
            self.forces[1].borrow()
        )
    }
}

pub struct Node3D {
    pub id: usize,
    pub coord: [Dtype; 3],
    pub disps: [RefCell<Dtype>; 3],
    pub forces: [RefCell<Dtype>; 3],
}

impl Node3D {
    pub fn new(id: usize, coord: [Dtype; 3]) -> Node3D {
        Node2D {
            id,
            coord,
            disps: [RefCell::new(0.0), RefCell::new(0.0), RefCell::new(0.0)],
            forces: [RefCell::new(0.0), RefCell::new(0.0), RefCell::new(0.0)],
        }
    }
}

impl fmt::Display for Node3D {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "\nNode_3D Info:\n\tId:    {}\n\tCoord: [{}, {}, {}]",
            self.id, self.coord[0], self.coord[1], self.coord[2]
        )
    }
}
