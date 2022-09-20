use std::fmt;

#[derive(Debug)]
pub struct Node1D {
    pub id: usize,
    pub coord: [f64; 1],
}

impl Node1D {
    pub fn new(id: usize, coord: [f64; 1]) -> Node1D {
        Node1D { id, coord }
    }
}

impl fmt::Display for Node1D {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "\nNode_1D Info:\n\tId:    {}\n\tCoord: {:?}",
            self.id, self.coord
        )
    }
}

#[derive(Debug)]
pub struct Node2D {
    pub id: usize,
    pub coord: [f64; 2],
}

impl Node2D {
    pub fn new(id: usize, coord: [f64; 2]) -> Node2D {
        Node2D { id, coord }
    }
}

impl fmt::Display for Node2D {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "\nNode_2D Info:\n\tId:    {}\n\tCoord: {:?}",
            self.id, self.coord
        )
    }
}

#[derive(Debug)]
pub struct Node3D {
    pub id: usize,
    pub coord: [f64; 3],
}

impl Node3D {
    pub fn new(id: usize, coord: [f64; 3]) -> Node3D {
        Node3D { id, coord }
    }
}

impl fmt::Display for Node3D {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "\nNode_3D Info:\n\tId:    {}\n\tCoord: {:?}",
            self.id, self.coord
        )
    }
}
