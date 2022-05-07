#[macro_use(c)]
extern crate cute;

mod elem;
mod node;

pub use elem::*;
pub use node::*;

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn gen_nodes() {
        let node1 = Node1D::new(1, [1.0]);
        let node2 = Node2D::new(2, [1.0, 2.0]);
        let node3 = Node3D::new(3, [1.0, 2.0, 3.0]);

        assert_eq!(1usize, node1.id);
        assert_ne!(3usize, node2.id);
        assert_eq!(3usize, node3.id);

        assert_eq!([1.0f64, 2.0f64], node2.coord);
        assert_eq!([1.0f64, 2.0f64, 3.0f64], node3.coord);
        assert_eq!([1.0f64], node1.coord);
    }

    #[test]
    fn gen_elements() {
        let node1 = Node2D::new(1, [0.0, 0.0]);
        let node2 = Node2D::new(2, [0.0, 1.0]);
        let node3 = Node2D::new(3, [1.0, 0.0]);
        let node4 = Node2D::new(4, [1.0, 1.0]);

        let tri1 = Triangle::new(1, [&node1, &node2, &node3]);
        let tri2 = Triangle::new(2, [&node4, &node2, &node3]);

        let rec1 = Rectangle::new(3, [&node1, &node2, &node3, &node4]);

        assert_eq!(1usize, tri1.id);
        assert_ne!(2usize, tri1.id);
        assert_eq!(3usize, rec1.id);

        assert_ne!([10.0f64, 11.0f64], tri2.nodes[0].coord);
        assert_eq!([1.0f64, 1.0f64], tri2.nodes[0].coord);

        assert_eq!(vec![0.0, 0.0, 1.0], tri1.xs());
        assert_eq!(vec![0.0, 0.0, 1.0, 1.0], rec1.get_xs());
        assert_ne!(vec![0.0, 0.0, 1.0], tri1.ys());
        assert_ne!(vec![0.0, 0.0, 1.0, 1.0], rec1.get_ys());
    }
}
