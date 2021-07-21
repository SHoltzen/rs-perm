pub mod perm;
pub mod group;
pub mod naive;

extern crate rand;

use perm::*;
use group::*;


fn main() {
    let mut g = Group::rubiks();
    g.random_schreier_sims();
    let sz = 3;
    // let sz : usize = (2 as usize).pow(27) * (3 as usize).pow(14) * (5 as usize).pow(3) * (7 as usize).pow(2) * 11;
    assert_eq!(g.sz(), sz);
}

