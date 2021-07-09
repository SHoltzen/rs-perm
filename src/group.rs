use super::perm::*;
use rand::Rng;

use std::collections::HashSet;



#[derive(Debug, Clone)]
pub struct Group {
    /// vector of generators
    gens: Vec<Perm>, 
     /// number of points
    n: usize,
    id: Perm,
}

impl Group {
    pub fn id(&self) -> &Perm {
        &self.id
    }

    pub fn num_gens(&self) -> usize {
        self.n
    }

    pub fn num_points(&self) -> usize {
        return self.n;
    }

    pub fn get_gens(&self) -> &Vec<Perm> {
        &self.gens
    }

    pub fn new(gens: Vec<Perm>) -> Group {
        assert!(gens.len() > 0);
        let l = gens[0].len();
        for i in 0..gens.len() {
            assert_eq!(l, gens[i].len());
        }
        return Group {
            gens: gens,
            n: l,
            id: Perm::id(l)
        }
    }


    /// applies one step of the product replacement algorithm
    // pub fn rand_elem(&self) -> Perm {
    //     // often 4*n is used as a number of steps; we will use this here
    //     let mut rng = rand::thread_rng();
    //     let i = rng.gen_range(0..(self.num_gens()));
    //     let j = {
    //       let v = rng.gen_range(0..(self.num_gens() - 1));
    //       if v > i { v + 1 } else { v }
    //     };
        
    //     // following notation of https://www.math.ucla.edu/~pak/papers/what9.pdf
    //     let chooseR = 
    // }

    pub fn orbit(&self, p: usize) -> HashSet<usize> {
        let mut frontier: Vec<usize> = Vec::new();
        let mut orbit: HashSet<usize> = HashSet::new();
        // invariant: frontier contains points that are not in the orbit yet
        frontier.push(p);
        while !frontier.is_empty() {
            let top = frontier.pop().unwrap();
            orbit.insert(top);
            // build out set of points not reached yet
            for perm in self.gens.iter() {
                let newpoint = perm.apply(top);
                if !orbit.contains(&newpoint) {
                    frontier.push(newpoint);
                }
            }
        }
        return orbit;
    }

    // Computes generators for the stabilizer subgroup for the point `point`
    // fn stabilizer(&self, point: usize) -> Group {
    //     let sv = SchreierVector::new(self, point);
    //     let coset_repr = sv.coset_repr();
    //     let r : Vec<Perm> = Vec::new();
    // }

    pub fn cyclic_group(n: usize) -> Group {
        let mut gen : Vec<usize> = (1..(n)).collect();
        gen.push(0);
        return Group {
            n,
            id: Perm::id(n),
            gens: vec![Perm::new(gen)]
        }
    }

    pub fn symmetric(n: usize) -> Group {
        let mut gens = Vec::new();
        for i in 0..n {
            for j in (i+1)..n {
                gens.push(Perm::new_cyc(n, vec![vec![i, j]]))
            }
        }
        return Group {
            n,
            id: Perm::id(n),
            gens
        }
    }

    pub fn rubiks() -> Group {
        let gens = vec![ 
            Perm::new_cyc(49, vec![vec![ 1, 3, 8, 6], vec![ 2, 5, 7, 4], vec![ 9,33,25,17], vec![10,34,26,18], vec![11,35,27,19]]),
            Perm::new_cyc(49, vec![vec![ 9,11,16,14], vec![10,13,15,12], vec![ 1,17,41,40], vec![ 4,20,44,37], vec![ 6,22,46,35]]),
            Perm::new_cyc(49, vec![vec![17,19,24,22], vec![18,21,23,20], vec![ 6,25,43,16], vec![ 7,28,42,13], vec![ 8,30,41,11]]),
            Perm::new_cyc(49, vec![vec![25,27,32,30], vec![26,29,31,28], vec![ 3,38,43,19], vec![ 5,36,45,21], vec![ 8,33,48,24]]),
            Perm::new_cyc(49, vec![vec![33,35,40,38], vec![34,37,39,36], vec![ 3, 9,46,32], vec![ 2,12,47,29], vec![ 1,14,48,27]]),
            Perm::new_cyc(49, vec![vec![41,43,48,46], vec![42,45,47,44], vec![14,22,30,38], vec![15,23,31,39], vec![16,24,32,40]])
            ];
        return Group {n: 49, id: Perm::id(49), gens}
    }
}

#[test]
fn test_orb_cyclic() {
    let g = Group::cyclic_group(5);
    assert_eq!(g.orbit(0).len(), 5);
}

#[test]
fn test_orb_symmetric() {
    let g = Group::symmetric(5);
    assert_eq!(g.orbit(0).len(), 5);
    assert_eq!(g.orbit(1).len(), 5);
    assert_eq!(g.orbit(2).len(), 5);
    assert_eq!(g.orbit(3).len(), 5);
    assert_eq!(g.orbit(4).len(), 5);
}
