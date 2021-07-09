use super::group::*;
use super::perm::*;
use super::naive;

/// A schreier vector is a data structure for traversing a group that facilitates two important
/// operations: 
/// 1. Computing the orbit of an element.
/// 2. Computing the set of coset representatives for an element.
#[derive(Debug, Clone)]
struct SchreierVector {
    point: usize,
    // label the i'th generator of the group G as g_i. Then,
    // if v[i] = -1, then i is not in the orbit of `point` under the action of G.
    // if v[i] = k, then k is the index of a generator such that g_k^-1(i) takes you one
    //   step closer to `point`.
    v: Vec<i64>,
}

impl SchreierVector {
    /// Constructs a Schreier vector about the point `point`
    fn new(g: &Group, point: usize) -> SchreierVector {
        let mut r : Vec<i64> = vec![-1; g.num_points()];

        // maintain a stack that holds the frontier
        let mut stack : Vec<usize> = Vec::with_capacity(g.num_points());
        r[point] = Some(g.id().clone());
        stack.push(point);

        // now work on the frontier until it is empty
        while !stack.is_empty() {
            let top = stack.pop().unwrap();
            let prevg = r[top].clone().unwrap();

            // apply all generators to this point and add them to the frontier
            // if they have not yet been explored
            let newpoints : Vec<(&Perm, usize)> = g.get_gens().iter().map(|g| (g, g.apply(top))).collect();
            for i in 0..newpoints.len() {
                let (g, point) = newpoints[i];
                if r[point].is_none() {
                    r[point] = Some(prevg.compose(g));
                    stack.push(point);
                }
            }
        }
        return SchreierVector {
            v: r
        }
    }

    /// Find all coset representatives
    fn coset_repr(&self) -> Vec<Perm> {
        return self.v.iter().filter(|i| i.is_some()).map(|x| x.clone().unwrap()).collect();
    }

    /// Given a point in the group g, find its representative
    fn find_repr(&self, g: &Perm) -> Perm {
        // first compute what it does to the point; this gives us where it is on
        // the orbit
        let curp = g.apply(self.point);
        self.v[curp].clone().unwrap_or(self.g.id().clone()).clone()
    }
}


/// Generate the sub-group of g that stabilizes p using Schreier's lemma
pub fn stab_group(g: &Group, p: usize) -> Group {
    let sv = SchreierVector::new(g, p);

    // apply schreier lemma
    let cosetr = sv.coset_repr();
    let mut gens = Vec::new();
    for p1 in cosetr.iter() {
        for p2 in g.get_gens().iter() {
            let prod = p1.compose(p2);
            let phi = sv.find_repr(&prod).inv();
            gens.push(phi.compose(&prod))
        }
    }
    return Group::new(gens);
}



pub struct StabilizerChain {
    chain: Vec<Group>,
    stab_points: Vec<usize>,
    schreier: Vec<SchreierVector>
}


impl StabilizerChain {
    /// Generate a stabilizer chain for the group g, by default using points [1, 2, ..., n]
    pub fn naive_stab_chain(g: &Group) -> StabilizerChain {
        let mut chain : Vec<Group> = Vec::new();
        let mut curg = g;
        chain.push((*g).clone());
        for i in 0..(g.num_gens()) {
            chain.push(naive::stab_group(curg, i));
            curg = &chain.last().unwrap();
        }
        StabilizerChain { chain, stab_points: (0..(g.num_gens())).collect()}
    }

    /// Compute the size of the group using the stabilizer chain
    pub fn size(&self) -> usize {
        // compute the orbit of each of the stabilizing points
        let mut sz : usize = 1;
        // skip the last stabilizer, which stabilizes all points
        for i in 0..(self.chain.len()-1) {
            sz *= self.chain[i].orbit(self.stab_points[i]).len();
        }
        sz
    }

    /// Checks if g is in the group
    pub fn in_group(&self, p: &Perm) -> bool {
        // Recall from the proof of Schreier's theorem that every element of G
        // can be written in terms of a product of coset representatives from
        // each stabilizer subroup: $ g = h_1 \cdot h_2 \cdot ... \cdot h_n, where h_i \in U_i.
        //
        // We basically rely directly on this proof here: we are deconstructing
        // each group element according to this decomposition.

        // check to see if, in each stabilizer sub-group, there exists a perm
        // that is compatible with p
        for i in 0..self.stab_points.len() {
            let curp = p.apply(stab_points[i]);
            // check if curp is in the orbit

        }
    }
}

#[test]
fn test_naive_stab() {
    let g = Group::cyclic_group(5);
    let chain = StabilizerChain::naive_stab_chain(&g);
    assert_eq!(chain.chain[1].orbit(chain.stab_points[0]).len(), 1);
}

#[test]
fn test_cyclic_naive() {
    let g = Group::cyclic_group(5);
    let chain = StabilizerChain::naive_stab_chain(&g);
    assert_eq!(chain.size(), 5);
}

#[test]
fn test_sym_naive() {
    let g = Group::symmetric(5);
    let chain = StabilizerChain::naive_stab_chain(&g);
    assert_eq!(chain.size(), 120);
}

#[test]
fn test_rubiks_naive() {
    let g = Group::rubiks();
    let chain = StabilizerChain::naive_stab_chain(&g);
    assert_eq!(chain.size(), 0);
}
