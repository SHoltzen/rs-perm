use super::perm::*;
use rand::Rng;

use std::collections::HashSet;

pub struct GroupRNG {
    gens: Vec<Perm>,
    n: usize,
    id: Perm,
    rng: rand::rngs::ThreadRng
}

impl GroupRNG {
    pub fn new(g: &Group) -> GroupRNG {
        GroupRNG {
            gens: g.gens.clone(),
            n: g.n,
            id: g.id.clone(),
            rng: rand::thread_rng()
        }
    } 

    /// applies one step of the product replacement algorithm
    fn step(&mut self) {
        // special case is needed for cyclic group: product replacement does not apply
        // here, since there is only 1 generator. In this case, we apply a 
        // random power to the basis point
        
        

        // following notation of https://www.math.ucla.edu/~pak/papers/what9.pdf
        let i = self.rng.gen_range(0..(self.gens.len()));
        let j = {
          let v = self.rng.gen_range(0..(self.gens.len() - 1));
          if v >= i { v + 1 } else { v }
        };
        let pi = &self.gens[i];
        let pj = &self.gens[j];
        
        let chooseR : bool = self.rng.gen_bool(0.5);
        let chooseplus : bool = self.rng.gen_bool(0.5);
        if chooseR {
            if chooseplus {
                self.gens[i] = pi.compose(&pj);
            } else {
                self.gens[i] = pi.compose(&pj.inv())
            }
        } else {
            if chooseplus {
                self.gens[i] = pj.compose(&pi);
            } else {
                self.gens[i] = pj.inv().compose(&pi);
            }
        }
    }

    fn draw(&mut self) -> Perm {
        // use 4*n as number of steps, this seems pretty standard and is what GAP uses
        for _ in 0..(4*self.n) {
            self.step();
        }

        // return a random generator
        let i = self.rng.gen_range(0..(self.gens.len()));
        self.gens[i].clone()
    }
}

#[derive(Debug, Clone)]
struct StabView<'a> {
    g: Vec<&'a Perm>,
    /// the set of stabilized points
    base: Vec<usize>,
    /// the next basis point (to be used for schreier)
    point: usize,
    id: &'a Perm,
    schreier: Vec<Option<&'a Perm>>
}

impl<'a> StabView<'a> {
    fn new(g: &'a Group, base: Vec<usize>, point: usize) -> StabView<'a> {
        // pick out all generators that stabilize the current base point
        let gens : Vec<&'a Perm> = g.gens.iter().filter(|g| 
            base.iter().all(|p| 
                g.apply(*p) == *p)
            ).collect();

        // now, build the schreier vector by inspecting how `gens` transforms `point`; this 
        // is essentially an orbit computation
        let mut s = vec![None; g.n];
        let mut frontier = Vec::new();
        s[point] = Some(g.id());
        frontier.push(point);
        while !frontier.is_empty() {
            let top = frontier.pop().unwrap();
            for g in gens.iter() {
                let newp = g.apply(top);
                if !s[newp].is_none() { continue }
                s[newp] = Some(g);
                frontier.push(newp);
            }
        }

        StabView { g: gens, base, schreier: s, point, id: g.id()}
    }

    pub fn add_gen(&mut self, g: &'a Perm) -> () {
        // first filter it: see if it fixes the base
        if !self.base.iter().all(|p| g.apply(*p) == *p) { return }

        // now we know g fixes the base; add it to the set of generators and
        // re-build schreier
        self.g.push(g);
        for i in 0..(self.schreier.len()) {
            self.schreier[i] = None;
        }
        let mut frontier = Vec::new();
        let point = self.point;
        self.schreier[point] = Some(self.id);
        frontier.push(point);
        while !frontier.is_empty() {
            let top = frontier.pop().unwrap();
            for g in self.g.iter() {
                let newp = g.apply(top);
                if !self.schreier[newp].is_none() { continue }
                self.schreier[newp] = Some(g);
                frontier.push(newp);
            }
        }
    }

    /// Find the representative of g
    pub fn repr(&self, g: &Perm) -> Perm {
        // apply g to the point and return that entry
        let mut curp = g.apply(self.point);
        if self.schreier[curp].is_none() {
            return self.id.clone();
        }

        let mut r = self.id.clone();
        // work backwards until we are back at the original point
        while curp != self.point {
            let curg = self.schreier[curp].unwrap();
            curp = curg.inv().apply(curp);
            // println!("curp: {:?}, curg: {:?}", curp, curg);
            r = curg.compose(&r);
        }
        return r;
    }


    pub fn base_orbit_size(&self) -> usize {
        // return # non-none entries in the schreier vector
        return self.schreier.iter().map(|x| if x.is_none() { 0 } else { 1 }).sum()
    }
}

#[derive(Debug, Clone)]
pub struct Group {
    /// vector of generators
    gens: Vec<Perm>, 
    /// number of points
    n: usize,
    /// identity element
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

    // strip the point according to the stabilizer chain
    // returns the residue
    fn strip(&self, stabchain: &Vec<StabView>, g: &Perm) -> Perm {
        let mut curg = g.clone();
        // printing out the trace of this is quite educational!
        for s in stabchain.iter() {
            let curp = curg.apply(s.point);
            if s.schreier[curp].is_none() {
                return curg;
            }
            curg = s.repr(&curg).inv().compose(&curg);
        }
        return curg;
    }

    fn gen_stab_chain(&self, base: &Vec<usize>) -> Vec<StabView> {
        let mut stab_chain = Vec::new();
        let mut sub_base = Vec::new();
        for i in 0..(base.len()) {
            stab_chain.push(StabView::new(self, sub_base.clone(), i));
            sub_base.push(i);
        }
        return stab_chain;
    }

    /// Expand the generating set of this group to include strong generators
    pub fn random_schreier_sims(&mut self) -> () {
        let mut count = 0;

        let mut grng = GroupRNG::new(self);
        let base : Vec<usize> = (0..self.n).collect();
        // generate new candidate schreier vectors until 4 in a row do not strip
        // to new generators
        while count < self.n {
            // make stabilizer chain
            let stab_chain = self.gen_stab_chain(&base);

            // get a random group element and attempt to sift it if it is
            // already in the stabilizer chain, it is not a new strong generator
            // and we continue; otherwise, we add it to the generating set of
            // the group
            let rnd = grng.draw();
            let stripped = self.strip(&stab_chain, &rnd);
            if stripped != *self.id() {
                self.gens.push(stripped);
                println!("new generator found, total generators: {}", self.gens.len());
                count = 0;
            } else {
                count += 1;
            }
        }
    }

    /// Compute the size of the group
    /// Assumes the group is generated by an SGS
    pub fn sz(&self) -> usize {
        let base : Vec<usize> = (0..self.n).collect();
        let chain = self.gen_stab_chain(&base);
        return chain.iter().map(|x| x.base_orbit_size()).product()
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
            id: Perm::id(l),
        }
    }


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

#[test]
fn test_repr() {
    // start with a group that is already an SGS; this is an SGS for S4
    let g = Group::new(vec![
        Perm::new_cyc(4, vec![vec![0, 1, 2, 3]]),
        Perm::new_cyc(4, vec![vec![1, 2, 3]]),
        Perm::new_cyc(4, vec![vec![2, 3]]),
    ]);

    // generate the stabilizer chain
    let s1 = StabView::new(&g, vec![], 0);

    // two perms are in the same coset if they both map the basis point
    // to the same point, so let's verify that fact; both these cycles map
    // the basis point 0 to 1 
    let rep1 = s1.repr(&Perm::new_cyc(4, vec![vec![0, 1, 2, 3]]));
    let rep2 = s1.repr(&Perm::new_cyc(4, vec![vec![0, 1]]));
    let rep3 = s1.repr(&Perm::new_cyc(4, vec![vec![0, 2]]));
    assert_eq!(rep1, rep2);
    assert_eq!(rep1.apply(0), 1);
    assert_ne!(rep1, rep3);
}

#[test]
fn test_strip() {
    // start with a group that is already an SGS; this is an SGS for S4
    let g = Group::new(vec![
        Perm::new_cyc(4, vec![vec![0, 1, 2, 3]]),
        Perm::new_cyc(4, vec![vec![1, 2, 3]]),
        Perm::new_cyc(4, vec![vec![2, 3]]),
    ]);

    // generate the stabilizer chain
    let chain = g.gen_stab_chain(&(0..3).collect());
    // let s1 = StabView::new(&g, vec![], 0);
    // let s2 = StabView::new(&g, vec![0], 1);
    // let s3 = StabView::new(&g, vec![0, 1], 2);
    // let s4 = StabView::new(&g, vec![0, 1, 2], 3);
    // let chain = vec![s1, s2, s3, s4];

    // now see if we can properly strip elements of S4
    assert_eq!(g.strip(&chain, &Perm::new_cyc(4, vec![vec![0, 3]])), *g.id());
    assert_eq!(g.strip(&chain, &Perm::new_cyc(4, vec![vec![1, 2, 3]])), *g.id());
    assert_eq!(g.strip(&chain, &Perm::new_cyc(4, vec![vec![0, 2, 3]])), *g.id());
    assert_eq!(g.strip(&chain, &Perm::new_cyc(4, vec![vec![3, 1, 2]])), *g.id());
}

#[test]
fn test_size() {
    // start with a group that is already an SGS; this is an SGS for S4
    let g = Group::new(vec![
        Perm::new_cyc(4, vec![vec![0, 1, 2, 3]]),
        Perm::new_cyc(4, vec![vec![1, 2, 3]]),
        Perm::new_cyc(4, vec![vec![2, 3]]),
    ]);
    assert_eq!(g.sz(), 4*3*2);
}

#[test]
fn test_schreier() {
    let mut g = Group::symmetric(12);
    // expand basis to include strong generators
    g.random_schreier_sims();
    assert_eq!(g.sz(), 12 * 11 * 10 * 9 * 8*7*6*5*4*3*2);
}

#[test]
fn test_schreier_rubiks() {
    let mut g = Group::rubiks();
    g.random_schreier_sims();
    let sz = 3;
    // let sz : usize = (2 as usize).pow(27) * (3 as usize).pow(14) * (5 as usize).pow(3) * (7 as usize).pow(2) * 11;
    assert_eq!(g.sz(), sz);
}