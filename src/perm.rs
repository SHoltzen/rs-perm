#[derive(Debug, Eq, Clone)]
pub struct Perm {
    p: Vec<usize>
}

impl Perm {
    pub fn apply(&self, point: usize) -> usize {
        return self.p[point]
    }

    /// Create a new permutation using matrix notation
    pub fn new(v: Vec<usize>) -> Perm {
        return Perm { p: v }
    }

    /// Create a new permutation using cycle representation
    pub fn new_cyc(sz: usize, cyc: Vec<Vec<usize>>) -> Perm {
        // convert it into a cycle
        let mut p : Vec<usize> = (0..sz).collect();
        for c in cyc {
            for i in 0..c.len() {
                p[c[i]] = c[(i+1) % c.len()]
            }
        }
        Perm {p}
    }

    pub fn id(sz: usize) -> Perm {
        return Perm::new((0..(sz)).collect())
    }

    pub fn len(&self) -> usize {
        return self.p.len()
    }

    pub fn inv(&self) -> Perm {
        let mut r = vec![0; self.p.len()];
        for i in 0..(self.p.len()) {
            r[self.p[i]] = i
        }
        return Perm::new(r)
    }

    /// compute self * b (apply the permutation b, then self)
    pub fn compose(&self, b: &Perm) -> Perm {
        assert_eq!(self.p.len(), b.p.len());
        let mut r = vec![0; self.p.len()];
        for i in 0..(self.p.len()) {
            r[i] = self.apply(b.apply(i))
        }
        return Perm::new(r);
    }
}

impl PartialEq for Perm {
    fn eq(&self, other: &Self) -> bool {
        self.p == other.p
    }
}



#[test]
fn compose_test() {
    let p1 = Perm::new(vec![0, 2, 1]);
    let p2 = p1.compose(&p1);
    assert_eq!(p2, Perm::id(3))
}

#[test]
fn inv_test() {
    let p1 = Perm::new(vec![0, 2, 1, 4, 7, 8, 3, 5, 6]);
    assert_eq!(p1.compose(&p1.inv()), Perm::id(9))
}

