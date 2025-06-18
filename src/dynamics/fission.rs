//  membrane Fission dynamics module

pub struct Fission {
    pub membrane: Membrane,
}

impl Fission {
    pub fn new(membrane: Membrane) -> Self {
        Self { membrane }
    }
}   