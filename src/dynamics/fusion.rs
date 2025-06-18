//  membrane Fusion dynamics module

pub struct Fusion {
    pub membrane: Membrane,
}

impl Fusion {
    pub fn new(membrane: Membrane) -> Self {
        Self { membrane }
    }
}       