// membrane Phase separation module

pub struct PhaseSeparation {
    pub membrane: Membrane,
}

impl PhaseSeparation {
    pub fn new(membrane: Membrane) -> Self {
        Self { membrane }
    }
}       