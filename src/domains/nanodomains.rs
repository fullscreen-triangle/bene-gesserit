//  membrane Nanodomains module

pub struct Nanodomains {
    pub membrane: Membrane,
}

impl Nanodomains {
    pub fn new(membrane: Membrane) -> Self {
        Self { membrane }
    }
}   