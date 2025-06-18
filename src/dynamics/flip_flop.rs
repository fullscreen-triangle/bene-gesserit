//  membrane Flip-Flop dynamics module

pub struct FlipFlop {
    pub membrane: Membrane,
}

impl FlipFlop {
    pub fn new(membrane: Membrane) -> Self {
        Self { membrane }
    }
}           