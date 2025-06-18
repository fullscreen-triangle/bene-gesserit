//! Membrane diffusion dynamics
//!
//! This module simulates the diffusion of lipids and proteins within
//! the membrane, including:
//! - Lateral diffusion of membrane components
//! - Brownian motion simulation
//! - Diffusion coefficients calculation
//! - ATP-dependent active transport

use crate::{
    constants::*,
    error::{MembraneError, Result},
    types::*,
};
use std::collections::HashMap;

/// Diffusion types in membranes
#[derive(Debug, Clone, PartialEq)]
pub enum DiffusionType {
    /// Free lateral diffusion
    Free,
    /// Hindered diffusion (obstacles)
    Hindered,
    /// Active transport (ATP-dependent)
    Active,
    /// Confined diffusion (lipid rafts)
    Confined,
}

/// Diffusing particle properties
#[derive(Debug, Clone)]
pub struct DiffusingParticle {
    /// Particle identifier
    pub id: String,
    /// Current position (nm)
    pub position: (f64, f64),
    /// Velocity (nm/s)
    pub velocity: (f64, f64),
    /// Diffusion coefficient (nm²/s)
    pub diffusion_coefficient: f64,
    /// Particle radius (nm)
    pub radius: f64,
    /// Particle type
    pub particle_type: ParticleType,
    /// Current diffusion mode
    pub diffusion_mode: DiffusionType,
}

/// Types of diffusing particles
#[derive(Debug, Clone, PartialEq)]
pub enum ParticleType {
    /// Lipid molecule
    Lipid(LipidType),
    /// Membrane protein
    Protein(ProteinType),
    /// Ion or small molecule
    Ion(IonType),
    /// Other membrane component
    Other(String),
}

/// Diffusion simulation manager
#[derive(Debug, Clone)]
pub struct DiffusionSimulator {
    /// Particles being tracked
    pub particles: HashMap<String, DiffusingParticle>,
    /// Simulation parameters
    pub parameters: DiffusionParameters,
    /// Membrane boundaries
    pub boundaries: MembraneGeometry,
    /// Diffusion statistics
    pub statistics: DiffusionStatistics,
}

/// Diffusion simulation parameters
#[derive(Debug, Clone)]
pub struct DiffusionParameters {
    /// Simulation timestep (s)
    pub timestep: f64,
    /// Temperature (K)
    pub temperature: f64,
    /// Membrane viscosity (Pa·s)
    pub viscosity: f64,
    /// Random seed for reproducibility
    pub random_seed: u64,
    /// Enable ATP-dependent transport
    pub atp_transport_enabled: bool,
}

/// Membrane geometry for diffusion boundaries
#[derive(Debug, Clone)]
pub struct MembraneGeometry {
    /// Membrane patch dimensions (nm)
    pub width: f64,
    pub height: f64,
    /// Periodic boundary conditions
    pub periodic_boundaries: bool,
    /// Obstacles and constraints
    pub obstacles: Vec<Obstacle>,
}

/// Obstacle in membrane diffusion
#[derive(Debug, Clone)]
pub struct Obstacle {
    /// Obstacle position (nm)
    pub position: (f64, f64),
    /// Obstacle radius (nm)
    pub radius: f64,
    /// Obstacle type
    pub obstacle_type: ObstacleType,
}

/// Types of diffusion obstacles
#[derive(Debug, Clone, PartialEq)]
pub enum ObstacleType {
    /// Protein complex
    ProteinComplex,
    /// Lipid raft
    LipidRaft,
    /// Cytoskeletal attachment
    Cytoskeleton,
    /// Temporary barrier
    Temporary,
}

/// Diffusion statistics
#[derive(Debug, Clone)]
pub struct DiffusionStatistics {
    /// Mean square displacement over time
    pub msd: HashMap<String, Vec<f64>>,
    /// Calculated diffusion coefficients
    pub calculated_diffusion_coeffs: HashMap<String, f64>,
    /// Particle trajectories
    pub trajectories: HashMap<String, Vec<(f64, f64)>>,
    /// Collision counts
    pub collision_counts: HashMap<String, u64>,
}

impl DiffusionSimulator {
    /// Create a new diffusion simulator
    pub fn new(geometry: MembraneGeometry) -> Self {
        Self {
            particles: HashMap::new(),
            parameters: DiffusionParameters::default(),
            boundaries: geometry,
            statistics: DiffusionStatistics::new(),
        }
    }
    
    /// Add a particle to the simulation
    pub fn add_particle(&mut self, particle: DiffusingParticle) -> Result<()> {
        // Validate particle position within boundaries
        if !self.is_position_valid(particle.position) {
            return Err(MembraneError::InvalidParameter(
                format!("Particle position {:?} outside boundaries", particle.position)
            ));
        }
        
        // Initialize statistics for this particle
        self.statistics.msd.insert(particle.id.clone(), Vec::new());
        self.statistics.trajectories.insert(particle.id.clone(), vec![particle.position]);
        self.statistics.collision_counts.insert(particle.id.clone(), 0);
        
        self.particles.insert(particle.id.clone(), particle);
        Ok(())
    }
    
    /// Check if position is within membrane boundaries
    fn is_position_valid(&self, position: (f64, f64)) -> bool {
        position.0 >= 0.0 && position.0 <= self.boundaries.width &&
        position.1 >= 0.0 && position.1 <= self.boundaries.height
    }
    
    /// Simulate one diffusion step
    pub fn step(&mut self, dt: f64) -> Result<()> {
        let mut new_positions = HashMap::new();
        
        // Calculate new positions for all particles
        for (particle_id, particle) in &self.particles {
            let new_position = self.calculate_new_position(particle, dt)?;
            new_positions.insert(particle_id.clone(), new_position);
        }
        
        // Apply new positions and handle collisions
        for (particle_id, new_position) in new_positions {
            if let Some(particle) = self.particles.get_mut(&particle_id) {
                let old_position = particle.position;
                
                // Check for collisions and boundaries
                let final_position = self.handle_boundaries_and_collisions(
                    &particle_id, old_position, new_position
                )?;
                
                particle.position = final_position;
                
                // Update velocity
                particle.velocity = (
                    (final_position.0 - old_position.0) / dt,
                    (final_position.1 - old_position.1) / dt,
                );
                
                // Update statistics
                self.update_particle_statistics(&particle_id, final_position);
            }
        }
        
        Ok(())
    }
    
    /// Calculate new position for a particle using Brownian motion
    fn calculate_new_position(&self, particle: &DiffusingParticle, dt: f64) -> Result<(f64, f64)> {
        match particle.diffusion_mode {
            DiffusionType::Free => self.calculate_free_diffusion(particle, dt),
            DiffusionType::Hindered => self.calculate_hindered_diffusion(particle, dt),
            DiffusionType::Active => self.calculate_active_transport(particle, dt),
            DiffusionType::Confined => self.calculate_confined_diffusion(particle, dt),
        }
    }
    
    /// Calculate free Brownian diffusion
    fn calculate_free_diffusion(&self, particle: &DiffusingParticle, dt: f64) -> Result<(f64, f64)> {
        // Random displacement based on diffusion coefficient
        let sigma = (2.0 * particle.diffusion_coefficient * dt).sqrt();
        
        // Generate random displacements (simplified)
        let dx = sigma * 0.1; // Placeholder for Gaussian random
        let dy = sigma * 0.1; // Placeholder for Gaussian random
        
        let new_x = particle.position.0 + dx;
        let new_y = particle.position.1 + dy;
        
        Ok((new_x, new_y))
    }
    
    /// Calculate hindered diffusion (reduced mobility)
    fn calculate_hindered_diffusion(&self, particle: &DiffusingParticle, dt: f64) -> Result<(f64, f64)> {
        // Reduce diffusion coefficient based on local environment
        let hindrance_factor = self.calculate_hindrance_factor(particle.position);
        let effective_diffusion = particle.diffusion_coefficient * hindrance_factor;
        
        let sigma = (2.0 * effective_diffusion * dt).sqrt();
        let dx = sigma * 0.1; // Placeholder
        let dy = sigma * 0.1; // Placeholder
        
        let new_x = particle.position.0 + dx;
        let new_y = particle.position.1 + dy;
        
        Ok((new_x, new_y))
    }
    
    /// Calculate active transport (biased random walk)
    fn calculate_active_transport(&self, particle: &DiffusingParticle, dt: f64) -> Result<(f64, f64)> {
        // Active transport has both random and directed components
        let active_velocity = self.calculate_active_velocity(particle)?;
        
        // Directed displacement
        let directed_dx = active_velocity.0 * dt;
        let directed_dy = active_velocity.1 * dt;
        
        // Random displacement (reduced compared to free diffusion)
        let sigma = (0.5 * particle.diffusion_coefficient * dt).sqrt();
        let random_dx = sigma * 0.1; // Placeholder
        let random_dy = sigma * 0.1; // Placeholder
        
        let new_x = particle.position.0 + directed_dx + random_dx;
        let new_y = particle.position.1 + directed_dy + random_dy;
        
        Ok((new_x, new_y))
    }
    
    /// Calculate confined diffusion within boundaries
    fn calculate_confined_diffusion(&self, particle: &DiffusingParticle, dt: f64) -> Result<(f64, f64)> {
        // Find confining region (e.g., lipid raft)
        let confinement = self.find_confinement_region(particle.position);
        
        // Normal diffusion within confined region
        let sigma = (2.0 * particle.diffusion_coefficient * dt).sqrt();
        let dx = sigma * 0.1; // Placeholder
        let dy = sigma * 0.1; // Placeholder
        
        let mut new_x = particle.position.0 + dx;
        let mut new_y = particle.position.1 + dy;
        
        // Apply confinement boundaries
        if let Some(conf) = confinement {
            // Reflect off confinement boundaries
            new_x = new_x.max(conf.0).min(conf.2);
            new_y = new_y.max(conf.1).min(conf.3);
        }
        
        Ok((new_x, new_y))
    }
    
    /// Calculate hindrance factor based on local obstacles
    fn calculate_hindrance_factor(&self, position: (f64, f64)) -> f64 {
        let mut hindrance = 1.0;
        
        for obstacle in &self.boundaries.obstacles {
            let distance = self.calculate_distance(position, obstacle.position);
            let influence_range = obstacle.radius * 2.0;
            
            if distance < influence_range {
                // Reduce mobility near obstacles
                let reduction = 1.0 - (influence_range - distance) / influence_range * 0.8;
                hindrance *= reduction;
            }
        }
        
        hindrance.max(0.1) // Minimum 10% mobility
    }
    
    /// Calculate active transport velocity
    fn calculate_active_velocity(&self, particle: &DiffusingParticle) -> Result<(f64, f64)> {
        // ATP-dependent transport velocity
        let base_velocity = 100.0; // nm/s
        let angle = 0.0; // Could be dynamic based on gradients
        
        Ok((
            base_velocity * angle.cos(),
            base_velocity * angle.sin(),
        ))
    }
    
    /// Find confinement region
    fn find_confinement_region(&self, _position: (f64, f64)) -> Option<(f64, f64, f64, f64)> {
        None // Placeholder
    }
    
    /// Handle boundary conditions and collisions
    fn handle_boundaries_and_collisions(&mut self, particle_id: &str, 
                                      old_pos: (f64, f64), 
                                      new_pos: (f64, f64)) -> Result<(f64, f64)> {
        let mut final_pos = new_pos;
        
        // Handle membrane boundaries
        if self.boundaries.periodic_boundaries {
            // Periodic boundary conditions (wrap around)
            final_pos.0 = ((final_pos.0 % self.boundaries.width) + self.boundaries.width) % self.boundaries.width;
            final_pos.1 = ((final_pos.1 % self.boundaries.height) + self.boundaries.height) % self.boundaries.height;
        } else {
            // Reflective boundary conditions
            if final_pos.0 < 0.0 {
                final_pos.0 = -final_pos.0;
            } else if final_pos.0 > self.boundaries.width {
                final_pos.0 = 2.0 * self.boundaries.width - final_pos.0;
            }
            
            if final_pos.1 < 0.0 {
                final_pos.1 = -final_pos.1;
            } else if final_pos.1 > self.boundaries.height {
                final_pos.1 = 2.0 * self.boundaries.height - final_pos.1;
            }
        }
        
        // Check for collisions with obstacles
        for obstacle in &self.boundaries.obstacles {
            let distance = self.calculate_distance(final_pos, obstacle.position);
            if distance < obstacle.radius {
                // Collision detected - reflect particle
                final_pos = self.reflect_from_obstacle(old_pos, final_pos, obstacle);
                
                // Count collision
                if let Some(count) = self.statistics.collision_counts.get_mut(particle_id) {
                    *count += 1;
                }
            }
        }
        
        Ok(final_pos)
    }
    
    /// Reflect particle from obstacle
    fn reflect_from_obstacle(&self, old_pos: (f64, f64), new_pos: (f64, f64), 
                           obstacle: &Obstacle) -> (f64, f64) {
        let dx = new_pos.0 - obstacle.position.0;
        let dy = new_pos.1 - obstacle.position.1;
        let distance = (dx * dx + dy * dy).sqrt();
        
        if distance == 0.0 {
            return old_pos;
        }
        
        let nx = dx / distance;
        let ny = dy / distance;
        
        let reflected_x = obstacle.position.0 + nx * obstacle.radius * 1.01;
        let reflected_y = obstacle.position.1 + ny * obstacle.radius * 1.01;
        
        (reflected_x, reflected_y)
    }
    
    /// Calculate Euclidean distance between two points
    fn calculate_distance(&self, pos1: (f64, f64), pos2: (f64, f64)) -> f64 {
        let dx = pos1.0 - pos2.0;
        let dy = pos1.1 - pos2.1;
        (dx * dx + dy * dy).sqrt()
    }
    
    /// Update statistics for a particle
    fn update_particle_statistics(&mut self, particle_id: &str, position: (f64, f64)) {
        if let Some(trajectory) = self.statistics.trajectories.get_mut(particle_id) {
            trajectory.push(position);
            
            if trajectory.len() > 10000 {
                trajectory.remove(0);
            }
        }
        
        self.update_msd(particle_id);
    }
    
    /// Update mean square displacement
    fn update_msd(&mut self, particle_id: &str) {
        if let Some(trajectory) = self.statistics.trajectories.get(particle_id) {
            if trajectory.len() < 2 {
                return;
            }
            
            let start_pos = trajectory[0];
            let current_pos = trajectory[trajectory.len() - 1];
            
            let dx = current_pos.0 - start_pos.0;
            let dy = current_pos.1 - start_pos.1;
            let msd = dx * dx + dy * dy;
            
            if let Some(msd_vec) = self.statistics.msd.get_mut(particle_id) {
                msd_vec.push(msd);
            }
        }
    }
    
    /// Get particle by ID
    pub fn get_particle(&self, particle_id: &str) -> Option<&DiffusingParticle> {
        self.particles.get(particle_id)
    }
}

impl Default for DiffusionParameters {
    fn default() -> Self {
        Self {
            timestep: 1e-6, // 1 μs
            temperature: PHYSIOLOGICAL_TEMPERATURE,
            viscosity: 1e-3, // Water viscosity at 37°C
            random_seed: 42,
            atp_transport_enabled: true,
        }
    }
}

impl DiffusionStatistics {
    pub fn new() -> Self {
        Self {
            msd: HashMap::new(),
            calculated_diffusion_coeffs: HashMap::new(),
            trajectories: HashMap::new(),
            collision_counts: HashMap::new(),
        }
    }
}

impl DiffusingParticle {
    /// Create a new diffusing particle
    pub fn new(id: String, position: (f64, f64), particle_type: ParticleType) -> Self {
        let diffusion_coefficient = Self::calculate_default_diffusion_coefficient(&particle_type);
        let radius = Self::calculate_particle_radius(&particle_type);
        
        Self {
            id,
            position,
            velocity: (0.0, 0.0),
            diffusion_coefficient,
            radius,
            particle_type,
            diffusion_mode: DiffusionType::Free,
        }
    }
    
    /// Calculate default diffusion coefficient based on particle type
    fn calculate_default_diffusion_coefficient(particle_type: &ParticleType) -> f64 {
        match particle_type {
            ParticleType::Lipid(_) => 1e12, // nm²/s
            ParticleType::Protein(protein_type) => {
                match protein_type {
                    ProteinType::NaKATPase => 0.1e12,
                    ProteinType::CaATPase => 0.15e12,
                    ProteinType::VGSC => 0.5e12,
                    ProteinType::VGKC => 0.8e12,
                    ProteinType::VGCC => 0.3e12,
                    _ => 0.2e12,
                }
            },
            ParticleType::Ion(_) => 1e12,
            ParticleType::Other(_) => 0.5e12,
        }
    }
    
    /// Calculate particle radius
    fn calculate_particle_radius(particle_type: &ParticleType) -> f64 {
        match particle_type {
            ParticleType::Lipid(_) => 0.5, // nm
            ParticleType::Protein(protein_type) => {
                match protein_type {
                    ProteinType::NaKATPase => 3.0,
                    ProteinType::CaATPase => 2.5,
                    ProteinType::VGSC => 2.0,
                    ProteinType::VGKC => 1.5,
                    ProteinType::VGCC => 2.0,
                    _ => 2.0,
                }
            },
            ParticleType::Ion(_) => 0.1,
            ParticleType::Other(_) => 1.0,
        }
    }
}   