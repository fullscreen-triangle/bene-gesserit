//! Protein clustering and spatial organization in membranes
//!
//! This module implements protein clustering dynamics including:
//! - Cluster formation and dissolution
//! - Spatial organization patterns
//! - Cooperative clustering effects
//! - ATP-dependent cluster maintenance

use crate::{
    constants::*,
    error::{MembraneError, Result},
    types::*,
};
use std::collections::HashMap;

/// Protein cluster structure
#[derive(Debug, Clone)]
pub struct ProteinCluster {
    /// Unique cluster identifier
    pub id: String,
    /// Proteins in this cluster
    pub proteins: Vec<String>,
    /// Cluster center position (nm)
    pub center: (f64, f64),
    /// Cluster radius (nm)
    pub radius: f64,
    /// Clustering energy (kJ/mol)
    pub binding_energy: f64,
    /// Cluster stability factor
    pub stability: f64,
    /// Time since cluster formation
    pub age: f64,
}

/// Clustering dynamics manager
#[derive(Debug, Clone)]
pub struct ClusteringManager {
    /// Active protein clusters
    pub clusters: HashMap<String, ProteinCluster>,
    /// Clustering parameters
    pub parameters: ClusteringParameters,
    /// Clustering statistics
    pub statistics: ClusteringStatistics,
}

/// Clustering parameters and thresholds
#[derive(Debug, Clone)]
pub struct ClusteringParameters {
    /// Maximum clustering distance (nm)
    pub max_cluster_distance: f64,
    /// Minimum proteins per cluster
    pub min_cluster_size: usize,
    /// Maximum proteins per cluster
    pub max_cluster_size: usize,
    /// Clustering energy threshold
    pub energy_threshold: f64,
    /// ATP cost per cluster maintenance
    pub atp_cost_per_cluster: f64,
}

/// Clustering statistics
#[derive(Debug, Clone)]
pub struct ClusteringStatistics {
    /// Total clusters formed
    pub total_clusters_formed: u64,
    /// Total clusters dissolved
    pub total_clusters_dissolved: u64,
    /// Average cluster size
    pub avg_cluster_size: f64,
    /// Average cluster lifetime
    pub avg_cluster_lifetime: f64,
    /// Total ATP consumed in clustering
    pub total_atp_consumed: f64,
}

impl ClusteringManager {
    /// Create a new clustering manager
    pub fn new() -> Self {
        Self {
            clusters: HashMap::new(),
            parameters: ClusteringParameters::default(),
            statistics: ClusteringStatistics::new(),
        }
    }
    
    /// Update clustering dynamics
    pub fn update_clustering(&mut self, proteins: &HashMap<String, crate::molecular::MembraneProtein>,
                           dt: f64, atp_pool: &mut AtpMeasurement) -> Result<()> {
        // Age existing clusters
        for cluster in self.clusters.values_mut() {
            cluster.age += dt;
        }
        
        // Check for new cluster formation
        self.detect_new_clusters(proteins)?;
        
        // Update existing clusters
        self.update_existing_clusters(proteins, dt, atp_pool)?;
        
        // Check for cluster dissolution
        self.check_cluster_dissolution(proteins)?;
        
        // Update statistics
        self.update_statistics();
        
        Ok(())
    }
    
    /// Detect formation of new clusters
    fn detect_new_clusters(&mut self, proteins: &HashMap<String, crate::molecular::MembraneProtein>) -> Result<()> {
        let mut unclustered_proteins: Vec<_> = proteins.keys()
            .filter(|id| !self.is_protein_clustered(id))
            .collect();
        
        while unclustered_proteins.len() >= self.parameters.min_cluster_size {
            if let Some(new_cluster) = self.form_cluster_from_proteins(&unclustered_proteins, proteins)? {
                // Remove clustered proteins from unclustered list
                unclustered_proteins.retain(|id| !new_cluster.proteins.contains(&id.to_string()));
                
                // Add new cluster
                self.clusters.insert(new_cluster.id.clone(), new_cluster);
                self.statistics.total_clusters_formed += 1;
            } else {
                break; // No more clusters can be formed
            }
        }
        
        Ok(())
    }
    
    /// Form a cluster from available proteins
    fn form_cluster_from_proteins(&self, available_proteins: &[&String],
                                proteins: &HashMap<String, crate::molecular::MembraneProtein>) -> Result<Option<ProteinCluster>> {
        if available_proteins.is_empty() {
            return Ok(None);
        }
        
        // Start with the first protein
        let seed_protein_id = available_proteins[0];
        let seed_protein = proteins.get(seed_protein_id)
            .ok_or_else(|| MembraneError::InvalidParameter(format!("Protein {} not found", seed_protein_id)))?;
        
        let mut cluster_proteins = vec![seed_protein_id.to_string()];
        let mut cluster_center = seed_protein.position;
        
        // Find nearby proteins to add to cluster
        for &protein_id in available_proteins.iter().skip(1) {
            if let Some(protein) = proteins.get(protein_id) {
                let distance = self.calculate_distance(cluster_center, protein.position);
                
                if distance <= self.parameters.max_cluster_distance && 
                   cluster_proteins.len() < self.parameters.max_cluster_size {
                    
                    // Check if this protein type can cluster with existing proteins
                    if self.can_cluster_together(&cluster_proteins, protein_id, proteins)? {
                        cluster_proteins.push(protein_id.to_string());
                        
                        // Update cluster center
                        cluster_center = self.calculate_cluster_center(&cluster_proteins, proteins)?;
                    }
                }
            }
        }
        
        // Only form cluster if minimum size is met
        if cluster_proteins.len() >= self.parameters.min_cluster_size {
            let cluster_id = format!("cluster_{}", self.clusters.len());
            let radius = self.calculate_cluster_radius(&cluster_proteins, proteins, cluster_center)?;
            let binding_energy = self.calculate_cluster_binding_energy(&cluster_proteins, proteins)?;
            
            Ok(Some(ProteinCluster {
                id: cluster_id,
                proteins: cluster_proteins,
                center: cluster_center,
                radius,
                binding_energy,
                stability: self.calculate_cluster_stability(binding_energy),
                age: 0.0,
            }))
        } else {
            Ok(None)
        }
    }
    
    /// Check if protein types can cluster together
    fn can_cluster_together(&self, existing_proteins: &[String], new_protein_id: &str,
                          proteins: &HashMap<String, crate::molecular::MembraneProtein>) -> Result<bool> {
        let new_protein = proteins.get(new_protein_id)
            .ok_or_else(|| MembraneError::InvalidParameter(format!("Protein {} not found", new_protein_id)))?;
        
        // Simple clustering rules based on protein types
        for existing_id in existing_proteins {
            if let Some(existing_protein) = proteins.get(existing_id) {
                match (&existing_protein.protein_type, &new_protein.protein_type) {
                    // Same protein types can always cluster
                    (a, b) if a == b => return Ok(true),
                    
                    // Compatible protein types
                    (ProteinType::NaKATPase, ProteinType::CaATPase) |
                    (ProteinType::CaATPase, ProteinType::NaKATPase) => return Ok(true),
                    
                    (ProteinType::VGSC, ProteinType::VGKC) |
                    (ProteinType::VGKC, ProteinType::VGSC) => return Ok(true),
                    
                    // Incompatible types
                    _ => continue,
                }
            }
        }
        
        Ok(false)
    }
    
    /// Calculate cluster center position
    fn calculate_cluster_center(&self, protein_ids: &[String],
                              proteins: &HashMap<String, crate::molecular::MembraneProtein>) -> Result<(f64, f64)> {
        let mut sum_x = 0.0;
        let mut sum_y = 0.0;
        let mut count = 0;
        
        for protein_id in protein_ids {
            if let Some(protein) = proteins.get(protein_id) {
                sum_x += protein.position.0;
                sum_y += protein.position.1;
                count += 1;
            }
        }
        
        if count > 0 {
            Ok((sum_x / count as f64, sum_y / count as f64))
        } else {
            Err(MembraneError::InvalidParameter("No proteins found for cluster center calculation".to_string()))
        }
    }
    
    /// Calculate cluster radius
    fn calculate_cluster_radius(&self, protein_ids: &[String],
                              proteins: &HashMap<String, crate::molecular::MembraneProtein>,
                              center: (f64, f64)) -> Result<f64> {
        let mut max_distance = 0.0;
        
        for protein_id in protein_ids {
            if let Some(protein) = proteins.get(protein_id) {
                let distance = self.calculate_distance(center, protein.position);
                max_distance = max_distance.max(distance);
            }
        }
        
        Ok(max_distance)
    }
    
    /// Calculate cluster binding energy
    fn calculate_cluster_binding_energy(&self, protein_ids: &[String],
                                      proteins: &HashMap<String, crate::molecular::MembraneProtein>) -> Result<f64> {
        let mut total_energy = 0.0;
        
        // Pairwise interactions within cluster
        for i in 0..protein_ids.len() {
            for j in (i+1)..protein_ids.len() {
                if let (Some(protein1), Some(protein2)) = (proteins.get(&protein_ids[i]), proteins.get(&protein_ids[j])) {
                    let distance = self.calculate_distance(protein1.position, protein2.position);
                    let interaction_energy = self.calculate_pairwise_interaction_energy(
                        &protein1.protein_type, &protein2.protein_type, distance
                    );
                    total_energy += interaction_energy;
                }
            }
        }
        
        Ok(total_energy)
    }
    
    /// Calculate pairwise interaction energy between proteins
    fn calculate_pairwise_interaction_energy(&self, type1: &ProteinType, type2: &ProteinType, distance: f64) -> f64 {
        // Simplified interaction model
        let r0 = 5e-9; // 5 nm characteristic distance
        let base_energy = match (type1, type2) {
            (ProteinType::NaKATPase, ProteinType::NaKATPase) => -10.0, // Strong self-interaction
            (ProteinType::CaATPase, ProteinType::CaATPase) => -8.0,
            (ProteinType::VGSC, ProteinType::VGKC) => -5.0, // Moderate cross-interaction
            _ => -2.0, // Weak default interaction
        };
        
        base_energy * (-distance / r0).exp()
    }
    
    /// Calculate cluster stability
    fn calculate_cluster_stability(&self, binding_energy: f64) -> f64 {
        // Stability based on binding energy relative to thermal energy
        let thermal_energy = BOLTZMANN_CONSTANT * PHYSIOLOGICAL_TEMPERATURE;
        (binding_energy.abs() / thermal_energy).min(1.0)
    }
    
    /// Update existing clusters
    fn update_existing_clusters(&mut self, proteins: &HashMap<String, crate::molecular::MembraneProtein>,
                              dt: f64, atp_pool: &mut AtpMeasurement) -> Result<()> {
        let mut clusters_to_update = Vec::new();
        
        for (cluster_id, cluster) in &self.clusters {
            // Check ATP requirement for cluster maintenance
            let atp_required = self.parameters.atp_cost_per_cluster * dt;
            let atp_molecules_required = (atp_required / ATP_ENERGY_PER_MOLECULE) as u64;
            
            if atp_pool.molecules >= atp_molecules_required {
                atp_pool.molecules -= atp_molecules_required;
                self.statistics.total_atp_consumed += atp_required;
                
                clusters_to_update.push(cluster_id.clone());
            }
        }
        
        // Update cluster properties
        for cluster_id in clusters_to_update {
            if let Some(cluster) = self.clusters.get_mut(&cluster_id) {
                // Recalculate cluster center and radius
                cluster.center = self.calculate_cluster_center(&cluster.proteins, proteins)?;
                cluster.radius = self.calculate_cluster_radius(&cluster.proteins, proteins, cluster.center)?;
                cluster.binding_energy = self.calculate_cluster_binding_energy(&cluster.proteins, proteins)?;
                cluster.stability = self.calculate_cluster_stability(cluster.binding_energy);
            }
        }
        
        Ok(())
    }
    
    /// Check for cluster dissolution
    fn check_cluster_dissolution(&mut self, proteins: &HashMap<String, crate::molecular::MembraneProtein>) -> Result<()> {
        let mut clusters_to_dissolve = Vec::new();
        
        for (cluster_id, cluster) in &self.clusters {
            // Check dissolution criteria
            if self.should_dissolve_cluster(cluster, proteins)? {
                clusters_to_dissolve.push(cluster_id.clone());
            }
        }
        
        // Remove dissolved clusters
        for cluster_id in clusters_to_dissolve {
            self.clusters.remove(&cluster_id);
            self.statistics.total_clusters_dissolved += 1;
        }
        
        Ok(())
    }
    
    /// Determine if cluster should dissolve
    fn should_dissolve_cluster(&self, cluster: &ProteinCluster,
                             proteins: &HashMap<String, crate::molecular::MembraneProtein>) -> Result<bool> {
        // Dissolve if stability is too low
        if cluster.stability < 0.1 {
            return Ok(true);
        }
        
        // Dissolve if cluster is too spread out
        if cluster.radius > 2.0 * self.parameters.max_cluster_distance {
            return Ok(true);
        }
        
        // Dissolve if too few proteins remain
        if cluster.proteins.len() < self.parameters.min_cluster_size {
            return Ok(true);
        }
        
        // Dissolve if cluster is too old (optional aging mechanism)
        if cluster.age > 3600.0 { // 1 hour
            return Ok(true);
        }
        
        Ok(false)
    }
    
    /// Check if protein is already in a cluster
    fn is_protein_clustered(&self, protein_id: &str) -> bool {
        self.clusters.values().any(|cluster| cluster.proteins.contains(&protein_id.to_string()))
    }
    
    /// Calculate distance between two points
    fn calculate_distance(&self, pos1: (f64, f64), pos2: (f64, f64)) -> f64 {
        let dx = pos1.0 - pos2.0;
        let dy = pos1.1 - pos2.1;
        (dx * dx + dy * dy).sqrt()
    }
    
    /// Update clustering statistics
    fn update_statistics(&mut self) {
        if !self.clusters.is_empty() {
            let total_proteins: usize = self.clusters.values().map(|c| c.proteins.len()).sum();
            self.statistics.avg_cluster_size = total_proteins as f64 / self.clusters.len() as f64;
            
            let total_age: f64 = self.clusters.values().map(|c| c.age).sum();
            self.statistics.avg_cluster_lifetime = total_age / self.clusters.len() as f64;
        }
    }
    
    /// Get cluster containing specific protein
    pub fn get_protein_cluster(&self, protein_id: &str) -> Option<&ProteinCluster> {
        self.clusters.values().find(|cluster| cluster.proteins.contains(&protein_id.to_string()))
    }
    
    /// Get all clusters
    pub fn get_all_clusters(&self) -> &HashMap<String, ProteinCluster> {
        &self.clusters
    }
}

impl Default for ClusteringParameters {
    fn default() -> Self {
        Self {
            max_cluster_distance: 20e-9, // 20 nm
            min_cluster_size: 3,
            max_cluster_size: 50,
            energy_threshold: -5.0, // kJ/mol
            atp_cost_per_cluster: 1e-21, // J/s per cluster
        }
    }
}

impl ClusteringStatistics {
    pub fn new() -> Self {
        Self {
            total_clusters_formed: 0,
            total_clusters_dissolved: 0,
            avg_cluster_size: 0.0,
            avg_cluster_lifetime: 0.0,
            total_atp_consumed: 0.0,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_clustering_manager_creation() {
        let manager = ClusteringManager::new();
        assert!(manager.clusters.is_empty());
        assert_eq!(manager.parameters.min_cluster_size, 3);
    }
    
    #[test]
    fn test_distance_calculation() {
        let manager = ClusteringManager::new();
        let distance = manager.calculate_distance((0.0, 0.0), (3.0, 4.0));
        assert!((distance - 5.0).abs() < 1e-10);
    }
    
    #[test]
    fn test_cluster_stability_calculation() {
        let manager = ClusteringManager::new();
        let stability = manager.calculate_cluster_stability(-50.0);
        assert!(stability > 0.0 && stability <= 1.0);
    }
}
