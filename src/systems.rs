//! Systems Level - Macro-scale physiological flows and gradients
//!
//! This module implements the highest level of the membrane dynamics hierarchy,
//! handling system-wide flows and gradients that affect membrane function:
//! - Oxygen uptake and metabolic flows
//! - Temperature gradients and thermal regulation
//! - Bulk ionic transport and osmotic balance
//! - Fluid dynamics and perfusion
//! - Homeostatic control and feedback loops

use crate::{
    constants::*,
    error::{MembraneError, Result},
    circuit_interface::CircuitInterface,
    types::*,
};
use std::collections::HashMap;

/// Systems-level physiological state and control
#[derive(Debug, Clone)]
pub struct SystemsLevel {
    /// Metabolic flows and energy balance
    pub metabolism: MetabolicFlows,
    /// Thermal dynamics and temperature control
    pub thermal: ThermalDynamics,
    /// Ionic transport and osmotic regulation
    pub ionic_transport: IonicTransport,
    /// Fluid dynamics and perfusion
    pub fluid_dynamics: FluidDynamics,
    /// Homeostatic control systems
    pub homeostasis: HomeostaticControl,
    /// System-wide ATP pool and energy distribution
    pub energy_distribution: EnergyDistribution,
    /// Environmental conditions
    pub environment: EnvironmentalConditions,
}

/// Metabolic flows and respiratory exchange
#[derive(Debug, Clone)]
pub struct MetabolicFlows {
    /// Oxygen uptake rate (mol/s)
    pub oxygen_uptake: f64,
    /// CO₂ production rate (mol/s)
    pub co2_production: f64,
    /// Glucose consumption rate (mol/s)
    pub glucose_consumption: f64,
    /// ATP production rate (mol/s)
    pub atp_production: f64,
    /// Metabolic rate (W)
    pub metabolic_rate: f64,
    /// Respiratory quotient (CO₂/O₂)
    pub respiratory_quotient: f64,
    /// Lactate production (anaerobic metabolism)
    pub lactate_production: f64,
}

/// Thermal dynamics and temperature regulation
#[derive(Debug, Clone)]
pub struct ThermalDynamics {
    /// Core body temperature (K)
    pub core_temperature: f64,
    /// Local tissue temperature (K)
    pub local_temperature: f64,
    /// Temperature gradient (K/m)
    pub temperature_gradient: f64,
    /// Heat production rate (W)
    pub heat_production: f64,
    /// Heat loss rate (W)
    pub heat_loss: f64,
    /// Thermal conductivity (W/m·K)
    pub thermal_conductivity: f64,
    /// Blood flow thermal transport (W)
    pub convective_heat_transport: f64,
}

/// Bulk ionic transport and osmotic regulation
#[derive(Debug, Clone)]
pub struct IonicTransport {
    /// Bulk ion fluxes (mol/s per ion type)
    pub bulk_fluxes: HashMap<IonType, f64>,
    /// Osmotic pressure gradients (Pa)
    pub osmotic_gradients: Vec<f64>,
    /// Total osmolarity (osmol/L)
    pub total_osmolarity: f64,
    /// Water flux (m³/s)
    pub water_flux: f64,
    /// Electrochemical potential differences (J/mol)
    pub electrochemical_potentials: HashMap<IonType, f64>,
    /// Ion pump coordination (system-wide)
    pub pump_coordination: PumpCoordination,
}

/// Fluid dynamics and perfusion
#[derive(Debug, Clone)]
pub struct FluidDynamics {
    /// Blood flow rate (m³/s)
    pub blood_flow: f64,
    /// Lymphatic flow rate (m³/s)
    pub lymphatic_flow: f64,
    /// Interstitial fluid pressure (Pa)
    pub interstitial_pressure: f64,
    /// Capillary pressure (Pa)
    pub capillary_pressure: f64,
    /// Permeability coefficients
    pub permeability: PermeabilityCoefficients,
    /// Fluid velocity field
    pub velocity_field: VelocityField,
}

/// Homeostatic control and regulatory mechanisms
#[derive(Debug, Clone)]
pub struct HomeostaticControl {
    /// pH regulation system
    pub ph_control: PhRegulation,
    /// Calcium homeostasis
    pub calcium_control: CalciumHomeostasis,
    /// ATP/ADP ratio control
    pub energy_control: EnergyHomeostasis,
    /// Temperature regulation
    pub thermal_control: ThermalRegulation,
    /// Feedback control gains
    pub control_gains: ControlGains,
}

/// System-wide energy distribution
#[derive(Debug, Clone)]
pub struct EnergyDistribution {
    /// Total ATP pool (mol)
    pub total_atp_pool: f64,
    /// ATP distribution by tissue type
    pub atp_allocation: HashMap<String, f64>,
    /// Energy demand by process
    pub energy_demands: HashMap<String, f64>,
    /// Energy supply rates
    pub energy_supplies: HashMap<String, f64>,
    /// Energy efficiency metrics
    pub system_efficiency: f64,
}

/// Environmental conditions affecting the system
#[derive(Debug, Clone)]
pub struct EnvironmentalConditions {
    /// Ambient temperature (K)
    pub ambient_temperature: f64,
    /// Atmospheric pressure (Pa)
    pub atmospheric_pressure: f64,
    /// Oxygen partial pressure (Pa)
    pub oxygen_partial_pressure: f64,
    /// CO₂ partial pressure (Pa)
    pub co2_partial_pressure: f64,
    /// Humidity (0.0 to 1.0)
    pub humidity: f64,
    /// External ion concentrations
    pub external_ion_concentrations: HashMap<IonType, f64>,
}

/// Coordination of ion pumps across the system
#[derive(Debug, Clone)]
pub struct PumpCoordination {
    /// Global pump synchronization factor
    pub synchronization: f64,
    /// Pump efficiency coordination
    pub efficiency_coordination: f64,
    /// ATP allocation to pumps
    pub pump_atp_allocation: f64,
    /// Pump activity modulation
    pub activity_modulation: HashMap<ProteinType, f64>,
}

/// Permeability coefficients for different substances
#[derive(Debug, Clone)]
pub struct PermeabilityCoefficients {
    /// Water permeability (m/s)
    pub water: f64,
    /// Ion permeabilities (m/s per ion type)
    pub ions: HashMap<IonType, f64>,
    /// Gas permeabilities (m/s)
    pub oxygen: f64,
    pub carbon_dioxide: f64,
    /// Glucose permeability (m/s)
    pub glucose: f64,
}

/// Velocity field for fluid transport
#[derive(Debug, Clone)]
pub struct VelocityField {
    /// Bulk fluid velocity (m/s)
    pub bulk_velocity: (f64, f64, f64),
    /// Convective transport coefficients
    pub convection_coeffs: HashMap<String, f64>,
    /// Diffusive transport coefficients
    pub diffusion_coeffs: HashMap<String, f64>,
}

/// pH regulation system
#[derive(Debug, Clone)]
pub struct PhRegulation {
    /// Current pH
    pub current_ph: f64,
    /// Target pH
    pub target_ph: f64,
    /// Buffer capacity (mol/pH unit)
    pub buffer_capacity: f64,
    /// CO₂/HCO₃⁻ buffer system
    pub bicarbonate_buffer: f64,
    /// Phosphate buffer system
    pub phosphate_buffer: f64,
    /// Protein buffer contribution
    pub protein_buffer: f64,
}

/// Calcium homeostasis control
#[derive(Debug, Clone)]
pub struct CalciumHomeostasis {
    /// Free Ca²⁺ concentration (M)
    pub free_calcium: f64,
    /// Bound Ca²⁺ (to proteins, buffers)
    pub bound_calcium: f64,
    /// Calcium buffering capacity
    pub buffering_capacity: f64,
    /// Ca²⁺ pump activity
    pub pump_activity: f64,
    /// Calcium release from stores
    pub store_release: f64,
}

/// Energy homeostasis control
#[derive(Debug, Clone)]
pub struct EnergyHomeostasis {
    /// ATP/ADP ratio
    pub atp_adp_ratio: f64,
    /// Energy charge ((ATP + 0.5*ADP)/(ATP + ADP + AMP))
    pub energy_charge: f64,
    /// Phosphocreatine buffer
    pub phosphocreatine: f64,
    /// Glycogen stores
    pub glycogen_stores: f64,
    /// Metabolic flux control
    pub flux_control: f64,
}

/// Thermal regulation system
#[derive(Debug, Clone)]
pub struct ThermalRegulation {
    /// Thermoregulatory setpoint (K)
    pub setpoint: f64,
    /// Heat production control
    pub heat_production_control: f64,
    /// Heat loss control (vasodilation/constriction)
    pub heat_loss_control: f64,
    /// Behavioral thermoregulation
    pub behavioral_control: f64,
}

/// Control system gains and parameters
#[derive(Debug, Clone)]
pub struct ControlGains {
    /// Proportional gains
    pub proportional: HashMap<String, f64>,
    /// Integral gains
    pub integral: HashMap<String, f64>,
    /// Derivative gains
    pub derivative: HashMap<String, f64>,
    /// Control time constants
    pub time_constants: HashMap<String, f64>,
}

impl SystemsLevel {
    /// Create new systems level with physiological defaults
    pub fn new(environment: EnvironmentalConditions) -> Self {
        Self {
            metabolism: MetabolicFlows::physiological(),
            thermal: ThermalDynamics::physiological(environment.ambient_temperature),
            ionic_transport: IonicTransport::physiological(),
            fluid_dynamics: FluidDynamics::physiological(),
            homeostasis: HomeostaticControl::physiological(),
            energy_distribution: EnergyDistribution::physiological(),
            environment,
        }
    }
    
    /// Update systems level state based on lower-level circuit interfaces
    pub fn update(&mut self, circuit_interfaces: &[CircuitInterface], dt: f64) -> Result<()> {
        // Update metabolic flows based on ATP consumption
        self.update_metabolism(circuit_interfaces, dt)?;
        
        // Update thermal dynamics
        self.update_thermal_dynamics(dt)?;
        
        // Update ionic transport
        self.update_ionic_transport(circuit_interfaces, dt)?;
        
        // Update fluid dynamics
        self.update_fluid_dynamics(dt)?;
        
        // Apply homeostatic control
        self.apply_homeostatic_control(dt)?;
        
        // Redistribute energy
        self.update_energy_distribution(circuit_interfaces, dt)?;
        
        Ok(())
    }
    
    /// Update metabolic flows based on ATP demand
    fn update_metabolism(&mut self, circuit_interfaces: &[CircuitInterface], dt: f64) -> Result<()> {
        // Calculate total ATP consumption from all circuit interfaces
        let total_atp_consumption: f64 = circuit_interfaces.iter()
            .map(|ci| ci.atp_dynamics.atp_consumption_rate)
            .sum();
        
        // Update ATP production to match demand (with efficiency)
        let production_efficiency = 0.38; // ~38% efficiency of cellular respiration
        self.metabolism.atp_production = total_atp_consumption / production_efficiency;
        
        // Calculate oxygen consumption from ATP production
        // 6 O₂ + C₆H₁₂O₆ → 6 CO₂ + 6 H₂O + ~38 ATP
        let atp_per_glucose = 38.0;
        let o2_per_glucose = 6.0;
        self.metabolism.glucose_consumption = self.metabolism.atp_production / atp_per_glucose;
        self.metabolism.oxygen_uptake = self.metabolism.glucose_consumption * o2_per_glucose;
        self.metabolism.co2_production = self.metabolism.glucose_consumption * 6.0;
        
        // Calculate respiratory quotient
        if self.metabolism.oxygen_uptake > 0.0 {
            self.metabolism.respiratory_quotient = self.metabolism.co2_production / self.metabolism.oxygen_uptake;
        }
        
        // Calculate metabolic rate (power)
        self.metabolism.metabolic_rate = self.metabolism.atp_production * 
            atp::ATP_HYDROLYSIS_ENERGY.abs() / production_efficiency;
        
        // Handle anaerobic metabolism if oxygen is insufficient
        let oxygen_deficit = self.metabolism.oxygen_uptake - self.get_oxygen_supply_rate();
        if oxygen_deficit > 0.0 {
            // Lactate fermentation: glucose → 2 lactate + 2 ATP
            self.metabolism.lactate_production = oxygen_deficit / 3.0; // Rough approximation
        }
        
        Ok(())
    }
    
    /// Update thermal dynamics
    fn update_thermal_dynamics(&mut self, dt: f64) -> Result<()> {
        // Heat production from metabolism
        let metabolic_heat = self.metabolism.metabolic_rate * 0.62; // ~62% as heat
        self.thermal.heat_production = metabolic_heat;
        
        // Heat loss through various mechanisms
        let temperature_diff = self.thermal.core_temperature - self.environment.ambient_temperature;
        let convective_loss = 10.0 * temperature_diff; // W, simplified
        let radiative_loss = 5.67e-8 * (self.thermal.core_temperature.powi(4) - 
                                       self.environment.ambient_temperature.powi(4));
        self.thermal.heat_loss = convective_loss + radiative_loss;
        
        // Update core temperature using heat balance
        let heat_capacity = 3500.0; // J/K for human body
        let net_heat = self.thermal.heat_production - self.thermal.heat_loss;
        let dtemp_dt = net_heat / heat_capacity;
        
        self.thermal.core_temperature += dtemp_dt * dt;
        
        // Update local temperature (simplified)
        self.thermal.local_temperature = self.thermal.core_temperature - 
            self.thermal.temperature_gradient * 0.01; // 1 cm depth
        
        Ok(())
    }
    
    /// Update ionic transport
    fn update_ionic_transport(&mut self, circuit_interfaces: &[CircuitInterface], dt: f64) -> Result<()> {
        // Aggregate ion fluxes from all circuit interfaces
        for ion_type in [IonType::Sodium, IonType::Potassium, IonType::Calcium, IonType::Chloride] {
            let mut total_flux = 0.0;
            
            for ci in circuit_interfaces {
                // Sum fluxes from all proteins in each circuit interface
                for protein in ci.membrane.proteins.values() {
                    match (&protein.protein_type, ion_type) {
                        (ProteinType::NaKATPase, IonType::Sodium) => {
                            total_flux -= protein.activity * 3.0 * 100.0; // 3 Na+ out per cycle
                        }
                        (ProteinType::NaKATPase, IonType::Potassium) => {
                            total_flux += protein.activity * 2.0 * 100.0; // 2 K+ in per cycle
                        }
                        (ProteinType::CaATPase, IonType::Calcium) => {
                            total_flux -= protein.activity * 1.0 * 50.0; // 1 Ca2+ out per cycle
                        }
                        _ => {}
                    }
                }
            }
            
            self.ionic_transport.bulk_fluxes.insert(ion_type, total_flux);
        }
        
        // Update osmotic gradients
        let total_osmolarity: f64 = self.ionic_transport.bulk_fluxes.values().sum::<f64>().abs();
        self.ionic_transport.total_osmolarity = total_osmolarity * 1e-3; // Convert to osmol/L
        
        // Calculate water flux from osmotic gradients
        let osmotic_pressure_diff = self.ionic_transport.total_osmolarity * 
            GAS_CONSTANT * self.thermal.local_temperature;
        self.ionic_transport.water_flux = osmotic_pressure_diff * 1e-12; // m³/s
        
        Ok(())
    }
    
    /// Update fluid dynamics
    fn update_fluid_dynamics(&mut self, dt: f64) -> Result<()> {
        // Blood flow affects heat and mass transport
        let metabolic_demand_factor = self.metabolism.metabolic_rate / 100.0; // Normalized
        self.fluid_dynamics.blood_flow = 5e-6 * metabolic_demand_factor; // m³/s, scaled with demand
        
        // Lymphatic flow for fluid balance
        self.fluid_dynamics.lymphatic_flow = self.ionic_transport.water_flux * 0.1; // 10% of water flux
        
        // Update pressures
        let hydrostatic_pressure = 101325.0 + 1000.0 * 9.81 * 0.1; // Pa, atmospheric + 10 cm H₂O
        self.fluid_dynamics.capillary_pressure = hydrostatic_pressure;
        self.fluid_dynamics.interstitial_pressure = hydrostatic_pressure - 2000.0; // 20 mmHg lower
        
        // Update convective heat transport
        let blood_heat_capacity = 3600.0; // J/kg·K
        let blood_density = 1060.0; // kg/m³
        self.thermal.convective_heat_transport = self.fluid_dynamics.blood_flow * 
            blood_density * blood_heat_capacity * 
            (self.thermal.core_temperature - self.thermal.local_temperature);
        
        Ok(())
    }
    
    /// Apply homeostatic control mechanisms
    fn apply_homeostatic_control(&mut self, dt: f64) -> Result<()> {
        // pH control
        let ph_error = self.homeostasis.ph_control.target_ph - self.homeostasis.ph_control.current_ph;
        let ph_correction = ph_error * 0.1 * dt; // Proportional control
        self.homeostasis.ph_control.current_ph += ph_correction;
        
        // Temperature control
        let temp_error = self.homeostasis.thermal_control.setpoint - self.thermal.core_temperature;
        let temp_gain = self.homeostasis.control_gains.proportional.get("temperature").copied().unwrap_or(1.0);
        let heat_adjustment = temp_error * temp_gain;
        self.thermal.heat_production += heat_adjustment;
        
        // Calcium homeostasis
        let ca_target = 1e-7; // 100 nM free Ca²⁺
        let ca_error = ca_target - self.homeostasis.calcium_control.free_calcium;
        let ca_correction = ca_error * 0.5 * dt;
        self.homeostasis.calcium_control.free_calcium += ca_correction;
        
        // Energy homeostasis
        let target_energy_charge = 0.9;
        let energy_error = target_energy_charge - self.homeostasis.energy_control.energy_charge;
        let energy_correction = energy_error * 0.2 * dt;
        self.homeostasis.energy_control.energy_charge += energy_correction;
        
        Ok(())
    }
    
    /// Update energy distribution across the system
    fn update_energy_distribution(&mut self, circuit_interfaces: &[CircuitInterface], dt: f64) -> Result<()> {
        // Calculate total ATP demand
        let total_demand: f64 = circuit_interfaces.iter()
            .map(|ci| ci.atp_dynamics.atp_consumption_rate)
            .sum();
        
        // Update ATP allocations
        for (i, ci) in circuit_interfaces.iter().enumerate() {
            let allocation_key = format!("circuit_{}", i);
            let demand_fraction = if total_demand > 0.0 {
                ci.atp_dynamics.atp_consumption_rate / total_demand
            } else {
                1.0 / circuit_interfaces.len() as f64
            };
            self.energy_distribution.atp_allocation.insert(allocation_key, demand_fraction);
        }
        
        // Update system efficiency
        let energy_input = self.metabolism.metabolic_rate;
        let useful_work = total_demand * atp::ATP_HYDROLYSIS_ENERGY.abs();
        self.energy_distribution.system_efficiency = if energy_input > 0.0 {
            useful_work / energy_input
        } else {
            0.0
        };
        
        Ok(())
    }
    
    /// Get oxygen supply rate based on respiratory and circulatory function
    fn get_oxygen_supply_rate(&self) -> f64 {
        // Simplified oxygen supply calculation
        let oxygen_capacity = self.fluid_dynamics.blood_flow * 0.2; // 20% O₂ capacity
        let oxygen_saturation = self.environment.oxygen_partial_pressure / 21000.0; // Normalized to air
        oxygen_capacity * oxygen_saturation.min(1.0)
    }
    
    /// Get systems-level circuit parameters for integration
    pub fn get_systems_circuit_parameters(&self) -> SystemsCircuitParams {
        SystemsCircuitParams {
            metabolic_node: MetabolicNode {
                oxygen_uptake: self.metabolism.oxygen_uptake,
                co2_production: self.metabolism.co2_production,
                atp_production: self.metabolism.atp_production,
                metabolic_impedance: self.metabolism.metabolic_rate.recip(),
            },
            thermal_node: ThermalNode {
                temperature: self.thermal.local_temperature,
                heat_flux: self.thermal.heat_production - self.thermal.heat_loss,
                thermal_resistance: 1.0 / self.thermal.thermal_conductivity,
                thermal_capacitance: 3500.0, // J/K
            },
            ionic_node: IonicNode {
                total_current: self.ionic_transport.bulk_fluxes.values().sum(),
                osmotic_pressure: self.ionic_transport.total_osmolarity * GAS_CONSTANT * self.thermal.local_temperature,
                ionic_conductance: 1e-3, // S, bulk conductance
            },
            fluid_node: FluidNode {
                flow_rate: self.fluid_dynamics.blood_flow,
                pressure: self.fluid_dynamics.capillary_pressure,
                hydraulic_resistance: 1e6, // Pa·s/m³
            },
            control_parameters: ControlParameters {
                homeostatic_gains: self.homeostasis.control_gains.proportional.clone(),
                stability_margins: self.calculate_stability_margins(),
                adaptation_rates: self.calculate_adaptation_rates(),
            },
        }
    }
    
    /// Calculate system stability margins
    fn calculate_stability_margins(&self) -> HashMap<String, f64> {
        let mut margins = HashMap::new();
        
        // Temperature stability
        let temp_margin = 1.0 - (self.thermal.core_temperature - temperature::PHYSIOLOGICAL).abs() / 10.0;
        margins.insert("temperature".to_string(), temp_margin.max(0.0));
        
        // pH stability
        let ph_margin = 1.0 - (self.homeostasis.ph_control.current_ph - 7.4).abs() / 0.5;
        margins.insert("ph".to_string(), ph_margin.max(0.0));
        
        // Energy stability
        margins.insert("energy".to_string(), self.homeostasis.energy_control.energy_charge);
        
        margins
    }
    
    /// Calculate adaptation rates for different systems
    fn calculate_adaptation_rates(&self) -> HashMap<String, f64> {
        let mut rates = HashMap::new();
        
        rates.insert("metabolic".to_string(), 0.1); // 10% per update
        rates.insert("thermal".to_string(), 0.05);  // 5% per update
        rates.insert("ionic".to_string(), 0.2);     // 20% per update
        rates.insert("fluid".to_string(), 0.15);    // 15% per update
        
        rates
    }
}

/// Systems-level circuit parameters for hierarchical integration
#[derive(Debug, Clone)]
pub struct SystemsCircuitParams {
    pub metabolic_node: MetabolicNode,
    pub thermal_node: ThermalNode,
    pub ionic_node: IonicNode,
    pub fluid_node: FluidNode,
    pub control_parameters: ControlParameters,
}

/// Metabolic circuit node
#[derive(Debug, Clone)]
pub struct MetabolicNode {
    pub oxygen_uptake: f64,
    pub co2_production: f64,
    pub atp_production: f64,
    pub metabolic_impedance: f64,
}

/// Thermal circuit node
#[derive(Debug, Clone)]
pub struct ThermalNode {
    pub temperature: f64,
    pub heat_flux: f64,
    pub thermal_resistance: f64,
    pub thermal_capacitance: f64,
}

/// Ionic transport circuit node
#[derive(Debug, Clone)]
pub struct IonicNode {
    pub total_current: f64,
    pub osmotic_pressure: f64,
    pub ionic_conductance: f64,
}

/// Fluid dynamics circuit node
#[derive(Debug, Clone)]
pub struct FluidNode {
    pub flow_rate: f64,
    pub pressure: f64,
    pub hydraulic_resistance: f64,
}

/// Control system parameters
#[derive(Debug, Clone)]
pub struct ControlParameters {
    pub homeostatic_gains: HashMap<String, f64>,
    pub stability_margins: HashMap<String, f64>,
    pub adaptation_rates: HashMap<String, f64>,
}

// Default implementations for physiological states
impl MetabolicFlows {
    fn physiological() -> Self {
        Self {
            oxygen_uptake: 3.5e-6,      // ~3.5 mL/min/kg * 70kg / 60s
            co2_production: 2.8e-6,     // RQ = 0.8
            glucose_consumption: 1e-7,   // mol/s
            atp_production: 3.8e-6,     // mol/s
            metabolic_rate: 100.0,      // W (basal metabolic rate)
            respiratory_quotient: 0.8,
            lactate_production: 0.0,
        }
    }
}

impl ThermalDynamics {
    fn physiological(ambient_temp: f64) -> Self {
        Self {
            core_temperature: PHYSIOLOGICAL_TEMPERATURE,
            local_temperature: PHYSIOLOGICAL_TEMPERATURE - 1.0,
            temperature_gradient: 100.0, // K/m
            heat_production: 100.0,      // W
            heat_loss: 100.0,           // W (balanced)
            thermal_conductivity: 0.5,   // W/m·K
            convective_heat_transport: 50.0, // W
        }
    }
}

impl IonicTransport {
    fn physiological() -> Self {
        let mut bulk_fluxes = HashMap::new();
        bulk_fluxes.insert(IonType::Sodium, -1e-12);    // mol/s, net efflux
        bulk_fluxes.insert(IonType::Potassium, 0.67e-12); // mol/s, net influx
        bulk_fluxes.insert(IonType::Calcium, -1e-15);   // mol/s, small efflux
        bulk_fluxes.insert(IonType::Chloride, 1e-12);   // mol/s, net influx
        
        Self {
            bulk_fluxes,
            osmotic_gradients: vec![0.0; 10],
            total_osmolarity: 0.3, // osmol/L
            water_flux: 1e-12,     // m³/s
            electrochemical_potentials: HashMap::new(),
            pump_coordination: PumpCoordination {
                synchronization: 0.8,
                efficiency_coordination: 0.9,
                pump_atp_allocation: 0.3,
                activity_modulation: HashMap::new(),
            },
        }
    }
}

impl FluidDynamics {
    fn physiological() -> Self {
        Self {
            blood_flow: 5e-6,           // m³/s (~5 mL/s local)
            lymphatic_flow: 1e-7,       // m³/s
            interstitial_pressure: 101325.0 - 2000.0, // Pa
            capillary_pressure: 101325.0 + 3000.0,    // Pa
            permeability: PermeabilityCoefficients::physiological(),
            velocity_field: VelocityField {
                bulk_velocity: (1e-3, 0.0, 0.0), // m/s
                convection_coeffs: HashMap::new(),
                diffusion_coeffs: HashMap::new(),
            },
        }
    }
}

impl PermeabilityCoefficients {
    fn physiological() -> Self {
        let mut ions = HashMap::new();
        ions.insert(IonType::Sodium, 1e-12);
        ions.insert(IonType::Potassium, 2e-12);
        ions.insert(IonType::Calcium, 1e-14);
        ions.insert(IonType::Chloride, 5e-12);
        
        Self {
            water: 2.4e-14,    // m/s
            ions,
            oxygen: 1e-9,      // m/s
            carbon_dioxide: 2e-9, // m/s
            glucose: 1e-11,    // m/s
        }
    }
}

impl HomeostaticControl {
    fn physiological() -> Self {
        let mut proportional = HashMap::new();
        proportional.insert("temperature".to_string(), 10.0);
        proportional.insert("ph".to_string(), 5.0);
        proportional.insert("calcium".to_string(), 100.0);
        proportional.insert("energy".to_string(), 2.0);
        
        Self {
            ph_control: PhRegulation {
                current_ph: 7.4,
                target_ph: 7.4,
                buffer_capacity: 0.03, // mol/pH
                bicarbonate_buffer: 0.024,
                phosphate_buffer: 0.003,
                protein_buffer: 0.003,
            },
            calcium_control: CalciumHomeostasis {
                free_calcium: 1e-7,    // 100 nM
                bound_calcium: 1e-3,   // 1 mM total
                buffering_capacity: 1000.0,
                pump_activity: 0.8,
                store_release: 0.0,
            },
            energy_control: EnergyHomeostasis {
                atp_adp_ratio: 10.0,
                energy_charge: 0.9,
                phosphocreatine: 20e-3, // 20 mM
                glycogen_stores: 100e-3, // 100 mM glucose equiv
                flux_control: 1.0,
            },
            thermal_control: ThermalRegulation {
                setpoint: PHYSIOLOGICAL_TEMPERATURE,
                heat_production_control: 1.0,
                heat_loss_control: 1.0,
                behavioral_control: 0.0,
            },
            control_gains: ControlGains {
                proportional,
                integral: HashMap::new(),
                derivative: HashMap::new(),
                time_constants: HashMap::new(),
            },
        }
    }
}

impl EnergyDistribution {
    fn physiological() -> Self {
        Self {
            total_atp_pool: 5e-3 * 1e-12, // 5 mM in 1 pL volume
            atp_allocation: HashMap::new(),
            energy_demands: HashMap::new(),
            energy_supplies: HashMap::new(),
            system_efficiency: 0.25, // 25% overall efficiency
        }
    }
}

impl Default for EnvironmentalConditions {
    fn default() -> Self {
        let mut external_ions = HashMap::new();
        external_ions.insert(IonType::Sodium, 0.145);    // 145 mM
        external_ions.insert(IonType::Potassium, 0.005); // 5 mM
        external_ions.insert(IonType::Calcium, 0.002);   // 2 mM
        external_ions.insert(IonType::Chloride, 0.110);  // 110 mM
        
        Self {
            ambient_temperature: 298.15,     // 25°C
            atmospheric_pressure: 101325.0,  // Pa
            oxygen_partial_pressure: 21000.0, // Pa (21% of atm)
            co2_partial_pressure: 40.0,      // Pa (0.04% of atm)
            humidity: 0.5,                   // 50% RH
            external_ion_concentrations: external_ions,
        }
    }
} 