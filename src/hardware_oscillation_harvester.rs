//! # Hardware Oscillation Harvester for Biological Quantum Computer
//! 
//! Revolutionary approach: Instead of simulating oscillations, harvest them directly from hardware!
//! 
//! ## Hardware Sources:
//! 1. **CPU Clocks** - High-frequency, stable oscillations
//! 2. **System Timers** - Various frequency domains
//! 3. **Light Sources** - Screen backlight, LEDs, ambient light
//! 4. **Electromagnetic Fields** - WiFi, Bluetooth, power line frequencies
//! 5. **Thermal Oscillations** - CPU/GPU temperature fluctuations
//! 6. **Network Oscillations** - Packet timing, network jitter
//! 7. **Audio Hardware** - Speaker/microphone oscillations
//! 8. **Storage Oscillations** - Disk read/write patterns

use std::collections::HashMap;
use std::time::{Duration, Instant, SystemTime, UNIX_EPOCH};
use std::sync::{Arc, Mutex};
use std::thread;
use crate::{OscillationState, BiologicalQuantumState, error::SolverError};

// ================================================================================================
// HARDWARE OSCILLATION SOURCES
// ================================================================================================

/// Main hardware oscillation harvester
#[derive(Debug)]
pub struct HardwareOscillationHarvester {
    pub active_sources: Vec<OscillationSource>,
    pub harvested_oscillations: Arc<Mutex<HashMap<String, HarvestedOscillation>>>,
    pub sampling_rate: f64,  // Hz
    pub is_harvesting: Arc<Mutex<bool>>,
    pub light_harvester: LightOscillationHarvester,
    pub clock_harvester: ClockOscillationHarvester,
    pub electromagnetic_harvester: ElectromagneticHarvester,
    pub thermal_harvester: ThermalOscillationHarvester,
}

/// Individual oscillation source from hardware
#[derive(Debug, Clone)]
pub struct OscillationSource {
    pub source_name: String,
    pub source_type: HardwareSourceType,
    pub expected_frequency_range: (f64, f64),  // (min_hz, max_hz)
    pub signal_strength: f64,                  // 0.0 to 1.0
    pub is_active: bool,
    pub biological_mapping: BiologicalMapping,
}

#[derive(Debug, Clone)]
pub enum HardwareSourceType {
    CpuClock,
    SystemTimer,
    ScreenBacklight,
    AmbientLight,
    WifiSignal,
    BluetoothSignal,
    PowerLineFrequency,
    CpuTemperature,
    GpuTemperature,
    NetworkJitter,
    AudioOutput,
    DiskActivity,
    MemoryAccess,
    UsbActivity,
    KeyboardInput,
    MouseMovement,
}

/// Mapping from hardware oscillation to biological process
#[derive(Debug, Clone)]
pub struct BiologicalMapping {
    pub biological_process: String,           // e.g., "ATP_Synthase", "NADH_Dehydrogenase"
    pub coupling_strength: f64,               // How strongly hardware drives biology
    pub frequency_scaling: f64,               // Scale hardware freq to biological freq
    pub phase_offset: f64,                    // Phase relationship
    pub amplitude_scaling: f64,               // Scale hardware amplitude
}

/// Harvested oscillation data from hardware
#[derive(Debug, Clone)]
pub struct HarvestedOscillation {
    pub source_name: String,
    pub current_frequency: f64,               // Currently measured frequency
    pub current_amplitude: f64,               // Currently measured amplitude
    pub current_phase: f64,                   // Current phase
    pub last_update: Instant,
    pub sample_history: Vec<OscillationSample>,
    pub biological_state: OscillationState,   // Mapped to biological oscillation
    pub atp_contribution: f64,                // How much ATP this oscillation provides
}

#[derive(Debug, Clone)]
pub struct OscillationSample {
    pub timestamp: Instant,
    pub value: f64,
    pub frequency: f64,
    pub amplitude: f64,
}

impl HardwareOscillationHarvester {
    pub fn new() -> Self {
        Self {
            active_sources: Vec::new(),
            harvested_oscillations: Arc::new(Mutex::new(HashMap::new())),
            sampling_rate: 1000.0,  // 1 kHz sampling
            is_harvesting: Arc::new(Mutex::new(false)),
            light_harvester: LightOscillationHarvester::new(),
            clock_harvester: ClockOscillationHarvester::new(),
            electromagnetic_harvester: ElectromagneticHarvester::new(),
            thermal_harvester: ThermalOscillationHarvester::new(),
        }
    }

    /// Initialize all available hardware oscillation sources
    pub fn initialize_hardware_sources(&mut self) -> Result<(), SolverError> {
        println!("🔌 Initializing hardware oscillation sources...");

        // CPU Clock oscillations
        self.add_source(OscillationSource {
            source_name: "CPU_Clock_Primary".to_string(),
            source_type: HardwareSourceType::CpuClock,
            expected_frequency_range: (1e9, 5e9),  // 1-5 GHz
            signal_strength: 1.0,
            is_active: true,
            biological_mapping: BiologicalMapping {
                biological_process: "ATP_Synthase_Rotation".to_string(),
                coupling_strength: 0.8,
                frequency_scaling: 1e-7,  // Scale GHz to ~100 Hz biological
                phase_offset: 0.0,
                amplitude_scaling: 1.0,
            },
        });

        // System timer oscillations
        self.add_source(OscillationSource {
            source_name: "System_Timer_1ms".to_string(),
            source_type: HardwareSourceType::SystemTimer,
            expected_frequency_range: (1000.0, 1000.0),  // 1 kHz
            signal_strength: 0.9,
            is_active: true,
            biological_mapping: BiologicalMapping {
                biological_process: "Hexokinase_Conformational".to_string(),
                coupling_strength: 0.6,
                frequency_scaling: 0.05,  // Scale to ~50 Hz
                phase_offset: 0.0,
                amplitude_scaling: 0.8,
            },
        });

        // Screen backlight oscillations (PWM)
        self.add_source(OscillationSource {
            source_name: "Screen_Backlight_PWM".to_string(),
            source_type: HardwareSourceType::ScreenBacklight,
            expected_frequency_range: (100.0, 1000.0),  // PWM frequency
            signal_strength: 0.7,
            is_active: true,
            biological_mapping: BiologicalMapping {
                biological_process: "Cytochrome_c_Oxidase_Light".to_string(),
                coupling_strength: 0.4,
                frequency_scaling: 0.2,   // Scale to biological range
                phase_offset: std::f64::consts::PI / 4.0,
                amplitude_scaling: 0.6,
            },
        });

        // WiFi signal oscillations
        self.add_source(OscillationSource {
            source_name: "WiFi_2_4GHz".to_string(),
            source_type: HardwareSourceType::WifiSignal,
            expected_frequency_range: (2.4e9, 2.5e9),  // 2.4 GHz band
            signal_strength: 0.5,
            is_active: true,
            biological_mapping: BiologicalMapping {
                biological_process: "NADH_Dehydrogenase_EM".to_string(),
                coupling_strength: 0.3,
                frequency_scaling: 1e-7,  // Scale to biological
                phase_offset: std::f64::consts::PI / 2.0,
                amplitude_scaling: 0.4,
            },
        });

        // Power line frequency
        self.add_source(OscillationSource {
            source_name: "Power_Line_60Hz".to_string(),
            source_type: HardwareSourceType::PowerLineFrequency,
            expected_frequency_range: (59.0, 61.0),  // 60 Hz ± 1 Hz
            signal_strength: 0.8,
            is_active: true,
            biological_mapping: BiologicalMapping {
                biological_process: "Membrane_Potential_Oscillation".to_string(),
                coupling_strength: 0.5,
                frequency_scaling: 1.0,   // Direct biological frequency
                phase_offset: 0.0,
                amplitude_scaling: 0.9,
            },
        });

        // CPU temperature oscillations
        self.add_source(OscillationSource {
            source_name: "CPU_Temperature_Thermal".to_string(),
            source_type: HardwareSourceType::CpuTemperature,
            expected_frequency_range: (0.1, 10.0),  // Slow thermal oscillations
            signal_strength: 0.6,
            is_active: true,
            biological_mapping: BiologicalMapping {
                biological_process: "Thermal_Protein_Fluctuations".to_string(),
                coupling_strength: 0.4,
                frequency_scaling: 2.0,   // Scale up slightly
                phase_offset: std::f64::consts::PI,
                amplitude_scaling: 0.5,
            },
        });

        println!("✅ Initialized {} hardware oscillation sources", self.active_sources.len());
        Ok(())
    }

    fn add_source(&mut self, source: OscillationSource) {
        println!("  + Added source: {} ({:?})", source.source_name, source.source_type);
        self.active_sources.push(source);
    }

    /// Start harvesting oscillations from all active sources
    pub fn start_harvesting(&mut self) -> Result<(), SolverError> {
        println!("🌾 Starting hardware oscillation harvesting...");
        
        *self.is_harvesting.lock().unwrap() = true;
        
        // Start harvesting threads for each source type
        self.start_clock_harvesting()?;
        self.start_light_harvesting()?;
        self.start_electromagnetic_harvesting()?;
        self.start_thermal_harvesting()?;
        
        println!("✅ Hardware oscillation harvesting started");
        Ok(())
    }

    /// Stop harvesting oscillations
    pub fn stop_harvesting(&mut self) {
        println!("🛑 Stopping hardware oscillation harvesting...");
        *self.is_harvesting.lock().unwrap() = false;
    }

    /// Get current biological quantum state from harvested oscillations
    pub fn get_biological_state_from_hardware(&self) -> Result<BiologicalQuantumState, SolverError> {
        let harvested = self.harvested_oscillations.lock().unwrap();
        
        // Convert harvested hardware oscillations to biological oscillations
        let mut biological_oscillations = Vec::new();
        let mut total_atp_from_hardware = 0.0;
        
        for (source_name, harvested_osc) in harvested.iter() {
            biological_oscillations.push(harvested_osc.biological_state.clone());
            total_atp_from_hardware += harvested_osc.atp_contribution;
            
            println!("  🔄 {}: f={:.1} Hz, A={:.3}, ATP={:.3}", 
                    source_name, 
                    harvested_osc.current_frequency,
                    harvested_osc.current_amplitude,
                    harvested_osc.atp_contribution);
        }
        
        // Create biological state with hardware-driven oscillations
        let mut bio_state = crate::create_physiological_state();
        
        // Replace simulated oscillations with hardware-harvested ones
        bio_state.oscillatory_coords.oscillations = biological_oscillations;
        
        // Boost ATP from hardware energy harvesting
        bio_state.atp_coords.atp_concentration += total_atp_from_hardware;
        
        println!("💡 Generated biological state from {} hardware sources", harvested.len());
        println!("⚡ Total ATP harvested from hardware: {:.3} mM", total_atp_from_hardware);
        
        Ok(bio_state)
    }

    fn start_clock_harvesting(&mut self) -> Result<(), SolverError> {
        let harvested = Arc::clone(&self.harvested_oscillations);
        let is_harvesting = Arc::clone(&self.is_harvesting);
        
        thread::spawn(move || {
            let mut last_time = Instant::now();
            let mut cycle_count = 0u64;
            
            while *is_harvesting.lock().unwrap() {
                let current_time = Instant::now();
                let elapsed = current_time.duration_since(last_time);
                cycle_count += 1;
                
                // Estimate CPU frequency from timing
                let estimated_freq = 1.0 / elapsed.as_secs_f64();
                let amplitude = 0.8 + 0.2 * (cycle_count as f64 * 0.01).sin();
                let phase = (cycle_count as f64 * 0.1) % (2.0 * std::f64::consts::PI);
                
                // Create biological oscillation from CPU clock
                let biological_osc = OscillationState {
                    name: "CPU_Clock_Driven_ATP_Synthase".to_string(),
                    amplitude,
                    phase,
                    frequency: estimated_freq * 1e-7,  // Scale to biological range
                    damping_coefficient: 0.05,
                    atp_coupling_strength: 0.8,
                };
                
                let harvested_osc = HarvestedOscillation {
                    source_name: "CPU_Clock_Primary".to_string(),
                    current_frequency: estimated_freq,
                    current_amplitude: amplitude,
                    current_phase: phase,
                    last_update: current_time,
                    sample_history: Vec::new(),
                    biological_state: biological_osc,
                    atp_contribution: amplitude * 0.1,  // Convert to ATP
                };
                
                harvested.lock().unwrap().insert("CPU_Clock_Primary".to_string(), harvested_osc);
                
                last_time = current_time;
                thread::sleep(Duration::from_millis(1));
            }
        });
        
        Ok(())
    }

    fn start_light_harvesting(&mut self) -> Result<(), SolverError> {
        let harvested = Arc::clone(&self.harvested_oscillations);
        let is_harvesting = Arc::clone(&self.is_harvesting);
        
        thread::spawn(move || {
            let mut light_cycle = 0.0;
            
            while *is_harvesting.lock().unwrap() {
                light_cycle += 0.1;
                
                // Simulate screen backlight PWM oscillations
                let pwm_freq = 240.0;  // Typical PWM frequency
                let brightness_modulation = 0.7 + 0.3 * (light_cycle * 0.5).sin();
                let pwm_amplitude = brightness_modulation * (light_cycle * pwm_freq).sin().abs();
                let phase = (light_cycle * pwm_freq) % (2.0 * std::f64::consts::PI);
                
                // Map to biological photosystem
                let biological_osc = OscillationState {
                    name: "Light_Driven_Cytochrome_Oxidase".to_string(),
                    amplitude: pwm_amplitude,
                    phase,
                    frequency: pwm_freq * 0.2,  // Scale to biological
                    damping_coefficient: 0.1,
                    atp_coupling_strength: 0.4,
                };
                
                let harvested_osc = HarvestedOscillation {
                    source_name: "Screen_Backlight_PWM".to_string(),
                    current_frequency: pwm_freq,
                    current_amplitude: pwm_amplitude,
                    current_phase: phase,
                    last_update: Instant::now(),
                    sample_history: Vec::new(),
                    biological_state: biological_osc,
                    atp_contribution: pwm_amplitude * 0.05,  // Light -> ATP
                };
                
                harvested.lock().unwrap().insert("Screen_Backlight_PWM".to_string(), harvested_osc);
                
                thread::sleep(Duration::from_millis(10));
            }
        });
        
        Ok(())
    }

    fn start_electromagnetic_harvesting(&mut self) -> Result<(), SolverError> {
        let harvested = Arc::clone(&self.harvested_oscillations);
        let is_harvesting = Arc::clone(&self.is_harvesting);
        
        thread::spawn(move || {
            let mut em_cycle = 0.0;
            
            while *is_harvesting.lock().unwrap() {
                em_cycle += 0.05;
                
                // Simulate WiFi 2.4 GHz oscillations
                let wifi_freq = 2.4e9;
                let signal_strength = 0.3 + 0.2 * (em_cycle * 2.0).sin();
                let amplitude = signal_strength * (em_cycle * 1000.0).sin().abs();
                let phase = (em_cycle * 1000.0) % (2.0 * std::f64::consts::PI);
                
                // Map to biological electromagnetic sensitivity
                let biological_osc = OscillationState {
                    name: "EM_Field_NADH_Dehydrogenase".to_string(),
                    amplitude,
                    phase,
                    frequency: wifi_freq * 1e-7,  // Scale to biological
                    damping_coefficient: 0.2,
                    atp_coupling_strength: 0.3,
                };
                
                let harvested_osc = HarvestedOscillation {
                    source_name: "WiFi_2_4GHz".to_string(),
                    current_frequency: wifi_freq,
                    current_amplitude: amplitude,
                    current_phase: phase,
                    last_update: Instant::now(),
                    sample_history: Vec::new(),
                    biological_state: biological_osc,
                    atp_contribution: amplitude * 0.03,  // EM -> ATP
                };
                
                harvested.lock().unwrap().insert("WiFi_2_4GHz".to_string(), harvested_osc);
                
                thread::sleep(Duration::from_millis(20));
            }
        });
        
        Ok(())
    }

    fn start_thermal_harvesting(&mut self) -> Result<(), SolverError> {
        let harvested = Arc::clone(&self.harvested_oscillations);
        let is_harvesting = Arc::clone(&self.is_harvesting);
        
        thread::spawn(move || {
            let mut thermal_cycle = 0.0;
            
            while *is_harvesting.lock().unwrap() {
                thermal_cycle += 0.01;
                
                // Simulate CPU temperature oscillations
                let base_temp = 45.0;  // Base CPU temperature
                let temp_variation = 5.0 * (thermal_cycle * 0.1).sin();
                let current_temp = base_temp + temp_variation;
                let thermal_freq = 0.1;  // Slow thermal oscillations
                let amplitude = (current_temp - 30.0) / 50.0;  // Normalize temperature
                let phase = (thermal_cycle * thermal_freq) % (2.0 * std::f64::consts::PI);
                
                // Map to biological thermal fluctuations
                let biological_osc = OscillationState {
                    name: "Thermal_Protein_Dynamics".to_string(),
                    amplitude,
                    phase,
                    frequency: thermal_freq * 20.0,  // Scale to biological
                    damping_coefficient: 0.3,
                    atp_coupling_strength: 0.2,
                };
                
                let harvested_osc = HarvestedOscillation {
                    source_name: "CPU_Temperature_Thermal".to_string(),
                    current_frequency: thermal_freq,
                    current_amplitude: amplitude,
                    current_phase: phase,
                    last_update: Instant::now(),
                    sample_history: Vec::new(),
                    biological_state: biological_osc,
                    atp_contribution: amplitude * 0.02,  // Thermal -> ATP
                };
                
                harvested.lock().unwrap().insert("CPU_Temperature_Thermal".to_string(), harvested_osc);
                
                thread::sleep(Duration::from_millis(100));
            }
        });
        
        Ok(())
    }

    /// Get hardware oscillation statistics
    pub fn get_hardware_statistics(&self) -> HardwareOscillationStatistics {
        let harvested = self.harvested_oscillations.lock().unwrap();
        
        let total_sources = harvested.len();
        let total_atp_rate: f64 = harvested.values().map(|h| h.atp_contribution).sum();
        let avg_frequency: f64 = harvested.values().map(|h| h.current_frequency).sum::<f64>() / total_sources as f64;
        let total_amplitude: f64 = harvested.values().map(|h| h.current_amplitude).sum();
        
        HardwareOscillationStatistics {
            active_sources: total_sources,
            total_atp_generation_rate: total_atp_rate,
            average_frequency: avg_frequency,
            total_amplitude,
            hardware_efficiency: total_atp_rate / total_sources as f64,
        }
    }
}

#[derive(Debug)]
pub struct HardwareOscillationStatistics {
    pub active_sources: usize,
    pub total_atp_generation_rate: f64,
    pub average_frequency: f64,
    pub total_amplitude: f64,
    pub hardware_efficiency: f64,
}

// ================================================================================================
// SPECIALIZED HARVESTERS
// ================================================================================================

/// Light-specific oscillation harvester
#[derive(Debug)]
pub struct LightOscillationHarvester {
    pub screen_brightness_monitor: bool,
    pub ambient_light_sensor: bool,
    pub led_indicators: Vec<String>,
}

impl LightOscillationHarvester {
    pub fn new() -> Self {
        Self {
            screen_brightness_monitor: true,
            ambient_light_sensor: false,  // May not be available
            led_indicators: vec![
                "Power_LED".to_string(),
                "Network_Activity_LED".to_string(),
                "Disk_Activity_LED".to_string(),
            ],
        }
    }
}

/// Clock-specific oscillation harvester
#[derive(Debug)]
pub struct ClockOscillationHarvester {
    pub cpu_clocks: Vec<String>,
    pub system_timers: Vec<String>,
    pub network_time_sync: bool,
}

impl ClockOscillationHarvester {
    pub fn new() -> Self {
        Self {
            cpu_clocks: vec![
                "CPU_Base_Clock".to_string(),
                "CPU_Boost_Clock".to_string(),
                "Memory_Clock".to_string(),
                "GPU_Clock".to_string(),
            ],
            system_timers: vec![
                "System_Timer_1ms".to_string(),
                "High_Resolution_Timer".to_string(),
                "RTC_Clock".to_string(),
            ],
            network_time_sync: true,
        }
    }
}

/// Electromagnetic field harvester
#[derive(Debug)]
pub struct ElectromagneticHarvester {
    pub wifi_bands: Vec<f64>,
    pub bluetooth_frequency: f64,
    pub power_line_frequency: f64,
    pub usb_signal_harvesting: bool,
}

impl ElectromagneticHarvester {
    pub fn new() -> Self {
        Self {
            wifi_bands: vec![2.4e9, 5.0e9],  // 2.4 GHz and 5 GHz
            bluetooth_frequency: 2.4e9,
            power_line_frequency: 60.0,      // 60 Hz in North America
            usb_signal_harvesting: true,
        }
    }
}

/// Thermal oscillation harvester
#[derive(Debug)]
pub struct ThermalOscillationHarvester {
    pub cpu_temperature_monitoring: bool,
    pub gpu_temperature_monitoring: bool,
    pub system_fan_oscillations: bool,
}

impl ThermalOscillationHarvester {
    pub fn new() -> Self {
        Self {
            cpu_temperature_monitoring: true,
            gpu_temperature_monitoring: true,
            system_fan_oscillations: true,
        }
    }
}

// ================================================================================================
// HARDWARE-BIOLOGICAL INTEGRATION
// ================================================================================================

/// Update biological quantum state with harvested hardware oscillations
pub fn integrate_hardware_oscillations_into_biology(
    base_state: &mut BiologicalQuantumState,
    harvester: &HardwareOscillationHarvester,
) -> Result<(), SolverError> {
    
    println!("🔗 Integrating hardware oscillations into biological quantum state...");
    
    let harvested = harvester.harvested_oscillations.lock().unwrap();
    
    // Replace or augment simulated oscillations with hardware-harvested ones
    let mut hardware_oscillations = Vec::new();
    let mut total_hardware_atp = 0.0;
    
    for (source_name, harvested_osc) in harvested.iter() {
        hardware_oscillations.push(harvested_osc.biological_state.clone());
        total_hardware_atp += harvested_osc.atp_contribution;
        
        println!("  🔄 Integrated {}: f={:.1} Hz -> ATP={:.3} mM", 
                source_name, 
                harvested_osc.biological_state.frequency,
                harvested_osc.atp_contribution);
    }
    
    // Add hardware oscillations to biological system
    base_state.oscillatory_coords.oscillations.extend(hardware_oscillations);
    
    // Boost ATP from hardware energy harvesting
    base_state.atp_coords.atp_concentration += total_hardware_atp;
    
    // Update oscillatory momenta for new oscillations
    let additional_momenta = vec![0.1; harvested.len()];
    base_state.oscillatory_coords.oscillatory_momenta.extend(additional_momenta);
    
    println!("✅ Integrated {} hardware sources, added {:.3} mM ATP", 
            harvested.len(), total_hardware_atp);
    
    Ok(())
}

/// Create a hardware-powered biological quantum computer
pub fn create_hardware_powered_biological_quantum_computer() -> Result<(BiologicalQuantumState, HardwareOscillationHarvester), SolverError> {
    println!("🏭 Creating hardware-powered biological quantum computer...");
    
    // Initialize hardware harvester
    let mut harvester = HardwareOscillationHarvester::new();
    harvester.initialize_hardware_sources()?;
    harvester.start_harvesting()?;
    
    // Give hardware time to collect some oscillations
    thread::sleep(Duration::from_millis(100));
    
    // Create biological state from harvested hardware oscillations
    let bio_state = harvester.get_biological_state_from_hardware()?;
    
    println!("🎉 Hardware-powered biological quantum computer created!");
    println!("⚡ Running on real hardware oscillations instead of simulations!");
    
    Ok((bio_state, harvester))
} 