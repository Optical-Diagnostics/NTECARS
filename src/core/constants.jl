

const μ_0::Float64 = ustrip(VacuumMagneticPermeability)
const ε_0::Float64 = ustrip(VacuumElectricPermittivity)
const ħ::Float64   = ustrip(ReducedPlanckConstant)
const c::Float64   = ustrip(SpeedOfLightInVacuum)
const k_B::Float64 = ustrip(BoltzmannConstant)
const c₂::Float64  = 1.438770
const h::Float64   = ustrip(PlanckConstant)

const cm⁻¹_to_J::Float64 = 100.0*h*c
const cm⁻¹_to_angular_frequency::Float64 = 2.0*π*c*100.0
const atm_to_Pa::Float64 = 101325.0
const cm⁻¹_to_m⁻¹::Float64 = 100.0