export System, Particle, Velocity, Force, sys_init, System_Type, No_Mol, Molecular

# Defining data structures.
abstract type System_Type end

abstract type Thermostat end 

abstract type Bonded_Interactions end
abstract type Non_Bonded_Interactions end
    
abstract type TrajectoryDump end
abstract type Logger end

abstract type Shear end

struct No_Mol <: System_Type end
struct Molecular <: System_Type end

struct No_Thermostat <: Thermostat end
    
struct No_Bonds <: Bonded_Interactions end
struct No_Non_Bonded <: Non_Bonded_Interactions end
    
struct No_Logger <: Logger end
struct No_TrajectoryDump <: TrajectoryDump end

struct No_Shear <: Shear end

mutable struct System
    system_type::System_Type

    # Special features:
    shear::Shear
    shear_freq::Int64
    j::Float64

    bonded_interactions::Bonded_Interactions
    non_bonded_interactions::Non_Bonded_Interactions
    thermostat::Thermostat

    init_velocity::AbstractString
    ∆t::Float64
    T::Float64
    box::Vector{Float64}
    
    n_steps::Int64
    traj_freq::Int64
    log_freq::Int64

    data_file::AbstractString
    traj_file::AbstractString
    log_file::AbstractString
    shear_file::AbstractString
    momentum_file::AbstractString
    restart_file::AbstractString
    
    trajectorydump::TrajectoryDump
    logger::Logger
    restart::Bool

    first_step::Bool
end

mutable struct Particle
    i_mol::Int64
    spec::Int64
    mass::Float64
    charge::Float64
    pos::Vector{Float64}
    image::Vector{Int64}
end

mutable struct Velocity
    velocity::Vector{Float64}
end
        
mutable struct Force
    force::Vector{Float64}
end
        
struct Bond
    type::Int64
    p_1::Int64
    p_2::Int64
end

struct Angle
    type::Int64
    p_1::Int64
    p_2::Int64
    p_3::Int64
end

# Initialize system. 
function sys_init()
    system_type = No_Mol()

    shear = No_Shear()
    shear_freq = 0
    j = 0

    bonded_interactions = No_Bonds()
    non_bonded_interactions = No_Non_Bonded()
    thermostat = No_Thermostat()

    init_velocity = "Z"
    ∆t = 0.0
    T = 0.0
    box = [0.0, 0.0, 0.0]

    n_steps = 0
    traj_freq = 0
    log_freq = 0 
    
    data_file = "input.data"
    traj_file = "trj.xyz"
    log_file = "logger_MD.txt"
    shear_file = "velocity_profile.txt"
    momentum_file = "momentum.txt"
    restart_file = "restart.data"
    
    trajectorydump = No_TrajectoryDump()
    logger = No_Logger()
    restart = false

    first_step = true
    return System(system_type, shear, shear_freq, j, bonded_interactions, non_bonded_interactions, thermostat, init_velocity, ∆t, T, box, n_steps, traj_freq, log_freq, data_file, traj_file, log_file, shear_file, momentum_file, restart_file, trajectorydump, logger, restart, first_step)
end