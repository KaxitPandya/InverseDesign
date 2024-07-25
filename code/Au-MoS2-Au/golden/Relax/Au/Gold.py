# -*- coding: utf-8 -*-
# -------------------------------------------------------------
# Bulk Configuration
# -------------------------------------------------------------

# Set up lattice
vector_a = [0.0, 4.07825, 0.0]*Angstrom
vector_b = [0.0, 0.0, 4.07825]*Angstrom
vector_c = [4.07825, 0.0, 0.0]*Angstrom
lattice = UnitCell(vector_a, vector_b, vector_c)

# Define elements
elements = [Gold, Gold, Gold, Gold]

# Define coordinates
fractional_coordinates = [[ 0. ,  0. ,  0. ],
                          [ 0.5,  0. ,  0.5],
                          [ 0. ,  0.5,  0.5],
                          [ 0.5,  0.5,  0. ]]

# Set up configuration
bulk_configuration = BulkConfiguration(
    bravais_lattice=lattice,
    elements=elements,
    fractional_coordinates=fractional_coordinates
    )

# -------------------------------------------------------------
# Calculator
# -------------------------------------------------------------
#----------------------------------------
# Basis Set
#----------------------------------------
basis_set = [
    GGABasis.Gold_DoubleZetaPolarized,
    ]

k_point_sampling = KpointDensity(
    density_a=4.0*Angstrom,
    )
numerical_accuracy_parameters = NumericalAccuracyParameters(
    k_point_sampling=k_point_sampling,
    )

iteration_control_parameters = IterationControlParameters(
    damping_factor=0.2,
    number_of_history_steps=5,
    preconditioner=Kerker(energy_q0=0.01*Hartree, energy_qmax=80.0*Hartree, maximum_damping=0.01),
    )

calculator = PlaneWaveCalculator(
    wave_function_cutoff=20.0*Hartree,
    basis_set=basis_set,
    numerical_accuracy_parameters=numerical_accuracy_parameters,
    iteration_control_parameters=iteration_control_parameters,
    )

bulk_configuration.setCalculator(calculator)
nlprint(bulk_configuration)
bulk_configuration.update()
nlsave('/home/adminpc/Tarun/ATK/Train-data-Kaxit/Au-MoS2-Au/golden/Relax/Au/Gold.hdf5', bulk_configuration)

# -------------------------------------------------------------
# Optimize Geometry
# -------------------------------------------------------------

constraints = [BravaisLatticeConstraint()]

bulk_configuration = OptimizeGeometry(
    bulk_configuration,
    max_forces=0.05*eV/Ang,
    max_stress=0.1*GPa,
    max_steps=200,
    max_step_length=0.2*Ang,
    constraints=constraints,
    trajectory_filename='/home/adminpc/Tarun/ATK/Train-data-Kaxit/Au-MoS2-Au/golden/Relax/Au/Gold_trajectory.hdf5',
    trajectory_interval=5.0*Minute,
    restart_strategy=RestartFromTrajectory(),
    optimizer_method=LBFGS(),
    enable_optimization_stop_file=True,
)
nlsave('/home/adminpc/Tarun/ATK/Train-data-Kaxit/Au-MoS2-Au/golden/Relax/Au/Gold.hdf5', bulk_configuration)
nlprint(bulk_configuration)
