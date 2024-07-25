# -*- coding: utf-8 -*-
# -------------------------------------------------------------
# Bulk Configuration
# -------------------------------------------------------------

# Set up lattice
vector_a = [3.1604, 0.0, 0.0]*Angstrom
vector_b = [0.0, 5.47397, 0.0]*Angstrom
vector_c = [0.0, 0.0, 25.0]*Angstrom
lattice = UnitCell(vector_a, vector_b, vector_c)

# Define elements
elements = [Molybdenum, Molybdenum, Sulfur, Sulfur, Sulfur, Sulfur]

# Define coordinates
fractional_coordinates = [[ 0.25          ,  0.416666615329,  0.5           ],
                          [ 0.75          ,  0.916666923354,  0.5           ],
                          [ 0.25          ,  0.083333076646,  0.4365578     ],
                          [ 0.75          ,  0.583333384671,  0.4365578     ],
                          [ 0.25          ,  0.083333076646,  0.5634422     ],
                          [ 0.75          ,  0.583333384671,  0.5634422     ]]

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
    GGABasis.Sulfur_DoubleZetaPolarized,
    GGABasis.Molybdenum_DoubleZetaPolarized,
    ]

k_point_sampling = KpointDensity(
    density_a=4.0*Angstrom,
    density_c=0.0*Angstrom,
    )
numerical_accuracy_parameters = NumericalAccuracyParameters(
    k_point_sampling=k_point_sampling,
    )

iteration_control_parameters = IterationControlParameters(
    damping_factor=0.2,
    number_of_history_steps=5,
    preconditioner=Kerker(energy_q0=0.01*Hartree, energy_qmax=80.0*Hartree, maximum_damping=0.01),
    )

#----------------------------------------
# Grimme DFTD3
#----------------------------------------
correction_extension = GrimmeDFTD3(
    exchange_correlation=GGA.PBE,
    maximum_neighbour_distance=30.0*Ang,
    include_three_body_term=False,
    )

calculator = PlaneWaveCalculator(
    wave_function_cutoff=20.0*Hartree,
    basis_set=basis_set,
    numerical_accuracy_parameters=numerical_accuracy_parameters,
    iteration_control_parameters=iteration_control_parameters,
    correction_extension=correction_extension,
    )

bulk_configuration.setCalculator(calculator)
nlprint(bulk_configuration)
bulk_configuration.update()
nlsave('Molybdenite.hdf5', bulk_configuration)

# -------------------------------------------------------------
# Optimize Geometry
# -------------------------------------------------------------

constraints = [FixStrain(x=False, y=False, z=True)]

bulk_configuration = OptimizeGeometry(
    bulk_configuration,
    max_forces=0.05*eV/Ang,
    max_stress=0.1*GPa,
    max_steps=200,
    max_step_length=0.2*Ang,
    constraints=constraints,
    trajectory_filename='Molybdenite_trajectory.hdf5',
    trajectory_interval=5.0*Minute,
    restart_strategy=RestartFromTrajectory(),
    optimizer_method=LBFGS(),
    enable_optimization_stop_file=True,
)
nlsave('Molybdenite.hdf5', bulk_configuration)
nlprint(bulk_configuration)
