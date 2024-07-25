# -*- coding: utf-8 -*-
# -------------------------------------------------------------
# Bulk Configuration
# -------------------------------------------------------------

# Set up lattice
vector_a = [-4.059121037287073e-19, 4.163069473384337, -2.5984133191693297e-19]*Angstrom
vector_b = [-8.082671082026988e-19, -2.5984133191693297e-19, 4.16306947338436]*Angstrom
vector_c = [4.1630694733843665, -4.059121037287073e-19, -8.082671082026988e-19]*Angstrom
lattice = UnitCell(vector_a, vector_b, vector_c)

# Define elements
elements = [Gold, Gold, Gold, Gold]

# Define coordinates
fractional_coordinates = [[ 0.            ,  0.000000000001,  0.            ],
                          [ 0.5           ,  0.000000000001,  0.5           ],
                          [ 0.            ,  0.499999999999,  0.5           ],
                          [ 0.5           ,  0.499999999999,  0.            ]]

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
    density_mesh_cutoff=30.0*Hartree,
    k_point_sampling=k_point_sampling,
    )

#----------------------------------------
# Grimme DFTD3
#----------------------------------------
correction_extension = GrimmeDFTD3(
    exchange_correlation=GGA.PBE,
    maximum_neighbour_distance=30.0*Ang,
    include_three_body_term=False,
    )

calculator = LCAOCalculator(
    basis_set=basis_set,
    numerical_accuracy_parameters=numerical_accuracy_parameters,
    correction_extension=correction_extension,
    )

bulk_configuration.setCalculator(calculator)
nlprint(bulk_configuration)
bulk_configuration.update()
nlsave('/home/adminpc/Tarun/ATK/Train-data-Kaxit/Au-MoS2-Au/golden/Relax/Au/GoldLCAO.hdf5', bulk_configuration)

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
    trajectory_filename='/home/adminpc/Tarun/ATK/Train-data-Kaxit/Au-MoS2-Au/golden/Relax/Au/GoldLCAO_trajectory.hdf5',
    trajectory_interval=5.0*Minute,
    restart_strategy=RestartFromTrajectory(),
    optimizer_method=LBFGS(),
    enable_optimization_stop_file=True,
)
nlsave('/home/adminpc/Tarun/ATK/Train-data-Kaxit/Au-MoS2-Au/golden/Relax/Au/GoldLCAO.hdf5', bulk_configuration)
nlprint(bulk_configuration)
