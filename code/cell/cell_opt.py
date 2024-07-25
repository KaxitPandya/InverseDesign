# -*- coding: utf-8 -*-
# -------------------------------------------------------------
# Bulk Configuration
# -------------------------------------------------------------

# Set up lattice
vector_a = [8.651274691122111, 0.0, 0.0]*Angstrom
vector_b = [0.0, 5.767516460748075, 0.0]*Angstrom
vector_c = [0.0, 0.0, 16.37325]*Angstrom
lattice = UnitCell(vector_a, vector_b, vector_c)

# Define elements
elements = [Gold, Gold, Gold, Gold, Gold, Gold, Gold, Gold, Gold, Gold, Gold,
            Gold, Sulfur, Sulfur, Sulfur, Sulfur, Sulfur, Sulfur, Molybdenum,
            Molybdenum, Molybdenum, Molybdenum, Molybdenum, Molybdenum, Sulfur,
            Sulfur, Sulfur, Sulfur, Sulfur, Sulfur]

# Define coordinates
fractional_coordinates = [[ 0.166666666667,  0.25          ,  0.124540027178],
                          [ 0.5           ,  0.25          ,  0.124540027178],
                          [ 0.833333333333,  0.25          ,  0.124540027178],
                          [ 0.166666666667,  0.75          ,  0.124540027178],
                          [ 0.5           ,  0.75          ,  0.124540027178],
                          [ 0.833333333333,  0.75          ,  0.124540027178],
                          [ 0.            ,  0.            ,  0.249080054357],
                          [ 0.333333333333,  0.            ,  0.249080054357],
                          [ 0.666666666667,  0.            ,  0.249080054357],
                          [ 0.            ,  0.5           ,  0.249080054357],
                          [ 0.333333333333,  0.5           ,  0.249080054357],
                          [ 0.666666666667,  0.5           ,  0.249080054357],
                          [ 0.            ,  0.333333333333,  0.33994136778 ],
                          [ 0.333333333333,  0.333333333333,  0.33994136778 ],
                          [ 0.666666666667,  0.333333333333,  0.33994136778 ],
                          [ 0.166666666667,  0.833333333333,  0.33994136778 ],
                          [ 0.5           ,  0.833333333333,  0.33994136778 ],
                          [ 0.833333333333,  0.833333333333,  0.33994136778 ],
                          [ 0.166666666667,  0.166666666667,  0.436810040768],
                          [ 0.5           ,  0.166666666667,  0.436810040768],
                          [ 0.833333333333,  0.166666666667,  0.436810040768],
                          [ 0.            ,  0.666666666667,  0.436810040768],
                          [ 0.333333333333,  0.666666666667,  0.436810040768],
                          [ 0.666666666667,  0.666666666667,  0.436810040768],
                          [ 0.            ,  0.333333333333,  0.533678713756],
                          [ 0.333333333333,  0.333333333333,  0.533678713756],
                          [ 0.666666666667,  0.333333333333,  0.533678713756],
                          [ 0.166666666667,  0.833333333333,  0.533678713756],
                          [ 0.5           ,  0.833333333333,  0.533678713756],
                          [ 0.833333333333,  0.833333333333,  0.533678713756]]

# Set up configuration
bulk_configuration = BulkConfiguration(
    bravais_lattice=lattice,
    elements=elements,
    fractional_coordinates=fractional_coordinates
    )

# Add tags
bulk_configuration.addTags('Left Interface 0',  [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11])
bulk_configuration.addTags('Right Interface 0', [12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,
                                                 25, 26, 27, 28, 29])

# -------------------------------------------------------------
# Calculator
# -------------------------------------------------------------
#----------------------------------------
# Basis Set
#----------------------------------------

# Basis set for Sulfur

sulfur_3p = ConfinedOrbital(
    principal_quantum_number=3,
    angular_momentum=1,
    radial_cutoff_radius=5.58618376198*Bohr,
    confinement_start_radius=4.46894700958*Bohr,
    additional_charge=0,
    confinement_strength=17.9013087039*Hartree,
    confinement_power=1,
    radial_step_size=0.001*Bohr,
    )

sulfur_3p = ConfinedOrbital(
    principal_quantum_number=3,
    angular_momentum=1,
    radial_cutoff_radius=5.58618376198*Bohr,
    confinement_start_radius=4.46894700958*Bohr,
    additional_charge=0,
    confinement_strength=17.9013087039*Hartree,
    confinement_power=1,
    radial_step_size=0.001*Bohr,
    )

sulfur_3p_polarization = PolarizationOrbital(sulfur_3p)

sulfur_3p = ConfinedOrbital(
    principal_quantum_number=3,
    angular_momentum=1,
    radial_cutoff_radius=5.58618376198*Bohr,
    confinement_start_radius=4.46894700958*Bohr,
    additional_charge=0,
    confinement_strength=17.9013087039*Hartree,
    confinement_power=1,
    radial_step_size=0.001*Bohr,
    )

sulfur_3p_split = AnalyticalSplit(sulfur_3p, split_norm=0.15)

sulfur_3s = ConfinedOrbital(
    principal_quantum_number=3,
    angular_momentum=0,
    radial_cutoff_radius=4.51478507146*Bohr,
    confinement_start_radius=3.61182805717*Bohr,
    additional_charge=0,
    confinement_strength=22.1494486265*Hartree,
    confinement_power=1,
    radial_step_size=0.001*Bohr,
    )

sulfur_3s = ConfinedOrbital(
    principal_quantum_number=3,
    angular_momentum=0,
    radial_cutoff_radius=4.51478507146*Bohr,
    confinement_start_radius=3.61182805717*Bohr,
    additional_charge=0,
    confinement_strength=22.1494486265*Hartree,
    confinement_power=1,
    radial_step_size=0.001*Bohr,
    )

sulfur_3s_split = AnalyticalSplit(sulfur_3s, split_norm=0.15)

SulfurBasis = BasisSet(
    element=PeriodicTable.Sulfur,
    orbitals=[sulfur_3s, sulfur_3p, sulfur_3s_split, sulfur_3p_split, sulfur_3p_polarization],
    occupations=[2.0, 4.0, 0.0, 0.0, 0.0],
    hubbard_u=[0.0, 0.0, 0.0, 0.0, 0.0]*eV,
    dft_half_parameters=Automatic,
    filling_method=SphericalSymmetric,
    onsite_spin_orbit_split=[0.0, 0.0, 0.0, 0.0, 0.0]*eV,
    pseudopotential=NormConservingPseudoPotential("normconserving/S.GGAPBE.zip"),
    )


# Basis set for Molybdenum

molybdenum_4d = ConfinedOrbital(
    principal_quantum_number=4,
    angular_momentum=2,
    radial_cutoff_radius=5.94363171581*Bohr,
    confinement_start_radius=4.75490537265*Bohr,
    additional_charge=0,
    confinement_strength=16.8247301955*Hartree,
    confinement_power=1,
    radial_step_size=0.001*Bohr,
    )

molybdenum_4d = ConfinedOrbital(
    principal_quantum_number=4,
    angular_momentum=2,
    radial_cutoff_radius=5.94363171581*Bohr,
    confinement_start_radius=4.75490537265*Bohr,
    additional_charge=0,
    confinement_strength=16.8247301955*Hartree,
    confinement_power=1,
    radial_step_size=0.001*Bohr,
    )

molybdenum_4d_split = AnalyticalSplit(molybdenum_4d, split_norm=0.15)

molybdenum_5s = ConfinedOrbital(
    principal_quantum_number=5,
    angular_momentum=0,
    radial_cutoff_radius=7.94602904419*Bohr,
    confinement_start_radius=6.35682323535*Bohr,
    additional_charge=0,
    confinement_strength=12.5849024014*Hartree,
    confinement_power=1,
    radial_step_size=0.001*Bohr,
    )

molybdenum_5s = ConfinedOrbital(
    principal_quantum_number=5,
    angular_momentum=0,
    radial_cutoff_radius=7.94602904419*Bohr,
    confinement_start_radius=6.35682323535*Bohr,
    additional_charge=0,
    confinement_strength=12.5849024014*Hartree,
    confinement_power=1,
    radial_step_size=0.001*Bohr,
    )

molybdenum_5s_polarization = PolarizationOrbital(molybdenum_5s)

molybdenum_5s = ConfinedOrbital(
    principal_quantum_number=5,
    angular_momentum=0,
    radial_cutoff_radius=7.94602904419*Bohr,
    confinement_start_radius=6.35682323535*Bohr,
    additional_charge=0,
    confinement_strength=12.5849024014*Hartree,
    confinement_power=1,
    radial_step_size=0.001*Bohr,
    )

molybdenum_5s_split = AnalyticalSplit(molybdenum_5s, split_norm=0.15)

MolybdenumBasis = BasisSet(
    element=PeriodicTable.Molybdenum,
    orbitals=[molybdenum_4d, molybdenum_5s, molybdenum_4d_split, molybdenum_5s_split, molybdenum_5s_polarization],
    occupations=[5.0, 1.0, 0.0, 0.0, 0.0],
    hubbard_u=[0.0, 0.0, 0.0, 0.0, 0.0]*eV,
    dft_half_parameters=Automatic,
    filling_method=SphericalSymmetric,
    onsite_spin_orbit_split=[0.0, 0.0, 0.0, 0.0, 0.0]*eV,
    pseudopotential=NormConservingPseudoPotential("normconserving/MO.GGAPBE.zip"),
    )


# Basis set for Gold

gold_5d = ConfinedOrbital(
    principal_quantum_number=5,
    angular_momentum=2,
    radial_cutoff_radius=4.97009209077*Bohr,
    confinement_start_radius=3.97607367262*Bohr,
    additional_charge=0,
    confinement_strength=20.1203515294*Hartree,
    confinement_power=1,
    radial_step_size=0.001*Bohr,
    )

gold_5d = ConfinedOrbital(
    principal_quantum_number=5,
    angular_momentum=2,
    radial_cutoff_radius=4.97009209077*Bohr,
    confinement_start_radius=3.97607367262*Bohr,
    additional_charge=0,
    confinement_strength=20.1203515294*Hartree,
    confinement_power=1,
    radial_step_size=0.001*Bohr,
    )

gold_5d_split = AnalyticalSplit(gold_5d, split_norm=0.15)

gold_6s = ConfinedOrbital(
    principal_quantum_number=6,
    angular_momentum=0,
    radial_cutoff_radius=6.80550732918*Bohr,
    confinement_start_radius=5.44440586334*Bohr,
    additional_charge=0,
    confinement_strength=14.6939816773*Hartree,
    confinement_power=1,
    radial_step_size=0.001*Bohr,
    )

gold_6s = ConfinedOrbital(
    principal_quantum_number=6,
    angular_momentum=0,
    radial_cutoff_radius=6.80550732918*Bohr,
    confinement_start_radius=5.44440586334*Bohr,
    additional_charge=0,
    confinement_strength=14.6939816773*Hartree,
    confinement_power=1,
    radial_step_size=0.001*Bohr,
    )

gold_6s_polarization = PolarizationOrbital(gold_6s)

gold_6s = ConfinedOrbital(
    principal_quantum_number=6,
    angular_momentum=0,
    radial_cutoff_radius=6.80550732918*Bohr,
    confinement_start_radius=5.44440586334*Bohr,
    additional_charge=0,
    confinement_strength=14.6939816773*Hartree,
    confinement_power=1,
    radial_step_size=0.001*Bohr,
    )

gold_6s_split = AnalyticalSplit(gold_6s, split_norm=0.15)

GoldBasis = BasisSet(
    element=PeriodicTable.Gold,
    orbitals=[gold_5d, gold_6s, gold_5d_split, gold_6s_split, gold_6s_polarization],
    occupations=[10.0, 1.0, 0.0, 0.0, 0.0],
    hubbard_u=[0.0, 0.0, 0.0, 0.0, 0.0]*eV,
    dft_half_parameters=Automatic,
    filling_method=SphericalSymmetric,
    onsite_spin_orbit_split=[0.0, 0.0, 0.0, 0.0, 0.0]*eV,
    pseudopotential=NormConservingPseudoPotential("normconserving/AU.GGAPBE.zip"),
    )

basis_set = [
    SulfurBasis,
    MolybdenumBasis,
    GoldBasis,
    ]

#----------------------------------------
# Exchange-Correlation
#----------------------------------------
exchange_correlation = GGA.PBE

k_point_sampling = MonkhorstPackGrid(
    na=3,
    nb=5,
    nc=1,
    symmetries=[
        ([[ 1., 0., 0.],
          [ 0., 1., 0.],
          [ 0., 0., 1.]], [ 0., 0., 0.]),
        ([[-1., 0., 0.],
          [ 0.,-1., 0.],
          [ 0., 0.,-1.]], [ 0., 0., 0.]),
        ],
    force_timereversal=True,
    shift_to_gamma=[True, True, True],
    )
exact_exchange_parameters = ExactExchangeParameters(
    aux_basis_tolerance=0.001,
    number_of_waves=1024,
    maximum_wave_vector=50.0,
    integral_tolerance=0.0001,
    relative_screening_tolerance=0.01,
    )
numerical_accuracy_parameters = NumericalAccuracyParameters(
    density_mesh_cutoff=60.0*Hartree,
    k_point_sampling=k_point_sampling,
    radial_step_size=0.001*Bohr,
    density_cutoff=1e-06,
    interaction_max_range=20.0*Angstrom,
    number_of_reciprocal_points=1024,
    reciprocal_energy_cutoff=1250.0*Hartree,
    bands_per_electron=1.2,
    occupation_method=FermiDirac(1000.0*Kelvin*boltzmann_constant),
    exact_exchange_parameters=exact_exchange_parameters,
    compensation_charge_mesh_cutoff=75.0*Hartree,
    )

iteration_control_parameters = IterationControlParameters(
    tolerance=0.0001,
    max_steps=100,
    algorithm=PulayMixer(),
    damping_factor=0.1,
    number_of_history_steps=20,
    start_mixing_after_step=0,
    mixing_variable=HamiltonianVariable,
    preconditioner=Preconditioner.Off,
    linear_dependence_threshold=0.0,
    max_exx_updates=50,
    non_convergence_behavior=ContinueCalculation(),
    enable_scf_stop_file=True,
    )

poisson_solver = FastFourierSolver()

density_matrix_method = DiagonalizationSolver(
    bands_above_fermi_level=Automatic,
    processes_per_kpoint=None,
    )
algorithm_parameters = AlgorithmParameters(
    density_matrix_method=density_matrix_method,
    store_grids=True,
    store_basis_on_grid=Automatic,
    store_energy_density_matrix=Automatic,
    scf_restart_step_length=0.1*Angstrom,
    use_symmetries=True,
    )

parallel_parameters = ParallelParameters(
    processes_per_neb_image=None,
    processes_per_individual=None,
    processes_per_bias_point=None,
    processes_per_saddle_search=1,
    )

#----------------------------------------
# Grimme DFTD3
#----------------------------------------
correction_extension = GrimmeDFTD3(
    exchange_correlation=GGA.PBE,
    global_scale_factor_6=1.0,
    damping_factor_6=1.217,
    global_scale_factor_8=0.722,
    damping_factor_8=1.0,
    maximum_neighbour_distance=30.0*Ang,
    include_three_body_term=False,
    )

checkpoint_handler = CheckpointHandler(
    time_interval=0.5*Hour,
    )

calculator = LCAOCalculator(
    basis_set=basis_set,
    exchange_correlation=exchange_correlation,
    numerical_accuracy_parameters=numerical_accuracy_parameters,
    iteration_control_parameters=iteration_control_parameters,
    poisson_solver=poisson_solver,
    algorithm_parameters=algorithm_parameters,
    parallel_parameters=parallel_parameters,
    correction_extension=correction_extension,
    checkpoint_handler=checkpoint_handler,
    charge=0.0,
    fixed_spin_moment=False,
    )

bulk_configuration.setCalculator(calculator)
nlprint(bulk_configuration)
bulk_configuration.update()
nlsave('cell.hdf5', bulk_configuration)

# -------------------------------------------------------------
# Optimize Geometry
# -------------------------------------------------------------
bulk_configuration = OptimizeGeometry(
    bulk_configuration,
    max_forces=0.05*eV/Ang,
    max_steps=200,
    max_step_length=0.2*Ang,
    trajectory_filename='cell_trajectory.hdf5',
    trajectory_interval=5.0*Minute,
    restart_strategy=RestartFromTrajectory(),
    disable_stress=True,
    optimizer_method=LBFGS(),
    enable_optimization_stop_file=True,
)
nlsave('cell.hdf5', bulk_configuration)
nlprint(bulk_configuration)
