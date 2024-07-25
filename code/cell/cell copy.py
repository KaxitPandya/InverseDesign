# Set up lattice
vector_a = [2.8837582303740374, 0.0, 0.0]*Angstrom
vector_b = [0.0, 2.8837582303740374, 0.0]*Angstrom
vector_c = [0.0, 0.0, 24.469499999999996]*Angstrom
lattice = UnitCell(vector_a, vector_b, vector_c)

# Define elements
elements = [Gold, Gold, Gold, Gold, Gold, Gold, Gold, Gold, Gold, Gold, Gold,
            Gold]

# Define coordinates
fractional_coordinates = [[ 0.5           ,  0.5           ,  0.083333333333],
                          [ 0.5           ,  0.5           ,  0.25          ],
                          [ 0.5           ,  0.5           ,  0.416666666667],
                          [ 0.5           ,  0.5           ,  0.583333333333],
                          [ 0.5           ,  0.5           ,  0.75          ],
                          [ 0.5           ,  0.5           ,  0.916666666667],
                          [ 0.            ,  0.            ,  0.166666666667],
                          [ 0.            ,  0.            ,  0.333333333333],
                          [ 0.            ,  0.            ,  0.5           ],
                          [ 0.            ,  0.            ,  0.666666666667],
                          [ 0.            ,  0.            ,  0.833333333333],
                          [ 0.            ,  0.            ,  1.            ]]

# Set up configuration
bulk_configuration_gold = BulkConfiguration(
    bravais_lattice=lattice,
    elements=elements,
    fractional_coordinates=fractional_coordinates
    )

vector_a = [3.1604, 0.0, 0.0]*Angstrom
vector_b = [0.0, 5.47397337224068, 0.0]*Angstrom
vector_c = [0.0, 0.0, 12.295]*Angstrom
lattice = UnitCell(vector_a, vector_b, vector_c)

# Define elements
elements = [Molybdenum, Sulfur, Sulfur, Molybdenum, Sulfur, Sulfur]

# Define coordinates
fractional_coordinates = [[-0.            ,  0.666666666667,  0.25          ],
                          [-0.            ,  0.333333333333,  0.121         ],
                          [-0.            ,  0.333333333333,  0.379         ],
                          [ 0.5           ,  0.166666666667,  0.25          ],
                          [ 0.5           ,  0.833333333333,  0.121         ],
                          [ 0.5           ,  0.833333333333,  0.379         ]]

# Set up configuration
bulk_configuration_mo = BulkConfiguration(
    bravais_lattice=lattice,
    elements=elements,
    fractional_coordinates=fractional_coordinates
    )

# Setup interface builder
interface_builder = InterfaceBuilder(bulk_configuration_gold, bulk_configuration_mo)

# Get all of the interface matches, specifying that the second surface (MoS2) should be strained.
matches = interface_builder.matches(strain_method=StrainBoth)

# Print out table of matches
interface_builder.nlprint()

# print('{:>20s} {:>10s}'.format('number of atoms', 'strain (%)'))
for match in matches:
    # Get the number of atoms and strain from the match.
    number_of_atoms = match.numberOfAtoms()
    mean_absolute_strain = 100.0 * match.meanAbsoluteStrain()

    # Print out information about the match.
    print('{:20d} {:10.3f}'.format(number_of_atoms, mean_absolute_strain))

    # Make the interface configuration and save it to a file.
    interface_configuration = match.makeInterface()
    nlsave('matches.hdf5', interface_configuration)

# Scatter plot of number of atoms vs strain.
import pylab
pylab.scatter(
    [100 * match.meanAbsoluteStrain() for match in matches],
    [match.numberOfAtoms() for match in matches],
)
pylab.xlabel('Mean Absolute Strain (%)')
pylab.ylabel('Number of Atoms')
pylab.xscale('log')
pylab.savefig('size_strain_plot.png')