vector_a = [14.418791151870188, 0.0, 0.0]*Angstrom
vector_b = [0.0, -8.651274691122111, 0.0]*Angstrom
vector_c = [0.0, 0.0, 36.7645]*Angstrom
lattice = UnitCell(vector_a, vector_b, vector_c)

elements_gold = [
Gold,
Gold,
Gold,
Gold,
Gold,
Gold,
Gold,
Gold,
Gold,
Gold,
Gold,
Gold,
]

fractional_coordinates_gold = [
    [0.5, 0.5, 0.07692307692307693],
    [0.0, 0.0, 0.15384615384615385],
    [0.5, 0.5, 0.23076923076923078],
    [0.0, 0.0, 0.3076923076923077],
    [0.5, 0.5, 0.38461538461538464],
    [0.0, 0.0, 0.46153846153846156],
    [0.5, 0.5, 0.5384615384615385],
    [0.0, 0.0, 0.6153846153846154],
    [0.5, 0.5, 0.6923076923076923],
    [0.0, 0.0, 0.7692307692307693],
    [0.5, 0.5, 0.8461538461538463],
    [0.0, 0.0, 0.9230769230769231],
]

bulk_configuration_gold = BulkConfiguration(bravais_lattice=lattice,elements=elements_gold,fractional_coordinates=fractional_coordinates_gold)

vector_a = [3.1604, 0.0, 0.0]*Angstrom
vector_b = [0.0, 5.47397337224068, 0.0]*Angstrom
vector_c = [0.0, 0.0, 12.295]*Angstrom
lattice = UnitCell(vector_a, vector_b, vector_c)

elements_mos2 = [
Molybdenum,
Sulfur,
Sulfur,
Molybdenum,
Sulfur,
Sulfur,
]

fractional_coordinates_mos2 = [
    [-0.0, 0.666666666667, 0.25],
    [-0.0, 0.333333333333, 0.121],
    [-0.0, 0.333333333333, 0.379],
    [0.5, 0.166666666667, 0.25],
    [0.5, 0.833333333333, 0.121],
    [0.5, 0.833333333333, 0.379],
]

bulk_configuration_mo = BulkConfiguration(bravais_lattice=lattice,elements=elements_mos2,fractional_coordinates=fractional_coordinates_mos2)

interface_builder = InterfaceBuilder(configuration_0=bulk_configuration_gold,configuration_1=bulk_configuration_mo)

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
