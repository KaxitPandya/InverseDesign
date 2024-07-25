# -*- coding: utf-8 -*-
# -------------------------------------------------------------
# Bulk Configuration
# -------------------------------------------------------------

# Set up lattice
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
bulk_configuration = BulkConfiguration(
    bravais_lattice=lattice,
    elements=elements,
    fractional_coordinates=fractional_coordinates
    )
