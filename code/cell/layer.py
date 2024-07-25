# Set up lattice
vector_a = [2.8837582303740374, 0.0, 0.0]*Angstrom
vector_b = [0.0, 5.767516460748075, 0.0]*Angstrom
vector_c = [0.0, 0.0, 22.490624999999994]*Angstrom
lattice = UnitCell(vector_a, vector_b, vector_c)

# Define elements
elements = [Gold, Gold, Gold, Gold, Gold, Gold, Gold, Gold, Gold, Gold, Sulfur,
            Sulfur, Molybdenum, Molybdenum, Sulfur, Sulfur]

# Define coordinates
fractional_coordinates = [[ 0.            ,  0.            ,  0.090665555092],
                          [ 0.            ,  0.5           ,  0.090665555092],
                          [ 0.5           ,  0.25          ,  0.181331110185],
                          [ 0.5           ,  0.75          ,  0.181331110185],
                          [ 0.            ,  0.            ,  0.271996665277],
                          [ 0.            ,  0.5           ,  0.271996665277],
                          [ 0.5           ,  0.25          ,  0.36266222037 ],
                          [ 0.5           ,  0.75          ,  0.36266222037 ],
                          [ 0.            ,  0.            ,  0.453327775462],
                          [ 0.            ,  0.5           ,  0.453327775462],
                          [ 0.            ,  0.2498049769  ,  0.52668958572 ],
                          [ 0.5           ,  0.7498049769  ,  0.52668958572 ],
                          [ 0.5           ,  0.083138310234,  0.597210302685],
                          [ 0.            ,  0.583138310234,  0.597210302685],
                          [ 0.            ,  0.2498049769  ,  0.66773101965 ],
                          [ 0.5           ,  0.7498049769  ,  0.66773101965 ]]

# Set up configuration
bulk_configuration = BulkConfiguration(
    bravais_lattice=lattice,
    elements=elements,
    fractional_coordinates=fractional_coordinates
    )

# Add tags
bulk_configuration.addTags('Left Interface 0',  [0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
bulk_configuration.addTags('Right Interface 0', [10, 11, 12, 13, 14, 15])
