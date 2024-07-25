vector_a = [-12.68095444386112, 0.0, 0.0]*Angstrom
vector_b = [0.0, -43.73381263022966, 0.0]*Angstrom
vector_c = [0.0, 0.0, 5.537625098599001]*Angstrom
left_electrode_lattice = UnitCell(vector_a, vector_b, vector_c)

left_electrode_elements = [
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
    Gold,
    Gold,
    Gold,
    Gold,
    Gold,
    Gold,
    Molybdenum,
    Molybdenum,
    Molybdenum,
    Molybdenum,
    Sulfur,
    Sulfur,
    Sulfur,
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
    Gold,
    Gold,
    Gold,
    Molybdenum,
    Molybdenum,
    Molybdenum,
    Molybdenum,
]

left_electrode_coordinates = [
    [0.0, -6.244604210066, 0.230704462761],
    [-4.226984814616, -6.244604210066, 0.230704462761],
    [-8.453969629245, -6.244604210066, 0.230704462761],
    [-2.113492407314, -8.326138946784, 0.230704462761],
    [-6.340477221931, -8.326138946784, 0.230704462761],
    [-10.567462036547, -8.326138946784, 0.230704462761],
    [0.0, -10.407673683459, 0.230704462761],
    [-4.226984814616, -10.407673683459, 0.230704462761],
    [-8.453969629245, -10.407673683459, 0.230704462761],
    [-2.113492407314, -12.489208420133, 0.230704462761],
    [-6.340477221931, -12.489208420133, 0.230704462761],
    [-10.567462036547, -12.489208420133, 0.230704462761],
    [0.0, -14.570743156851, 0.230704462761],
    [-4.226984814616, -14.570743156851, 0.230704462761],
    [-8.453969629245, -14.570743156851, 0.230704462761],
    [-2.113492407314, -16.652277893525, 0.230704462761],
    [-6.340477221931, -16.652277893525, 0.230704462761],
    [-10.567462036547, -16.652277893525, 0.230704462761],
    [-2.113492407314, -25.499999999982, 0.230704462761],
    [-6.340477221931, -25.499999999982, 0.230704462761],
    [-10.567462036547, -25.499999999982, 0.230704462761],
    [0.0, -27.5815347367, 0.230704462761],
    [-4.226984814616, -27.5815347367, 0.230704462761],
    [-8.453969629245, -27.5815347367, 0.230704462761],
    [-2.113492407314, -29.663069473374, 0.230704462761],
    [-6.340477221931, -29.663069473374, 0.230704462761],
    [-10.567462036547, -29.663069473374, 0.230704462761],
    [-2.113492407314, -6.244604210066, 2.307313874738],
    [-6.340477221931, -6.244604210066, 2.307313874738],
    [-10.567462036547, -6.244604210066, 2.307313874738],
    [0.0, -8.326138946784, 2.307313874738],
    [-4.226984814616, -8.326138946784, 2.307313874738],
    [-8.453969629245, -8.326138946784, 2.307313874738],
    [-2.113492407314, -10.407673683459, 2.307313874738],
    [-6.340477221931, -10.407673683459, 2.307313874738],
    [-10.567462036547, -10.407673683459, 2.307313874738],
    [0.0, -12.489208420133, 2.307313874738],
    [-4.226984814616, -12.489208420133, 2.307313874738],
    [-8.453969629245, -12.489208420133, 2.307313874738],
    [-2.113492407314, -14.570743156851, 2.307313874738],
    [-6.340477221931, -14.570743156851, 2.307313874738],
    [-10.567462036547, -14.570743156851, 2.307313874738],
    [0.0, -16.652277893525, 2.307313874738],
    [-4.226984814616, -16.652277893525, 2.307313874738],
    [-8.453969629245, -16.652277893525, 2.307313874738],
    [0.0, -25.499999999982, 2.307313874738],
    [-4.226984814616, -25.499999999982, 2.307313874738],
    [-8.453969629245, -25.499999999982, 2.307313874738],
    [-2.113492407314, -27.5815347367, 2.307313874738],
    [-6.340477221931, -27.5815347367, 2.307313874738],
    [-10.567462036547, -27.5815347367, 2.307313874738],
    [0.0, -29.663069473374, 2.307313874738],
    [-4.226984814616, -29.663069473374, 2.307313874738],
    [-8.453969629245, -29.663069473374, 2.307313874738],
    [-2.376514508405, -21.233893995098, 2.538025721855],
    [-5.54675311937, -21.233893995098, 2.538025721855],
    [-8.716991730336, -21.233893995098, 2.538025721855],
    [-11.887230341301, -21.233893995098, 2.538025721855],
    [-3.961558571876, -22.825909323972, 3.461533033405],
    [-7.131797182841, -22.825909323972, 3.461533033405],
    [-10.302035793806, -22.825909323972, 3.461533033405],
    [0.0, -6.244604210066, 4.383923286714],
    [-4.226984814616, -6.244604210066, 4.383923286714],
    [-8.453969629245, -6.244604210066, 4.383923286714],
    [-2.113492407314, -8.326138946784, 4.383923286714],
    [-6.340477221931, -8.326138946784, 4.383923286714],
    [-10.567462036547, -8.326138946784, 4.383923286714],
    [0.0, -10.407673683459, 4.383923286714],
    [-4.226984814616, -10.407673683459, 4.383923286714],
    [-8.453969629245, -10.407673683459, 4.383923286714],
    [-2.113492407314, -12.489208420133, 4.383923286714],
    [-6.340477221931, -12.489208420133, 4.383923286714],
    [-10.567462036547, -12.489208420133, 4.383923286714],
    [0.0, -14.570743156851, 4.383923286714],
    [-4.226984814616, -14.570743156851, 4.383923286714],
    [-8.453969629245, -14.570743156851, 4.383923286714],
    [-2.113492407314, -16.652277893525, 4.383923286714],
    [-6.340477221931, -16.652277893525, 4.383923286714],
    [-10.567462036547, -16.652277893525, 4.383923286714],
    [-2.113492407314, -25.499999999982, 4.383923286714],
    [-6.340477221931, -25.499999999982, 4.383923286714],
    [-10.567462036547, -25.499999999982, 4.383923286714],
    [0.0, -27.5815347367, 4.383923286714],
    [-4.226984814616, -27.5815347367, 4.383923286714],
    [-8.453969629245, -27.5815347367, 4.383923286714],
    [-2.113492407314, -29.663069473374, 4.383923286714],
    [-6.340477221931, -29.663069473374, 4.383923286714],
    [-10.567462036547, -29.663069473374, 4.383923286714],
    [-0.793784767511, -21.23385061496, 5.306920635838],
    [-3.964023378476, -21.23385061496, 5.306920635838],
    [-7.134261989441, -21.23385061496, 5.306920635838],
    [-10.304500600407, -21.23385061496, 5.306920635838],
]

left_electrode = BulkConfiguration(bravais_lattice=left_electrode_lattice,elements=left_electrode_elements,cartesian_coordinates=left_electrode_coordinates)

vector_a = [-12.68095444386112, 0.0, 0.0]
vector_b = [0.0, -43.73381263022966, 0.0]
vector_c = [0.0, 0.0, 4.153218823954]
right_electrode_lattice = UnitCell(vector_a, vector_b, vector_c)

right_electrode_elements = [
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
    Gold,
    Gold,
    Gold,
    Molybdenum,
    Molybdenum,
    Molybdenum,
    Molybdenum,
    Sulfur,
    Sulfur,
    Sulfur,
    Sulfur,
    Sulfur,
    Sulfur,
    Sulfur,
    Sulfur,
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
    Gold,
    Gold,
    Gold,
    Molybdenum,
    Molybdenum,
    Molybdenum,
    Molybdenum,
]

right_electrode_coordinates = [
    [0.0, -6.244604210066, 0.230704462769],
    [-4.226984814616, -6.244604210066, 0.230704462769],
    [-8.453969629245, -6.244604210066, 0.230704462769],
    [-2.113492407314, -8.326138946784, 0.230704462769],
    [-6.340477221931, -8.326138946784, 0.230704462769],
    [-10.567462036547, -8.326138946784, 0.230704462769],
    [0.0, -10.407673683459, 0.230704462769],
    [-4.226984814616, -10.407673683459, 0.230704462769],
    [-8.453969629245, -10.407673683459, 0.230704462769],
    [-2.113492407314, -12.489208420133, 0.230704462769],
    [-6.340477221931, -12.489208420133, 0.230704462769],
    [-10.567462036547, -12.489208420133, 0.230704462769],
    [0.0, -14.570743156851, 0.230704462769],
    [-4.226984814616, -14.570743156851, 0.230704462769],
    [-8.453969629245, -14.570743156851, 0.230704462769],
    [-2.113492407314, -16.652277893525, 0.230704462769],
    [-6.340477221931, -16.652277893525, 0.230704462769],
    [-10.567462036547, -16.652277893525, 0.230704462769],
    [-2.113492407314, -25.499999999982, 0.230704462769],
    [-6.340477221931, -25.499999999982, 0.230704462769],
    [-10.567462036547, -25.499999999982, 0.230704462769],
    [0.0, -27.5815347367, 0.230704462769],
    [-4.226984814616, -27.5815347367, 0.230704462769],
    [-8.453969629245, -27.5815347367, 0.230704462769],
    [-2.113492407314, -29.663069473374, 0.230704462769],
    [-6.340477221931, -29.663069473374, 0.230704462769],
    [-10.567462036547, -29.663069473374, 0.230704462769],
    [-2.376514508405, -21.233893995098, 1.153619447218],
    [-5.54675311937, -21.233893995098, 1.153619447218],
    [-8.716991730336, -21.233893995098, 1.153619447218],
    [-11.887230341301, -21.233893995098, 1.153619447218],
    [-0.791293867526, -19.641661826939, 2.077100396344],
    [-3.961532478491, -19.641661826939, 2.077100396344],
    [-7.131771089456, -19.641661826939, 2.077100396344],
    [-10.302009700422, -19.641661826939, 2.077100396344],
    [-0.79131996091, -22.825909323972, 2.077126758768],
    [-3.961558571876, -22.825909323972, 2.077126758768],
    [-7.131797182841, -22.825909323972, 2.077126758768],
    [-10.302035793806, -22.825909323972, 2.077126758768],
    [-2.113492407314, -6.244604210066, 2.307313874746],
    [-6.340477221931, -6.244604210066, 2.307313874746],
    [-10.567462036547, -6.244604210066, 2.307313874746],
    [0.0, -8.326138946784, 2.307313874746],
    [-4.226984814616, -8.326138946784, 2.307313874746],
    [-8.453969629245, -8.326138946784, 2.307313874746],
    [-2.113492407314, -10.407673683459, 2.307313874746],
    [-6.340477221931, -10.407673683459, 2.307313874746],
    [-10.567462036547, -10.407673683459, 2.307313874746],
    [0.0, -12.489208420133, 2.307313874746],
    [-4.226984814616, -12.489208420133, 2.307313874746],
    [-8.453969629245, -12.489208420133, 2.307313874746],
    [-2.113492407314, -14.570743156851, 2.307313874746],
    [-6.340477221931, -14.570743156851, 2.307313874746],
    [-10.567462036547, -14.570743156851, 2.307313874746],
    [0.0, -16.652277893525, 2.307313874746],
    [-4.226984814616, -16.652277893525, 2.307313874746],
    [-8.453969629245, -16.652277893525, 2.307313874746],
    [0.0, -25.499999999982, 2.307313874746],
    [-4.226984814616, -25.499999999982, 2.307313874746],
    [-8.453969629245, -25.499999999982, 2.307313874746],
    [-2.113492407314, -27.5815347367, 2.307313874746],
    [-6.340477221931, -27.5815347367, 2.307313874746],
    [-10.567462036547, -27.5815347367, 2.307313874746],
    [0.0, -29.663069473374, 2.307313874746],
    [-4.226984814616, -29.663069473374, 2.307313874746],
    [-8.453969629245, -29.663069473374, 2.307313874746],
    [-0.793784767511, -21.23385061496, 3.922514361185],
    [-3.964023378476, -21.23385061496, 3.922514361185],
    [-7.134261989441, -21.23385061496, 3.922514361185],
    [-10.304500600407, -21.23385061496, 3.922514361185],
]

right_electrode = BulkConfiguration(bravais_lattice=right_electrode_lattice,elements=right_electrode_elements,cartesian_coordinates=right_electrode_coordinates)

vector_a = [-12.68095444386112, 0.0, 0.0]
vector_b = [0.0, -43.73381263022966, 0.0]
vector_c = [0.0, 0.0, 16.6128752958063]
central_region_lattice = UnitCell(vector_a, vector_b, vector_c)

central_region_elements = [
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
    Gold,
    Gold,
    Gold,
    Sulfur,
    Sulfur,
    Sulfur,
    Sulfur,
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
    Gold,
    Gold,
    Gold,
    Molybdenum,
    Molybdenum,
    Molybdenum,
    Molybdenum,
    Sulfur,
    Sulfur,
    Sulfur,
    Sulfur,
    Sulfur,
    Sulfur,
    Sulfur,
    Sulfur,
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
    Gold,
    Gold,
    Gold,
    Molybdenum,
    Molybdenum,
    Molybdenum,
    Molybdenum,
    Sulfur,
    Sulfur,
    Sulfur,
    Sulfur,
    Sulfur,
    Sulfur,
    Sulfur,
    Sulfur,
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
    Gold,
    Gold,
    Gold,
    Molybdenum,
    Molybdenum,
    Molybdenum,
    Molybdenum,
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
    Gold,
    Gold,
    Gold,
    Sulfur,
    Sulfur,
    Sulfur,
    Sulfur,
    Sulfur,
    Sulfur,
    Sulfur,
    Sulfur,
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
    Gold,
    Gold,
    Gold,
    Molybdenum,
    Molybdenum,
    Molybdenum,
    Molybdenum,
    Sulfur,
    Sulfur,
    Sulfur,
    Sulfur,
    Sulfur,
    Sulfur,
    Sulfur,
    Sulfur,
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
    Gold,
    Gold,
    Gold,
    Molybdenum,
    Molybdenum,
    Molybdenum,
    Molybdenum,
    Sulfur,
    Sulfur,
    Sulfur,
    Sulfur,
    Sulfur,
    Sulfur,
    Sulfur,
    Sulfur,
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
    Gold,
    Gold,
    Gold,
    Molybdenum,
    Molybdenum,
    Molybdenum,
    Molybdenum,
]

central_region_coordinates = [
    [0.0, -6.244604210066, 0.230704462761],
    [-4.226984814616, -6.244604210066, 0.230704462761],
    [-8.453969629245, -6.244604210066, 0.230704462761],
    [-2.113492407314, -8.326138946784, 0.230704462761],
    [-6.340477221931, -8.326138946784, 0.230704462761],
    [-10.567462036547, -8.326138946784, 0.230704462761],
    [0.0, -10.407673683459, 0.230704462761],
    [-4.226984814616, -10.407673683459, 0.230704462761],
    [-8.453969629245, -10.407673683459, 0.230704462761],
    [-2.113492407314, -12.489208420133, 0.230704462761],
    [-6.340477221931, -12.489208420133, 0.230704462761],
    [-10.567462036547, -12.489208420133, 0.230704462761],
    [0.0, -14.570743156851, 0.230704462761],
    [-4.226984814616, -14.570743156851, 0.230704462761],
    [-8.453969629245, -14.570743156851, 0.230704462761],
    [-2.113492407314, -16.652277893525, 0.230704462761],
    [-6.340477221931, -16.652277893525, 0.230704462761],
    [-10.567462036547, -16.652277893525, 0.230704462761],
    [-2.113492407314, -25.499999999982, 0.230704462761],
    [-6.340477221931, -25.499999999982, 0.230704462761],
    [-10.567462036547, -25.499999999982, 0.230704462761],
    [0.0, -27.5815347367, 0.230704462761],
    [-4.226984814616, -27.5815347367, 0.230704462761],
    [-8.453969629245, -27.5815347367, 0.230704462761],
    [-2.113492407314, -29.663069473374, 0.230704462761],
    [-6.340477221931, -29.663069473374, 0.230704462761],
    [-10.567462036547, -29.663069473374, 0.230704462761],
    [-2.378992903244, -22.826072019441, 0.691548894542],
    [-5.54923151421, -22.826072019441, 0.691548894542],
    [-8.719470125175, -22.826072019441, 0.691548894542],
    [-11.88970873614, -22.826072019441, 0.691548894542],
    [-2.113492407314, -6.244604210066, 2.307313874738],
    [-6.340477221931, -6.244604210066, 2.307313874738],
    [-10.567462036547, -6.244604210066, 2.307313874738],
    [0.0, -8.326138946784, 2.307313874738],
    [-4.226984814616, -8.326138946784, 2.307313874738],
    [-8.453969629245, -8.326138946784, 2.307313874738],
    [-2.113492407314, -10.407673683459, 2.307313874738],
    [-6.340477221931, -10.407673683459, 2.307313874738],
    [-10.567462036547, -10.407673683459, 2.307313874738],
    [0.0, -12.489208420133, 2.307313874738],
    [-4.226984814616, -12.489208420133, 2.307313874738],
    [-8.453969629245, -12.489208420133, 2.307313874738],
    [-2.113492407314, -14.570743156851, 2.307313874738],
    [-6.340477221931, -14.570743156851, 2.307313874738],
    [-10.567462036547, -14.570743156851, 2.307313874738],
    [0.0, -16.652277893525, 2.307313874738],
    [-4.226984814616, -16.652277893525, 2.307313874738],
    [-8.453969629245, -16.652277893525, 2.307313874738],
    [0.0, -25.499999999982, 2.307313874738],
    [-4.226984814616, -25.499999999982, 2.307313874738],
    [-8.453969629245, -25.499999999982, 2.307313874738],
    [-2.113492407314, -27.5815347367, 2.307313874738],
    [-6.340477221931, -27.5815347367, 2.307313874738],
    [-10.567462036547, -27.5815347367, 2.307313874738],
    [0.0, -29.663069473374, 2.307313874738],
    [-4.226984814616, -29.663069473374, 2.307313874738],
    [-8.453969629245, -29.663069473374, 2.307313874738],
    [-2.376514508405, -21.233893995098, 2.538025721855],
    [-5.54675311937, -21.233893995098, 2.538025721855],
    [-8.716991730336, -21.233893995098, 2.538025721855],
    [-11.887230341301, -21.233893995098, 2.538025721855],
    [-0.791293867526, -19.641661826939, 3.461506670981],
    [-3.961532478491, -19.641661826939, 3.461506670981],
    [-7.131771089456, -19.641661826939, 3.461506670981],
    [-10.302009700422, -19.641661826939, 3.461506670981],
    [-0.79131996091, -22.825909323972, 3.461533033405],
    [-3.961558571876, -22.825909323972, 3.461533033405],
    [-7.131797182841, -22.825909323972, 3.461533033405],
    [-10.302035793806, -22.825909323972, 3.461533033405],
    [0.0, -6.244604210066, 4.383923286714],
    [-4.226984814616, -6.244604210066, 4.383923286714],
    [-8.453969629245, -6.244604210066, 4.383923286714],
    [-2.113492407314, -8.326138946784, 4.383923286714],
    [-6.340477221931, -8.326138946784, 4.383923286714],
    [-10.567462036547, -8.326138946784, 4.383923286714],
    [0.0, -10.407673683459, 4.383923286714],
    [-4.226984814616, -10.407673683459, 4.383923286714],
    [-8.453969629245, -10.407673683459, 4.383923286714],
    [-2.113492407314, -12.489208420133, 4.383923286714],
    [-6.340477221931, -12.489208420133, 4.383923286714],
    [-10.567462036547, -12.489208420133, 4.383923286714],
    [0.0, -14.570743156851, 4.383923286714],
    [-4.226984814616, -14.570743156851, 4.383923286714],
    [-8.453969629245, -14.570743156851, 4.383923286714],
    [-2.113492407314, -16.652277893525, 4.383923286714],
    [-6.340477221931, -16.652277893525, 4.383923286714],
    [-10.567462036547, -16.652277893525, 4.383923286714],
    [-2.113492407314, -25.499999999982, 4.383923286714],
    [-6.340477221931, -25.499999999982, 4.383923286714],
    [-10.567462036547, -25.499999999982, 4.383923286714],
    [0.0, -27.5815347367, 4.383923286714],
    [-4.226984814616, -27.5815347367, 4.383923286714],
    [-8.453969629245, -27.5815347367, 4.383923286714],
    [-2.113492407314, -29.663069473374, 4.383923286714],
    [-6.340477221931, -29.663069473374, 4.383923286714],
    [-10.567462036547, -29.663069473374, 4.383923286714],
    [-0.793784767511, -21.23385061496, 5.306920635838],
    [-3.964023378476, -21.23385061496, 5.306920635838],
    [-7.134261989441, -21.23385061496, 5.306920635838],
    [-10.304500600407, -21.23385061496, 5.306920635838],
    [-2.378992903244, -22.826072019441, 6.229173993158],
    [-5.54923151421, -22.826072019441, 6.229173993158],
    [-8.719470125175, -22.826072019441, 6.229173993158],
    [-11.88970873614, -22.826072019441, 6.229173993158],
    [-2.378971914578, -19.641488000909, 6.229207776436],
    [-5.549210525543, -19.641488000909, 6.229207776436],
    [-8.719449136509, -19.641488000909, 6.229207776436],
    [-11.889687747474, -19.641488000909, 6.229207776436],
    [-2.113492407314, -6.244604210066, 6.460532698691],
    [-6.340477221931, -6.244604210066, 6.460532698691],
    [-10.567462036547, -6.244604210066, 6.460532698691],
    [0.0, -8.326138946784, 6.460532698691],
    [-4.226984814616, -8.326138946784, 6.460532698691],
    [-8.453969629245, -8.326138946784, 6.460532698691],
    [-2.113492407314, -10.407673683459, 6.460532698691],
    [-6.340477221931, -10.407673683459, 6.460532698691],
    [-10.567462036547, -10.407673683459, 6.460532698691],
    [0.0, -12.489208420133, 6.460532698691],
    [-4.226984814616, -12.489208420133, 6.460532698691],
    [-8.453969629245, -12.489208420133, 6.460532698691],
    [-2.113492407314, -14.570743156851, 6.460532698691],
    [-6.340477221931, -14.570743156851, 6.460532698691],
    [-10.567462036547, -14.570743156851, 6.460532698691],
    [0.0, -16.652277893525, 6.460532698691],
    [-4.226984814616, -16.652277893525, 6.460532698691],
    [-8.453969629245, -16.652277893525, 6.460532698691],
    [0.0, -25.499999999982, 6.460532698691],
    [-4.226984814616, -25.499999999982, 6.460532698691],
    [-8.453969629245, -25.499999999982, 6.460532698691],
    [-2.113492407314, -27.5815347367, 6.460532698691],
    [-6.340477221931, -27.5815347367, 6.460532698691],
    [-10.567462036547, -27.5815347367, 6.460532698691],
    [0.0, -29.663069473374, 6.460532698691],
    [-4.226984814616, -29.663069473374, 6.460532698691],
    [-8.453969629245, -29.663069473374, 6.460532698691],
    [-2.376514508405, -21.233893995098, 8.075650820454],
    [-5.54675311937, -21.233893995098, 8.075650820454],
    [-8.716991730336, -21.233893995098, 8.075650820454],
    [-11.887230341301, -21.233893995098, 8.075650820454],
    [0.0, -6.244604210066, 8.537142110668],
    [-4.226984814616, -6.244604210066, 8.537142110668],
    [-8.453969629245, -6.244604210066, 8.537142110668],
    [-2.113492407314, -8.326138946784, 8.537142110668],
    [-6.340477221931, -8.326138946784, 8.537142110668],
    [-10.567462036547, -8.326138946784, 8.537142110668],
    [0.0, -10.407673683459, 8.537142110668],
    [-4.226984814616, -10.407673683459, 8.537142110668],
    [-8.453969629245, -10.407673683459, 8.537142110668],
    [-2.113492407314, -12.489208420133, 8.537142110668],
    [-6.340477221931, -12.489208420133, 8.537142110668],
    [-10.567462036547, -12.489208420133, 8.537142110668],
    [0.0, -14.570743156851, 8.537142110668],
    [-4.226984814616, -14.570743156851, 8.537142110668],
    [-8.453969629245, -14.570743156851, 8.537142110668],
    [-2.113492407314, -16.652277893525, 8.537142110668],
    [-6.340477221931, -16.652277893525, 8.537142110668],
    [-10.567462036547, -16.652277893525, 8.537142110668],
    [-2.113492407314, -25.499999999982, 8.537142110668],
    [-6.340477221931, -25.499999999982, 8.537142110668],
    [-10.567462036547, -25.499999999982, 8.537142110668],
    [0.0, -27.5815347367, 8.537142110668],
    [-4.226984814616, -27.5815347367, 8.537142110668],
    [-8.453969629245, -27.5815347367, 8.537142110668],
    [-2.113492407314, -29.663069473374, 8.537142110668],
    [-6.340477221931, -29.663069473374, 8.537142110668],
    [-10.567462036547, -29.663069473374, 8.537142110668],
    [-0.791293867526, -19.641661826939, 8.999131769597],
    [-3.961532478491, -19.641661826939, 8.999131769597],
    [-7.131771089456, -19.641661826939, 8.999131769597],
    [-10.302009700422, -19.641661826939, 8.999131769597],
    [-0.79131996091, -22.825909323972, 8.999158132021],
    [-3.961558571876, -22.825909323972, 8.999158132021],
    [-7.131797182841, -22.825909323972, 8.999158132021],
    [-10.302035793806, -22.825909323972, 8.999158132021],
    [-2.113492407314, -6.244604210066, 10.613751522645],
    [-6.340477221931, -6.244604210066, 10.613751522645],
    [-10.567462036547, -6.244604210066, 10.613751522645],
    [0.0, -8.326138946784, 10.613751522645],
    [-4.226984814616, -8.326138946784, 10.613751522645],
    [-8.453969629245, -8.326138946784, 10.613751522645],
    [-2.113492407314, -10.407673683459, 10.613751522645],
    [-6.340477221931, -10.407673683459, 10.613751522645],
    [-10.567462036547, -10.407673683459, 10.613751522645],
    [0.0, -12.489208420133, 10.613751522645],
    [-4.226984814616, -12.489208420133, 10.613751522645],
    [-8.453969629245, -12.489208420133, 10.613751522645],
    [-2.113492407314, -14.570743156851, 10.613751522645],
    [-6.340477221931, -14.570743156851, 10.613751522645],
    [-10.567462036547, -14.570743156851, 10.613751522645],
    [0.0, -16.652277893525, 10.613751522645],
    [-4.226984814616, -16.652277893525, 10.613751522645],
    [-8.453969629245, -16.652277893525, 10.613751522645],
    [0.0, -25.499999999982, 10.613751522645],
    [-4.226984814616, -25.499999999982, 10.613751522645],
    [-8.453969629245, -25.499999999982, 10.613751522645],
    [-2.113492407314, -27.5815347367, 10.613751522645],
    [-6.340477221931, -27.5815347367, 10.613751522645],
    [-10.567462036547, -27.5815347367, 10.613751522645],
    [0.0, -29.663069473374, 10.613751522645],
    [-4.226984814616, -29.663069473374, 10.613751522645],
    [-8.453969629245, -29.663069473374, 10.613751522645],
    [-0.793784767511, -21.23385061496, 10.844545734438],
    [-3.964023378476, -21.23385061496, 10.844545734438],
    [-7.134261989441, -21.23385061496, 10.844545734438],
    [-10.304500600407, -21.23385061496, 10.844545734438],
    [-2.378992903244, -22.826072019441, 11.766799091757],
    [-5.54923151421, -22.826072019441, 11.766799091757],
    [-8.719470125175, -22.826072019441, 11.766799091757],
    [-11.88970873614, -22.826072019441, 11.766799091757],
    [-2.378971914578, -19.641488000909, 11.766832875036],
    [-5.549210525543, -19.641488000909, 11.766832875036],
    [-8.719449136509, -19.641488000909, 11.766832875036],
    [-11.889687747474, -19.641488000909, 11.766832875036],
    [0.0, -6.244604210066, 12.690360934622],
    [-4.226984814616, -6.244604210066, 12.690360934622],
    [-8.453969629245, -6.244604210066, 12.690360934622],
    [-2.113492407314, -8.326138946784, 12.690360934622],
    [-6.340477221931, -8.326138946784, 12.690360934622],
    [-10.567462036547, -8.326138946784, 12.690360934622],
    [0.0, -10.407673683459, 12.690360934622],
    [-4.226984814616, -10.407673683459, 12.690360934622],
    [-8.453969629245, -10.407673683459, 12.690360934622],
    [-2.113492407314, -12.489208420133, 12.690360934622],
    [-6.340477221931, -12.489208420133, 12.690360934622],
    [-10.567462036547, -12.489208420133, 12.690360934622],
    [0.0, -14.570743156851, 12.690360934622],
    [-4.226984814616, -14.570743156851, 12.690360934622],
    [-8.453969629245, -14.570743156851, 12.690360934622],
    [-2.113492407314, -16.652277893525, 12.690360934622],
    [-6.340477221931, -16.652277893525, 12.690360934622],
    [-10.567462036547, -16.652277893525, 12.690360934622],
    [-2.113492407314, -25.499999999982, 12.690360934622],
    [-6.340477221931, -25.499999999982, 12.690360934622],
    [-10.567462036547, -25.499999999982, 12.690360934622],
    [0.0, -27.5815347367, 12.690360934622],
    [-4.226984814616, -27.5815347367, 12.690360934622],
    [-8.453969629245, -27.5815347367, 12.690360934622],
    [-2.113492407314, -29.663069473374, 12.690360934622],
    [-6.340477221931, -29.663069473374, 12.690360934622],
    [-10.567462036547, -29.663069473374, 12.690360934622],
    [-2.376514508405, -21.233893995098, 13.61327591907],
    [-5.54675311937, -21.233893995098, 13.61327591907],
    [-8.716991730336, -21.233893995098, 13.61327591907],
    [-11.887230341301, -21.233893995098, 13.61327591907],
    [-0.791293867526, -19.641661826939, 14.536756868197],
    [-3.961532478491, -19.641661826939, 14.536756868197],
    [-7.131771089456, -19.641661826939, 14.536756868197],
    [-10.302009700422, -19.641661826939, 14.536756868197],
    [-0.79131996091, -22.825909323972, 14.53678323062],
    [-3.961558571876, -22.825909323972, 14.53678323062],
    [-7.131797182841, -22.825909323972, 14.53678323062],
    [-10.302035793806, -22.825909323972, 14.53678323062],
    [-2.113492407314, -6.244604210066, 14.766970346598],
    [-6.340477221931, -6.244604210066, 14.766970346598],
    [-10.567462036547, -6.244604210066, 14.766970346598],
    [0.0, -8.326138946784, 14.766970346598],
    [-4.226984814616, -8.326138946784, 14.766970346598],
    [-8.453969629245, -8.326138946784, 14.766970346598],
    [-2.113492407314, -10.407673683459, 14.766970346598],
    [-6.340477221931, -10.407673683459, 14.766970346598],
    [-10.567462036547, -10.407673683459, 14.766970346598],
    [0.0, -12.489208420133, 14.766970346598],
    [-4.226984814616, -12.489208420133, 14.766970346598],
    [-8.453969629245, -12.489208420133, 14.766970346598],
    [-2.113492407314, -14.570743156851, 14.766970346598],
    [-6.340477221931, -14.570743156851, 14.766970346598],
    [-10.567462036547, -14.570743156851, 14.766970346598],
    [0.0, -16.652277893525, 14.766970346598],
    [-4.226984814616, -16.652277893525, 14.766970346598],
    [-8.453969629245, -16.652277893525, 14.766970346598],
    [0.0, -25.499999999982, 14.766970346598],
    [-4.226984814616, -25.499999999982, 14.766970346598],
    [-8.453969629245, -25.499999999982, 14.766970346598],
    [-2.113492407314, -27.5815347367, 14.766970346598],
    [-6.340477221931, -27.5815347367, 14.766970346598],
    [-10.567462036547, -27.5815347367, 14.766970346598],
    [0.0, -29.663069473374, 14.766970346598],
    [-4.226984814616, -29.663069473374, 14.766970346598],
    [-8.453969629245, -29.663069473374, 14.766970346598],
    [-0.793784767511, -21.23385061496, 16.382170833037],
    [-3.964023378476, -21.23385061496, 16.382170833037],
    [-7.134261989441, -21.23385061496, 16.382170833037],
    [-10.304500600407, -21.23385061496, 16.382170833037],
]

central_region = BulkConfiguration(bravais_lattice=central_region_lattice,elements=central_region_elements,    cartesian_coordinates=central_region_coordinates)

