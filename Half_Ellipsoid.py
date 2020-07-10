#!/usr/bin/env python
import os
import numpy as np
import time
import placentagen as pg
from opencmiss.iron import iron

## Definition of the sh4ape of the mesh
volume = 4.7e5  # mm^3 #volume of placenta
thickness = 22.1  # mm #thickness of placenta (z-axis dimension)
ellipticity = 1.00  # no units #ellipticity of placenta - ratio of y to x axis dimensions

# define elements in x-, y- and z-
size_el = 10  # mm #element size
circle_prop = 0.85  # maximum proportion of z surface taken up by circular mesh - helps to get good curvature in sides of mesh
squareSizeRatio = 0.50  # proportion of cross section comprised of aquare elements

el_type = 2  # Put 1 for linear element, 2 for quadratic elements....

# Options
# export generated mesh for visualisation
export_mesh = False
filename_mesh = 'expected-results/placenta_mesh'

export_results = True
export_directory = 'output'


dv_bc_type = 'velocity'  # or pressure
dv_bc_value = 1  # cm per sec outlet vel of decidual vein

sa_bc_type = 'velocity'  # or pressure
sa_bc_value = 10  # cm per sec inlet vel of spiral artery

wall_bc_type = 'no_slip'  # or no_penetration

debug = False
# Generate an ellipodial mesh for simulation, w
pl_mesh = pg.gen_half_ellipsoid_structured(size_el, volume, thickness, ellipticity, squareSizeRatio, circle_prop, el_type,
                                         debug)

# Optional - export mesh prior to solution
if (export_mesh):
    pg.export_ex_coords(pl_mesh['nodes'], 'placenta', filename_mesh, 'exnode')
    if (el_type == 1):
        pg.export_exelem_3d_linear(pl_mesh['placental_el_con'], 'placenta', filename_mesh)  # Use this for linear
    elif (el_type == 2):
        pg.export_exelem_3d_quadratic(pl_mesh['elems'], 'placenta', filename_mesh)  # Use this for quadratic
    else:
        print('Cannot export elements for visualisation as element type does not exist')


numberOfDimensions = 3
numberOfComponents = 1

# Material properties of this example.....
porosity = 0.4
perm = 0.0052 #cm2 permeability
viscosity = 3.36*1e-2 #poise...
perm_over_vis = perm/viscosity  # permeability over viscosity is
gamma = 1/(1+2.5*(1-porosity))
beta = viscosity/gamma
initial_velocity = 0.0

#Spiral artery properties....
# Inlet locations
stem_file = 'stem_xy.txt'
spiral_rad = 0.23#cm
dec_rad = 0.2
offset = 2.0#cm
cx1 = offset
cy1 = 0.0
cx2 = -offset
cy2 = 0.0

# Set problem parameters
(coordinateSystemUserNumber,
 regionUserNumber,
 basisUserNumber,
 generatedMeshUserNumber,
 meshUserNumber,
 decompositionUserNumber,
 geometricFieldUserNumber,
 equationsSetFieldUserNumber,
 dependentFieldUserNumber,
 materialFieldUserNumber,
 equationsSetUserNumber,
 problemUserNumber) = range(1, 13)

number_of_dimensions = 3
number_of_mesh_components = 1
total_number_of_elements = len(pl_mesh["elems"])
total_number_of_nodes = len(pl_mesh["nodes"])
mesh_component_number = 1

# Create a RC coordinate system
coordinateSystem = iron.CoordinateSystem()
coordinateSystem.CreateStart(coordinateSystemUserNumber)
coordinateSystem.dimension = 3
coordinateSystem.CreateFinish()

# Create a region
region = iron.Region()
region.CreateStart(regionUserNumber, iron.WorldRegion)
region.label = "DarcyRegion"
region.coordinateSystem = coordinateSystem
region.CreateFinish()

# Create a tri-quadratic Lagrange basis
basis = iron.Basis()
basis.CreateStart(basisUserNumber)
basis.TypeSet(iron.BasisTypes.LAGRANGE_HERMITE_TP)
basis.numberOfXi = 3
if (el_type == 2):
    basis.interpolationXi = [iron.BasisInterpolationSpecifications.QUADRATIC_LAGRANGE] * 3
    numberOfNodesXi = 3
    numberOfGaussXi = 3
elif (el_type == 1):
    basis.interpolationXi = [iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE] * 3
    numberOfNodesXi = 2
    numberOfGaussXi = 2
basis.QuadratureNumberOfGaussXiSet([numberOfGaussXi] * 3)
basis.quadratureLocalFaceGaussEvaluate = True
basis.CreateFinish()

# Start the creation of the imported mesh in the region
mesh = iron.Mesh()
mesh.CreateStart(meshUserNumber, region, number_of_dimensions)
mesh.NumberOfComponentsSet(number_of_mesh_components)
mesh.NumberOfElementsSet(total_number_of_elements)

# Define nodes for the mesh
nodes = iron.Nodes()
nodes.CreateStart(region, total_number_of_nodes)
# Refers to nodes by their user number as described in the original mesh
nodes.UserNumbersAllSet(pl_mesh['node_list'])
nodes.CreateFinish()

elements = iron.MeshElements()
elements.CreateStart(mesh, mesh_component_number, basis)

# Set the nodes pertaining to each element
for idx, elem_details in enumerate(pl_mesh['elems']):
    elem_details = elem_details.astype('int32')
    element_nodes_array = elem_details[1:28].astype('int32')
    elements.NodesSet(idx + 1, elem_details[1:28].astype('int32'))

# Refers to elements by their user number as described in the original mesh
elements.UserNumbersAllSet(pl_mesh['elems'][:, 0])
elements.CreateFinish()

mesh.CreateFinish()

# Dedine number of computational nodes and mesh properties
numberOfComputationalNodes = iron.ComputationalNumberOfNodesGet()
computationalNodeNumber = iron.ComputationalNodeNumberGet()
# Create a decomposition for the mesh
decomposition = iron.Decomposition()
decomposition.CreateStart(decompositionUserNumber, mesh)
decomposition.type = iron.DecompositionTypes.CALCULATED
decomposition.numberOfDomains = numberOfComputationalNodes
decomposition.CalculateFacesSet(True)
decomposition.CreateFinish()

# Create a field for the geometry
geometricField = iron.Field()
geometricField.CreateStart(geometricFieldUserNumber, region)
geometricField.meshDecomposition = decomposition
# Set the scaling to use
geometricField.ScalingTypeSet(iron.FieldScalingTypes.NONE)
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 1, 1)
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 2, 1)
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 3, 1)
geometricField.CreateFinish()

# Update the geometric field parameters
geometricField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)

for idx, node_locs in enumerate(pl_mesh["nodes"]):
    node_num = pl_mesh['node_list'][idx]
    [x, y, z] = node_locs[1:4]

    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, 1,
                                            int(node_num), 1, x)
    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, 1,
                                            int(node_num), 2, y)
    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, 1,
                                            int(node_num), 3, z)

geometricField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)

# Create standard Darcy equations set
equationsSetField = iron.Field()
equationsSet = iron.EquationsSet()
equationsSetSpecification = [iron.EquationsSetClasses.FLUID_MECHANICS,
                             iron.EquationsSetTypes.DARCY_EQUATION,
                             iron.EquationsSetSubtypes.DARCY_BRINKMAN]
equationsSet.CreateStart(equationsSetUserNumber, region, geometricField,
                         equationsSetSpecification, equationsSetFieldUserNumber, equationsSetField)
equationsSet.CreateFinish()

# Create dependent field
dependentField = iron.Field()
equationsSet.DependentCreateStart(dependentFieldUserNumber, dependentField)
dependentField.VariableLabelSet(iron.FieldVariableTypes.U, "Dependent")
dependentField.DOFOrderTypeSet(iron.FieldVariableTypes.U, iron.FieldDOFOrderTypes.SEPARATED)
dependentField.DOFOrderTypeSet(iron.FieldVariableTypes.DELUDELN, iron.FieldDOFOrderTypes.SEPARATED)
equationsSet.DependentCreateFinish()

# Create material field
materialField = iron.Field()
equationsSet.MaterialsCreateStart(materialFieldUserNumber, materialField)
materialField.VariableLabelSet(iron.FieldVariableTypes.U, "Material")
equationsSet.MaterialsCreateFinish()

materialField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, porosity)
materialField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 2,
                                          perm_over_vis)

for knode in range(0, len(pl_mesh['node_list'])):
    materialField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, 1,
                                           int(knode + 1), 1, porosity)  # set porosity
    materialField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, 1,
                                           int(knode + 1), 2, perm_over_vis)  # set perm_over_vis
    materialField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, 1,
                                           int(knode + 1), 3, beta)

# Initialise dependent field

dependentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1,
                                           initial_velocity)

# Create equations
equations = iron.Equations()
equationsSet.EquationsCreateStart(equations)
equations.sparsityType = iron.EquationsSparsityTypes.SPARSE
equations.outputType = iron.EquationsOutputTypes.NONE
equationsSet.EquationsCreateFinish()

# Create Darcy equation problem
problem = iron.Problem()
problemSpecification = [iron.ProblemClasses.FLUID_MECHANICS,
                        iron.ProblemTypes.DARCY_EQUATION,
                        iron.ProblemSubtypes.DARCY_BRINKMAN]
problem.CreateStart(problemUserNumber, problemSpecification)
problem.CreateFinish()

# Create control loops
problem.ControlLoopCreateStart()
problem.ControlLoopCreateFinish()

# Create problem solver
solver = iron.Solver()
problem.SolversCreateStart()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE], 1, solver)
# solver.outputType = iron.SolverOutputTypes.SOLVER
solver.outputType = iron.SolverOutputTypes.NONE
solver.linearType = iron.LinearSolverTypes.ITERATIVE
solver.linearIterativeAbsoluteTolerance = 1.0E-12
solver.linearIterativeRelativeTolerance = 1.0E-12
problem.SolversCreateFinish()

## Create solver equations and add equations set to solver equations
solver = iron.Solver()
solverEquations = iron.SolverEquations()
problem.SolverEquationsCreateStart()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE], 1, solver)
solver.SolverEquationsGet(solverEquations)
solverEquations.sparsityType = iron.SolverEquationsSparsityTypes.SPARSE
equationsSetIndex = solverEquations.EquationsSetAdd(equationsSet)
problem.SolverEquationsCreateFinish()

# Create boundary conditions
boundaryConditions = iron.BoundaryConditions()
solverEquations.BoundaryConditionsCreateStart(boundaryConditions)

Total_Nodes = nodes.NumberOfNodesGet()

# INLET BC
spiral = []
decidual1 = []
decidual2 = []
base_nodes = []
surface_nodes = []
#collect the nodes of spiral arteries and decidual veins....
for i in range(0,Total_Nodes):
    z_val = geometricField.ParameterSetGetNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,i,3)
    x_val = geometricField.ParameterSetGetNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,i,1)
    y_val = geometricField.ParameterSetGetNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,i,2)
    if(z_val == 0.0):
        r   = np.sqrt((x_val**2)+(y_val**2))
        rd1 = np.sqrt(((x_val-cx1)**2)+((y_val-cy1)**2))
        rd2 = np.sqrt(((x_val - cx2) ** 2) + ((y_val - cy2) ** 2))
        if(r<=spiral_rad):
            spiral.append(i)
        elif(rd1<=dec_rad):
            decidual1.append(i)
        elif(rd2<=dec_rad):
            decidual2.append(i)
        else:
            base_nodes.append(i)


# Assign the boundary conditions....
#Spiral artery.....
for i in range(0, len(spiral)):
    boundaryConditions.SetNode(dependentField, iron.FieldVariableTypes.U, 1, 1, spiral[i], 3,
                               iron.BoundaryConditionsTypes.FIXED, sa_bc_value)
    boundaryConditions.SetNode(dependentField, iron.FieldVariableTypes.U, 1, 1, spiral[i], 2,
                               iron.BoundaryConditionsTypes.FIXED, 0.0)
    boundaryConditions.SetNode(dependentField, iron.FieldVariableTypes.U, 1, 1, spiral[i], 1,
                               iron.BoundaryConditionsTypes.FIXED, 0.0)
dec = []
dec = decidual1+decidual2
for i in range(0,len(dec)):
    boundaryConditions.SetNode(dependentField, iron.FieldVariableTypes.U, 1, 1, dec[i], 3,
                               iron.BoundaryConditionsTypes.FIXED, -dv_bc_value)
    boundaryConditions.SetNode(dependentField, iron.FieldVariableTypes.U, 1, 1, dec[i], 2,
                               iron.BoundaryConditionsTypes.FIXED, 0.0)
    boundaryConditions.SetNode(dependentField, iron.FieldVariableTypes.U, 1, 1, dec[i], 1,
                               iron.BoundaryConditionsTypes.FIXED, 0.0)
for i in range(0,len(base_nodes)):
    boundaryConditions.SetNode(dependentField, iron.FieldVariableTypes.U, 1, 1, base_nodes[i], 3,
                               iron.BoundaryConditionsTypes.FIXED, 0.0)
    boundaryConditions.SetNode(dependentField, iron.FieldVariableTypes.U, 1, 1, base_nodes[i], 2,
                               iron.BoundaryConditionsTypes.FIXED, 0.0)
    boundaryConditions.SetNode(dependentField, iron.FieldVariableTypes.U, 1, 1, base_nodes[i], 1,
                               iron.BoundaryConditionsTypes.FIXED, 0.0)



solverEquations.BoundaryConditionsCreateFinish()
start_time = time.time()
# Solve the problem
problem.Solve()
end_time = time.time()
print('Total time for solve = ' + str((end_time - start_time) / 60.0) + ' mins')

if (export_results):
    if not os.path.exists(export_directory):
        os.makedirs(export_directory)
        ## Export results
    export_file = export_directory + '/StaticDarcy'
    fields = iron.Fields()
    fields.CreateRegion(region)
    fields.NodesExport(export_file, "FORTRAN")
    fields.ElementsExport(export_file, "FORTRAN")
    fields.Finalise()

iron.Finalise()
raise SystemExit