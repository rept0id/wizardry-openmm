from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout

# Load in the PDB strucure
pdb = PDBFile('villin.pdb')

# Specifiy the forcefield
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

# Combine the molecular topology and the forcefield
system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME,
        nonbondedCutoff=1*nanometer, constraints=HBonds)

# Create the integrator to use for advacing the equations of motion.
# It specifies a Langevin integrator.
# The paramters set are temperature, friction coefficient, and timestep.
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)

# Combines the molecular topology, system, and integrator
# to begin a new simulation.
simulation = Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)

# Perform local energy minimization
print("Minimizing energy...")
simulation.minimizeEnergy(maxIterations=100)

# Write the trajectory to a file called "res.pdb"
simulation.reporters.append(PDBReporter('res.pdb', 1000))

# Report infomation to the screen as the simulation runs
simulation.reporters.append(StateDataReporter(stdout, 100, step=True,
        potentialEnergy=True, temperature=True))

# Run the simulation for 1000 timsteps
print("Running simulation...")
simulation.step(1000)