"""
    ConfigurationStruct

Encapsulates the structural and physical description of a single polymer or molecular species
used in the simulation. This includes its sequence, allowed states, topology, and associated
physics model.

Fields:
- `species_model::AbstractSpecies`: Physics model used to evaluate partition functions,
  densities, and chemical potentials (e.g., `freely_jointed`, `monomerbead`).

- `sequence::Vector{Int}`: Linear list of monomer indices, referencing entries from
  `used_monomers`. Each entry represents the base monomer type for that bead.

- `states::Vector{Int}`: Number of internal states for each bead. `1` indicates the monomer
  is fixed (no protonation or switching); `>1` indicates participation in a state family
  (e.g., acid/base transitions).

- `topology::TopologyStruct`: Graph-like description of how beads are connected,
  parsed from a SMILES-like or shorthand string. Supports branched and linear architectures.

- `state_family::Vector{Vector{Int}}`: For each bead, the list of allowable monomer
  states by index. These are derived from protonation/deprotonation or other chemical
  transformation networks, and correspond to entries in the full monomer list.

This struct is processed during system setup and passed to density and free energy
evaluators throughout the simulation.
"""


"""
    AbstractSpecies

An abstract type that serves as the parent for all species-level chain models 
(e.g., `freely_jointed`, `monomerbead`, etc.). Concrete species models should 
subtype this to provide customized behavior for chain configurations.
"""

abstract type AbstractSpecies end


struct ConfigurationStruct
    species_model :: AbstractSpecies        # Physics model for the chain (e.g., freely_jointed, monomerbead)
    sequence      :: Vector{Int}            # Monomer indices for each bead, from used_monomers
    states        :: Vector{Int}            # Number of states per bead (1 = no switching, >1 = e.g. protonation)
    topology      :: TopologyStruct         # Parsed topology from sequence string (e.g., linear, branched)
    state_family  :: Vector{Vector{Int}}     # For each bead, list of possible monomer states (["A", "B", "C"])
end