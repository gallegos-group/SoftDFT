"""
    TopologyStruct

Represents the hierarchical connectivity and linearization of a polymer's architecture,
parsed from a structured sequence string (e.g., linear, branched, or nested topologies).

This struct encodes both the tree-like structure and the linear indexing required
for downstream calculations such as propagators, density evaluation, and bond enumeration.

# Fields
- `monomers::String`: Sequence of monomer names in linear order (e.g., `"AAABBB"`).
- `indices::Vector{Int}`: Mapping from tree-based node order to linear bead index.
- `children::Vector{Vector{Int}}`: For each node, list of child segment indices (tree topology).
- `parents::Vector{Vector{Int}}`: For each node, list of parent segment indices.
- `levels::Vector{Vector{Int}}`: Grouping of segment indices by level (distance from root).
- `sequ_bonds::Vector{Tuple{Int, Int}}`: List of sequential (tangent) bonds between beads,
  specified as (parent, child) index pairs in the linear sequence.

This structure supports both linear and complex branched chains and is constructed during
molecular system setup from user-defined configuration strings.
"""


struct TopologyStruct
    monomers :: String                      # Monomer names by linear index
    indices :: Vector{Int}                  # Tree -> linear index
    children :: Vector{Vector{Int}}         # Children of each node
    parents :: Vector{Vector{Int}}          # Parents of each node
    levels :: Vector{Vector{Int}}           # Nodes at each level
    sequ_bonds :: Vector{Tuple{Int, Int}}   # Sequential (tangent) bonds
end