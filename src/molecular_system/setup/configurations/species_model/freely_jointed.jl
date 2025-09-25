"""
    freely_jointed <: AbstractSpecies

Represents a coarse-grained freely-jointed chain species_model. Each bead has a diameter, and the
maximum allowed distance between any two segments is computed from pairwise diameter sums
along the shortest path connecting them in the bonding topology.
"""

struct freely_jointed <: AbstractSpecies
    max_distance :: Matrix{Float64}

    function freely_jointed(numbeads)
        new(zeros(Float64, (numbeads, numbeads)))
    end
end

"""
    bonding_constraint(config, properties, species_model::freely_jointed)

For a freely-jointed chain, computes the maximum allowed pairwise distances between all 
beads based on their sequential or branched connectivity. Uses a bidirectional breadth-first 
search (BFS) on the polymer topology to determine shortest-path distances between beads.
The result is stored in `species_model.max_distance`.
"""

function bonding_constraint(config, properties, species_model :: freely_jointed)
    if !haskey(properties.monomers, :diameters)
        error("Missing 'diameters' in properties.monomers. This is required for freely_jointed species_model.")
    end

    diameters = properties.monomers[:diameters]
    distance_matrix = max_distance_matrix(config, diameters)
    @. species_model.max_distance = distance_matrix
end

"""
    bidirectional_bfs_max_distances(start, config, diameters) -> Vector{Float64}

Performs a bidirectional breadth-first search on the topology of the chain starting from a
given segment. Computes cumulative distances from `start` to all other segments using average 
diameter-based bond lengths. Returns a vector of maximum distances from `start` to each bead.
"""

# Function to perform a bidirectional BFS to capture max distances from `start` to other segments
function bidirectional_bfs_max_distances(start::Int, config, diameters)
    topology = config.topology
    sequence = config.sequence

    n = length(sequence)
    distances = fill(0.0, n)  # Initialize distances from `start` to all other segments
    queue = [(start, 0.0)]    # Queue for BFS with initial distance set to 0.0
    visited = Set([start])    # Track visited segments

    while !isempty(queue)
        (current, dist) = popfirst!(queue)
        
        # Traverse to each child of the current segment
        for child in topology.children[current]
            if child ∉ visited && child <= n
                v = sequence[current]
                v1 = sequence[child]
            
                dij = (diameters[v] + diameters[v1]) / 2.0
                
                new_dist = dist + dij
                distances[child] = max(distances[child], new_dist)
                push!(queue, (child, new_dist))
                push!(visited, child)
            end
        end

        # Traverse to each parent of the current segment
        for parent in topology.parents[current]
            if parent ∉ visited && parent <= n
                v = sequence[current]
                v1 = sequence[parent]
            
                dij = (diameters[v] + diameters[v1]) / 2.0
                
                new_dist = dist + dij
                distances[parent] = max(distances[parent], new_dist)
                push!(queue, (parent, new_dist))
                push!(visited, parent)
            end
        end
    end
    return distances
end

"""
    max_distance_matrix(config, diameters) -> Matrix{Float64}

Constructs a full matrix of maximum allowed distances between all pairs of segments in the chain.
Each entry (i, j) is the maximum path length between bead i and bead j, calculated using 
`bidirectional_bfs_max_distances`.
"""

# Function to build the maximum distance matrix for all pairs of segments
function max_distance_matrix(config, diameters)
    n = length(config.sequence)
    distance_matrix = zeros(Float64, n, n)

    for i in 1:n
        distances_from_i = bidirectional_bfs_max_distances(i, config, diameters)
        distance_matrix[i, :] .= distances_from_i  # Store distances from segment i to all others
    end

    for i in 1:n
        distance_matrix[i, i] = 0.0
    end

    return distance_matrix
end

