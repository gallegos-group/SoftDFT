"""
    process_bonds(configurations::Vector{ConfigurationStruct}) -> Vector{Tuple{Int, Int}}

Generates a sorted list of unique bonded bead-type pairs across all configurations
and their state variants. Each bond is a tuple `(i, j)` with `i ≤ j`, where `i` and `j` 
are bead type indices corresponding to specific chemical states.
"""
function process_bonds(configurations::Vector{ConfigurationStruct})::Vector{Tuple{Int, Int}}
    bond_set = Set{Tuple{Int, Int}}()

    for config in configurations
        @unpack sequence, state_family = config
        @unpack parents, levels = config.topology

        for level in levels, child in level
            for j in state_family[child]
                for parent in parents[child]
                    if !isempty(parent)
                        for j1 in state_family[parent]
                            push!(bond_set, (min(j, j1), max(j, j1)))
                        end
                    end
                end
            end
        end
    end

    return sort(collect(bond_set))
end

