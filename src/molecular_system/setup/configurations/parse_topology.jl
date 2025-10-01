"""
    parse_topology(sequence::String) -> TopologyStruct

Parses a monomer sequence string into a hierarchical topology tree. Handles:

- Linear sequences (e.g., `"ABC"`)
- Branched sequences (e.g., `"A(BC)"`)
- Repetition via multiplier notation (e.g., `"A(BC)2"` → `"ABCBC"`)

Returns a `TopologyStruct` with parent-child relationships, level indexing,
and sequential bond information for use in chain modeling.
"""



function parse_topology(sequence::String)
    # Initialize data containers
    indices = Int[]
    levels = Vector{Vector{Int}}()
    children = Vector{Vector{Int}}()
    parents = Vector{Vector{Int}}()
    monomer_list  = String[]
    sequ_bonds = Tuple{Int, Int}[]  # <-- NEW

    current_index = Ref(0)
    last_linear_index = Ref(0)

    function recurse(subseq::String, current_level::Int, parent_index::Int)
        char = subseq[1]

        if isletter(char)
            current_index[] += 1
            push!(monomer_list , string(char))

            if current_level > length(levels)
                push!(levels, Int[])
            end

            push!(indices, current_index[])
            push!(levels[current_level], current_index[])

            # Add children and parents structure
            if parent_index == 0
                push!(children, [])
                push!(parents, [])
            else
                if parent_index > length(children)
                    push!(children, [current_index[]])
                else
                    push!(children[parent_index], current_index[])
                end

                if current_index[] > length(parents)
                    push!(parents, [parent_index])
                else
                    push!(parents[current_index[]], parent_index)
                end
            end

            # Sequential bond tracking
            if last_linear_index[] != 0
                push!(sequ_bonds, (last_linear_index[], current_index[]))
            end
            last_linear_index[] = current_index[]

            # Recurse forward
            if length(subseq) > 1
                recurse(subseq[2:end], current_level + 1, current_index[])
            else
                push!(children, [])
            end

        elseif char == '('
            end_index = find_matching_parenthesis(subseq)
            group = subseq[2:end_index-1]
            remaining = subseq[end_index+1:end]
            multiplier = 1

            if !isempty(remaining) && isdigit(remaining[1])
                digit_chars = Char[]
                i = 1
                while i <= lastindex(remaining) && isdigit(remaining[i])
                    push!(digit_chars, remaining[i])
                    i += 1
                end
                multiplier = parse(Int, String(digit_chars))
                remaining = remaining[i:end]
            end

            for _ in 1:multiplier
                save_parent = parent_index
                save_level = current_level
                recurse(group, current_level, parent_index)
                parent_index = save_parent
                current_level = save_level
            end

            if !isempty(remaining)
                recurse(remaining, current_level, parent_index)
            end

        else
            error("Not a recognized character: $char")
        end
    end

    function find_matching_parenthesis(seq::String)
        balance = 1
        for j in firstindex(seq)+1:lastindex(seq)
            if seq[j] == '('
                balance += 1
            elseif seq[j] == ')'
                balance -= 1
            end
            if balance == 0
                return j
            end
        end
        error("Unmatched parenthesis in: $seq")
    end

    recurse(sequence, 1, 0)

    monomer_list  = join(monomer_list , "")

    return TopologyStruct(monomer_list , indices, children, parents, levels, sequ_bonds)
end

