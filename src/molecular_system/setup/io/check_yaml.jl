"""
    check_yaml(input_file::String)

Reads a YAML file and checks for duplicate monomer keys under the `monomers:` section.
Throws an error if duplicates are found.
"""
function check_yaml(input_file :: String)
    yaml_text = read(input_file, String)
    dupes = check_duplicate_keys_in_section(yaml_text, "monomers")
    if !isempty(dupes)
        error("Duplicate monomer keys in 'monomers' section: " * join(dupes, ", "))
    end
end

"""
    check_duplicate_keys_in_section(yaml_text::String, section::String) -> Vector{String}

Parses a raw YAML string and returns a list of duplicate keys in the specified section.
Handles indentation to ensure keys are correctly grouped under the section.
"""
function check_duplicate_keys_in_section(yaml_text::String, section::String)
    lines = split(yaml_text, '\n')
    in_section = false
    key_indent = nothing
    seen_keys = Set{String}()
    duplicates = String[]

    for line in lines
        if occursin(r"^\s*$", line) || startswith(strip(line), "#")
            continue  # skip empty lines or comments
        end

        if !in_section
            if strip(line) == section * ":"
                in_section = true
            end
            continue
        end

        if in_section
            # detect indentation and key
            m = match(r"^(\s*)([^:\s]+):", line)
            if m !== nothing
                indent, key = m.captures
                if key_indent === nothing
                    key_indent = indent
                elseif length(indent) < length(key_indent)
                    break  # out of section
                elseif length(indent) > length(key_indent)
                    continue  # deeper indentation; value lines
                end

                if key in seen_keys
                    push!(duplicates, key)
                else
                    push!(seen_keys, key)
                end
            else
                continue  # not a key line
            end
        end
    end

    return duplicates
end
