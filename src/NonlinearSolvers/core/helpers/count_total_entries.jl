# Fallback method (throws error for unsupported types)
count_total_entries(var) = error("Unsupported structure in solution_variables of type $(typeof(var))")

# For scalars
count_total_entries(var::Number) = 1

# For flat arrays of numbers
count_total_entries(var::AbstractArray{<:Number}) = length(var)

# For empty arrays
count_total_entries(var::AbstractArray) = 0

# For nested arrays (e.g., tuples of arrays, arrays of arrays)
count_total_entries(var::AbstractVector) = sum(count_total_entries.(var))
count_total_entries(var::Tuple) = sum(count_total_entries.(var))

# Need to revisit this one
if isdefined(Main, :ContinuationVariable)
    @eval count_total_entries(var::ContinuationVariable) = 1
end
