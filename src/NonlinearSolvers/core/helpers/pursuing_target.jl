function pursuing_target(x, start, target; inclusive=true, check_direction=true)
    if target === nothing
        return true
    end

    in_interval(x, a, b; inclusive=true) = inclusive ?
                (x >= min(a, b) && x <= max(a, b)) :
                (x >  min(a, b) && x <  max(a, b))

    # Check if value is in the interval between start and target
    in_range = in_interval(x, start, target; inclusive=inclusive)

    # Check direction if requested
    if check_direction
        dir_expected = sign(target - start)
        dir_actual   = sign(x - start)
        correct_direction = dir_actual == dir_expected
        return in_range && correct_direction
    else
        return in_range
    end
end
