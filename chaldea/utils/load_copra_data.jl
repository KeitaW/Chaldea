function load_copra_data(rf; drop=5)
    lines = open(rf, "r") do fp
    readlines(fp)
    end
    data = []
    for line in lines
        tmp = [parse(Int, l) for l in split(line, ' ')[1:(end-1)]]
        if length(tmp) >= drop
            push!(data, tmp)
        end
    end
    return(data)
end
