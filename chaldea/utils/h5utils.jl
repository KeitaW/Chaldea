function reconstructSparseMatrix(nrow, ncol, colptr, rowval, nzval)
    #=
    Reconstruct CSC matrix
    =#
    return(SparseMatrixCSC(nrow, ncol, colptr, rowval, nzval))
end

function initializeHDF5(sourcedir, filename)
    savefile = joinpath(sourcedir, filename)
    h5open(savefile, "w") do wf
        attrs(wf)["CreationDate"] = string(Dates.now())
        attrs(wf)["GitHEAD_ID"] = chomp(readall(`git rev-parse HEAD`))
    end
end

function saveSparseMatrixToHDF5(;h5file=nothing, groupName=nothing, cscmat=nothing, otherParams=nothing)
    #=
    groupName: Name of a matrix will be saved (should be full-path)
    otherParams: should be dict key: name of param value: value of param
    =#
    h5open(h5file, "r+") do wf
        local absolutePath = "/"
        local group
        for gn in unshift!(collect(splitdir(groupName)), "/")
            if !exists(wf[absolutePath], gn)
                group = g_create(wf, joinpath(absolutePath, gn))
                attrs(group)["CreationDate"] = string(Dates.now())
                attrs(group)["GitHEAD_ID"] = chomp(readall(`git rev-parse HEAD`))
            end
            absolutePath = joinpath(absolutePath, gn)
        end
        group = wf[absolutePath]
        group["colptr"] = cscmat.colptr
        group["rowval"] = cscmat.rowval
        group["nzval"] = cscmat.nzval
        group["nrow"] = size(cscmat)[1]
        group["ncol"] = size(cscmat)[2]
        if otherParams != nothing
            for (name, param) in otherParams
                attrs(group)[name] = param
            end
        end
    end
end

function readSparseMatrixFromHDF5(;h5file=nothing, groupName=nothing)
    #=
    groupName: usually, it should be the name of a matrix
    =#
    h5open(h5file, "r") do rf
        g = rf[groupName]
        return SparseMatrixCSC(read(g, "nrow"), read(g, "ncol"), read(g, "colptr"),
            read(g, "rowval"), read(g, "nzval"))
    end
end
