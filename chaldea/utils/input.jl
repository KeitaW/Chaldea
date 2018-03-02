# Get input from stdin. The script is written because of the lack of such a function in julia.
# Reference: http://stackoverflow.com/questions/17479782/julia-request-user-input-from-script
# At the mean time, you can parse the input as follows: x = parse(Int, input())
function input(prompt::AbstractString="")
    print(prompt)
    return chomp(readline())
end
