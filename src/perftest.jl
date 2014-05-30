using IAMF
include("sneasy.jl")

function doruns(model)
    for i=1:30000
        run(566,model)
    end
end

m = getsneasy()
run(566,m)

@time doruns(m)
