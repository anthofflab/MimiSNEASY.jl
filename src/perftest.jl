using IAMF
include("sneasy.jl")

function doruns(model)
    for i=1:30000
        run(model)
    end
end

m = getsneasy()
run(m)

@time doruns(m)
