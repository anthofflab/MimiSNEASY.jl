using Mimi
include("sneasy.jl")

function doruns(model)
    for i=1:30000
        run(model)
    end
end

m = get_model()
run(m)

@time doruns(m)
