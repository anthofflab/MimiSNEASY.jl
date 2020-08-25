using Mimi
using MimiSNEASY

function doruns(model)
    for i=1:30000
        run(model)
    end
end

m = MimiSNEASY.get_model()
run(m)

@time doruns(m)
