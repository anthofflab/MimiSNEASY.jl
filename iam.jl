module Iam

abstract ComponentParameters
abstract ComponentVariables

type ComponentInstance
	p::ComponentParameters
	v::ComponentVariables
end

import Base.show
show(io::IO, a::ComponentInstance) = print(io, "A component")

function run(nsteps::Int, components::Vector{ComponentInstance})
	ncomponents = size(components,1)
	for i=1:ncomponents
		init(components[i].ComponentParameters,components[i].ComponentVariables)
	end

	for t=1:nsteps
		for i=1:ncomponents
			timestep(components[i].ComponentParameters,components[i].ComponentVariables,t)
		end
	end
end

end