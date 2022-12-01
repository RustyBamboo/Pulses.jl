module magnus
commu(x,y) = x * y - y * x

c_1 = 0.5 - sqrt(3) / 6
c_2 = 0.5 + sqrt(3) / 6


# Need to define a function to compute the time-dependent parameter at a time t
function solve(H::Array{Matrix}, pulses::Array{Real}, time::Array{Real})

end

end

a = [1 1; 0 1]
b = [1 0; 0 1]

using .magnus

print(magnus.commu(a, b))