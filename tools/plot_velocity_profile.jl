using Plots
using CSV
using DataFrames

function plot_vel(x, y)
    
    #plot(x, y, label = "40_2000")
    plot!(x, y)
    xaxis!("v_x")
    yaxis!("z")
    title!("Velocity profile")

end

println("Enter name of input file.")
file1 = chomp(readline())
println("Enter name of output file.")
file0 = chomp(readline())

data = CSV.read(file1, DataFrame)
plot_vel(data.z, data.v_x)
savefig(file0)