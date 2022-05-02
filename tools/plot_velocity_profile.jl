using Plots
using CSV
using DataFrames

println("Enter name of input file.")
file1 = chomp(readline())
println("Enter name of output file.")
file0 = chomp(readline())

function plot_vel(x, y)
    
    #plot(x, y, label = "40_2000")
    plot!(x, y)
    xaxis!("z")
    yaxis!("v_x")
    title!("Velocity profile")

end

data = CSV.read(file1, DataFrame)
plot_vel(data.z, data.v_x)
savefig(file0)