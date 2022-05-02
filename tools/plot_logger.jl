using Plots
using CSV
using DataFrames

in_file = "logger_MD.txt"
println("What would you like to plot, buddy?")
println("timestep, Temp, PotentialEnergy, KineticEnergy, Momentum.")
println("First x.")
x = chomp(readline())
println("Now y.")
y = chomp(readline())

function plot_vel(x, y)
    
    #plot(x, y, label = "40_2000")
    plot!(x, y)
    xaxis!("$x")
    yaxis!("$y")
    title!("Logger")

end

data = CSV.read(in_file, DataFrame)
plot_vel(data.x, data.y)
savefig("Logger.png")