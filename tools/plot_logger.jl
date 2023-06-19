using Plots
using CSV
using DataFrames

function plot_log(x, y)
    
    #plot(x, y, label = "40_2000")
    plot!(x, y)
    xaxis!("x")
    yaxis!("y")
    title!("Logger")

end

file1 = "log.txt"
file0 = "out.png"

#println("Enter name of input file.")
#file1 = chomp(readline())
#println("Enter name of output file.")
#file0 = chomp(readline())
#in_file = "logger_MD.txt"
#println("What would you like to plot, buddy?")
#println("t (timestep), T (temperature), V (potential energy, K (kinetic energy), p (momentum).")
#println("First x.")
#x_data = chomp(readline())
#println("Now y.")
#y_data = chomp(readline())

df = CSV.read(file1, DataFrame, header = ["t", "T", "V", "K", "p"], skipto = 2)
plot_log(df.t, df.p)
savefig(file0)