using MakiePublication
using CairoMakie
using LaTeXStrings
using CSV
using DataFrames

function plot_rdf()
    
    #plot(x, y, label = "40_2000")
    fig = Figure(figure_padding=(0, 6, 0, 4))
    ax = Axis(fig, xlabel =L"r", ylabel="g(r)")
    file1 = "rdf.txt"
    data = CSV.read(file1, DataFrame)
    lines!(data.r, data.g, label=L"RDF")
    #xaxis!("r")
    #yaxis!("g(r)")
    #title!("Radial distribution function")

    fig[1,1] = ax
    return fig
end

#println("Enter name of input file.")
#file1 = chomp(readline())
#println("Enter name of output file.")
#file0 = chomp(readline())
file0 = "rdf.pdf"


fig = with_theme(plot_rdf, theme_acs())
savefig(file0, fig)
