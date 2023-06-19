using Plots
using LinearAlgebra
using ProgressMeter
using Statistics
println("Enter input file (XYZ format).")
traj_file = chomp(readline())
println("Enter name of text out-file (e.g. acf.txt).")
out_file = chomp(readline())
println("Enter name of plot output file (eg:plotacf.png).")
file0 = chomp(readline())
println("How many frames would you like to analyze?")
n_frames = parse(Int64, chomp(readline()))
println("Enter chain length?")
len_chn = parse(Int64, chomp(readline()))
counter = 0
t=0
fh = 2 # File header.
N = 0 # Number of particles.
totalno_frames=0
# 3 
# 1005
#box = [xlen, ylen, zlen]
open(traj_file) do file
    for line in eachline(file)
        global counter += 1
        a = split(line)
        if counter == 1
            global N = parse(Int64, a[1])
            println("$(N/len_chn) chains in the system.")
        end
        break
    end
end
n_chain=Int(N/len_chn)
xframe=zeros(N)
yframe=zeros(N)
zframe=zeros(N)
x=[]
y=[]
z=[]



open(traj_file) do file
    global t = 0
    for line in eachline(file)
        t+=1
        m=mod(t,fh+N)
        if (m!=1 && m!=2 )
            a = split(line)
            ind=m-fh
            
            if(m==0)
                ind=N
            end
            xframe[ind]=parse(Float64,a[2])
            yframe[ind]=parse(Float64,a[3])
            zframe[ind]=parse(Float64,a[4])        
        elseif m==1   
            push!(x,copy(xframe))
            push!(y,copy(yframe))
            push!(z,copy(zframe))
            global totalno_frames+=1
        end

                   
    end
end
push!(x,copy(xframe))
push!(y,copy(yframe))
push!(z,copy(zframe))

xe2e=zeros(n_frames,n_chain)
ye2e=zeros(n_frames, n_chain)
ze2e=zeros(n_frames,n_chain)

start_frame=totalno_frames-n_frames 
for i in 1:n_chain
    b=(i-1)*len_chn + 1
    c=i*len_chn
    for frame in 1:n_frames
        xe2e[frame,i]=(x[start_frame+frame+1][c]-x[start_frame+frame+1][b])
        ye2e[frame,i]=(y[start_frame+frame+1][c]-y[start_frame+frame+1][b])
        ze2e[frame,i]=(z[start_frame+frame+1][c]-z[start_frame+frame+1][b])
        mag=(xe2e[frame,i]^2 + ye2e[frame,i]^2 + ze2e[frame,i]^2)^0.5 #magnitude
        xe2e[frame,i]/=mag #normalising
        ye2e[frame,i]/=mag
        ze2e[frame,i]/=mag
    end
end
l_acf=n_frames÷2
e2eacf=zeros(l_acf)
for i in 1:(n_frames-l_acf)
    l=0
    for j in i+1:(i+l_acf)
        l+=1
        for k in 1:n_chain
            e2eacf[l]+=(xe2e[j,k]*xe2e[i,k]) + (ye2e[j,k]*ye2e[i,k]) +(ze2e[j,k]*ze2e[i,k])
        end
    end         
end
e2eacf/=(n_chain*(n_frames-l_acf))
time = collect(1:l_acf)
time*=1000
plot(time,e2eacf)
xaxis!("t")
yaxis!("acf")
title!("ACF(End to End) vs t")
savefig(file0)
out_file = open(out_file, "w")
string_out = "t acf\n"
write(out_file, string_out)
for i in 1:length(e2eacf)
    global string_out = string(time[i]) * " " * string(e2eacf[i]) * "\n"
    write(out_file, string_out)
end