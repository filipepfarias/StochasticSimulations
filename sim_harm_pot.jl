include("LangevinDynamics1D.jl");

begin
    fig = Figure(size=(500,500))
    ax1 = Axis(fig[1,1])
    Us1 = U.(xs,[1])
    Us2 = U.(xs,[2])

    lines!(ax1,xs,getindex.(Us1,1), color=Makie.wong_colors()[1],label=L"U(x),k=1")
    lines!(ax1,xs,getindex.(Us1,2),linestyle=:dash, color=Makie.wong_colors()[1],label=L"F(x), k=1")

    lines!(ax1,xs,getindex.(Us2,1), color=Makie.wong_colors()[2],label=L"U(x),k=2")
    lines!(ax1,xs,getindex.(Us2,2),linestyle=:dash, color=Makie.wong_colors()[2],label=L"F(x), k=2")
    axislegend()
end
fig

#this function returns the energy and force on a particle from a harmonic potential
function U(x,k=1,x0=0)
    #calculate the energy on force on the right hand side of the equal signs
    E = 0.5*k*(x-x0)^2
    F = -k*(x-x0)
    return E, F
end
my_k = 2
my_max_time = 100
initial_position = .1
initial_velocity = .5

my_gamma=10
my_kBT=0.25
my_dt=0.01

times, positions, velocities, total_energies = langevin_baoab(x -> U(x,my_k),
                                                              my_max_time, 
                                                              my_dt, 
                                                              my_gamma, 
                                                              my_kBT,
                                                              initial_position, 
                                                              initial_velocity)

begin
    fig = Figure(size=(500,500))
    ax1 = Axis(fig[1,1])
    ax2 = Axis(fig[2,1])

    lines!(ax1,times,velocities,label=L"v")
    lines!(ax1,times,positions,label=L"r")
    axislegend(ax1)
    lines!(ax2,times,total_energies,label=L"E")
    axislegend(ax2)
end
fig