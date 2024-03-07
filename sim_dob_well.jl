using ForwardDiff

include("LangevinDynamics1D.jl");

#this function returns the energy and force on a particle from a double well
function U(x,k=1,a=3)
    #calculate the energy on force on the right hand side of the equal signs
    E = 0.25*k*(x-a)^2 * (x+a)^2
    # E(x) = 0.25*k*(x-a)^2 * (x+a)^2
    F = -k*x*(x-a)^1*(x+a)^1
    # F(x) = ForwardDiff.derivative(E,x)
    return E, F
    # return E.(x), F.(x)
end
my_k = .1
my_a = 4
my_gamma=0.1
my_kBT=1.0
my_dt=0.05

begin
    xs = -6:.2:6
    fig = Figure(size=(500,500))
    ax1 = Axis(fig[1,1])
    Us1 = U.(xs,[my_k],[my_a])

    lines!(ax1,xs,getindex.(Us1,1), color=Makie.wong_colors()[1],label=L"U(x),k=1")
    lines!(ax1,xs,getindex.(Us1,2),linestyle=:dash, color=Makie.wong_colors()[1],label=L"F(x), k=1")
    hlines!(ax1,my_kBT,color=:black)
    limits!((-6,6),(-10,10))
    axislegend()
end
fig


my_max_time = 10_000
initial_position = .1
initial_velocity = .5

my_initial_position = my_a
my_initial_velocity = 1

times, positions, velocities, total_energies = langevin_baoab(x -> U(x,my_k,my_a),
                                                              my_max_time, 
                                                              my_dt, 
                                                              my_gamma, 
                                                              my_kBT,
                                                              my_initial_position, 
                                                              my_initial_velocity)

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