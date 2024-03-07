using CairoMakie

include("LangevinDynamics2D.jl");

#this function returns the energy and force on a particle from a harmonic potential
function U(x,k=1,x0=[0 0])
    #calculate the energy on force on the right hand side of the equal signs
    E = 0.5 .* k .* ((x-x0)*(x-x0)')[1] 
    F = -k .* (x-x0)
    return E, F
end

begin
    xs = -5:.1:5
    ys = -5:.1:5
    fig = Figure(size=(500,500))
    ax1 = Axis(fig[1,1])

    heatmap!(ax1,xs,ys,(x,y) -> U([x y])[1])
    Colorbar(fig[1,2])
    # lines!(ax1,xs,getindex.(Us1,2),linestyle=:dash, color=Makie.wong_colors()[1],label=L"F(x), k=1")

    # axislegend()
end
fig

my_k = .00

initial_position = [.5 .5]
initial_velocity = [.0 .0]

my_gamma= #100000
my_kBT=10000
my_dt=0.0001
my_max_time = 100 * my_dt

times, positions, velocities, total_energies = langevin_baoab(x -> U(x,my_k),
                                                              my_max_time, 
                                                              my_dt, 
                                                              my_gamma, 
                                                              my_kBT,
                                                              initial_position, 
                                                              initial_velocity)

begin
    fig = Figure(size=(500,500))
    ax1 = Axis(fig[1,1], limits=(0,1,0,1))
    ax2 = Axis(fig[2,1])

    # lines!(ax1,velocities[:,1],velocities[:,2],label=L"v")
    lines!(ax1,positions[:,1],positions[:,2],label=L"r")
    axislegend(ax1)
    # lines!(ax2,times,total_energies,label=L"E")
    # axislegend(ax2)
end
fig

begin
    set_theme!(theme_latexfonts())

    fig = Figure(size=(500,500))
    ax1 = Axis(fig[1,1], limits=(0,1,0,1), aspect = 1)
    # ax2 = Axis(fig[2,1])

    # lines!(ax1,velocities[:,1],velocities[:,2],label=L"v")
    lines!(ax1,positions[:,1],positions[:,2],color=:black)
    scatter!(ax1,positions[end,1],positions[end,2],
    marker=Circle,
    markerspace=:data,
    color=:gray,
    strokewidth=1,
    markersize=2 * .025)
    # axislegend(ax1)
    # lines!(ax2,times,total_energies,label=L"E")
    # axislegend(ax2)
end
fig
