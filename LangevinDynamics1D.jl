using CairoMakie
using LinearAlgebra

# Adapted from https://hockygroup.hosting.nyu.edu/exercise/langevin-dynamics.html

#this function returns the energy and force on a particle from a harmonic potential
function U(x,k=1,x0=0)
    #calculate the energy on force on the right hand side of the equal signs
    E = 0.5*k*(x-x0)^2
    F = -k*(x-x0)
    return E, F
end

#this function will plot the energy and force
#it is very general since it uses a special python trick of taking arbitrary named arguments (**kwargs) 
#and passes them on to a specified input function
# function plot_energy_force(fun::Function, xmin=-3,xmax=3,spacing=0.1)
#     xs = range(xmin,xmax+spacing,spacing)
#     energies, forces = fun(xs)
#     label = "U(x)"
#     for arg in kwargs:
#         label=label+", %s=%s"%(arg,str(kwargs[arg]))
#     p = plt.plot(xs,energies,label=label)
#     plt.plot(xs,forces,label="",color=p[0].get_color(),linestyle=:dash)
#     plt.legend(loc=0)
# end
    
#we can plot the energy (solid) and forces (dashed) to see if it looks right
# plot_energy_force(harmonic_oscillator_energy_force,k=1)
# plot_energy_force(harmonic_oscillator_energy_force,k=2)

function update_position(x,v,dt)
   return x + v .* dt/2
end

function update_velocity(v,F,dt)
    return v + F .* dt/2
end

# Ornstein-Uhlenbeck process
function update_velocity_OU(v,gamma,kBT,dt)
    R = randn()
    c1 = exp(-gamma*dt)
    c2 = sqrt(1-c1*c1)*sqrt(kBT)
    return c1*v + R*c2
end

function langevin_baoab(potential, max_time, dt, gamma, kBT, initial_position, initial_velocity,save_frequency=3)
    x = initial_position
    v = initial_velocity
    t = 0
    step_number = 0
    positions = Vector{Float64}()
    velocities = Vector{Float64}()
    total_energies = Vector{Float64}()
    save_times = Vector{Float64}()
    
    while t<max_time
        # B
        Ux, F = potential(x)
        v = update_velocity(v,F,dt)

        #A
        x = update_position(x,v,dt)

        #O
        v = update_velocity_OU(v,gamma,kBT,dt)
        
        #A
        x = update_position(x,v,dt)
        
        # B
        Ux, F = potential(x)
        v = update_velocity(v,F,dt)
        
        if ((step_number%save_frequency == 0) || (step_number>0))
            e_total = .5*v*v + Ux

            push!(positions,x)
            push!(velocities,v)
            push!(total_energies,e_total)
            push!(save_times,t)
        end
        t = t+dt
        step_number = step_number + 1
    end
    return save_times, positions, velocities, total_energies
end