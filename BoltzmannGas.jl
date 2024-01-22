## Boltzmann Distribution by Gas Simulation
## Source: https://github.com/lukepolson/youtube_channel/blob/main/Python%20GPU/multibody_boltzmann.ipynb

using GLMakie
using LinearAlgebra

n_particles = 2048;

box = (0,1,0,1);

r = rand(n_particles,2);

i_r_r = r[:,1] .>  .5;
i_r_l = r[:,1] .<= .5;

i_r = 1:n_particles;

# Setting velocities
v = zero(r);
v[i_r_r,1] .= -500;
v[i_r_l,1] .=  500;

# Calculating distances
i_pairs = vcat([[x y] for x in i_r for y in x+1:n_particles]...,)

# Verifying collision
radius = 0.0005;

# Calculating new velocities
function update_v(v,r,i)
    i1,i2 = i[:,1],i[:,2];
    vnew1 = sum((v[i1,:] - v[i2,:]) .* (r[i1,:] - r[i2,:]), dims=2); 
    vnew1 ./= sum(x -> x^2,r[i1,:] - r[i2,:],dims=2); 
    vnew1 = v[i1,:] - vnew1 .* (r[i1,:] - r[i2,:])

    vnew2 = sum((v[i1,:] - v[i2,:]) .* (r[i1,:] - r[i2,:]), dims=2) 
    vnew2 ./= sum(x -> x^2,r[i2,:] - r[i1,:],dims=2);
    vnew2 = v[i2,:] - vnew2 .* (r[i2,:] - r[i1,:])

    return vnew1, vnew2
end


function gas_motion(r::Observable, v, i_pairs, ts, dt, d_cutoff)

    f = Figure();
    ax = Axis(f[1,1], limits = (0,1,0,1));
    ax2 = Axis(f[1,2], limits = ((0,1500),nothing));

    vo = Observable(v);

    i_r_r = r[][:,1] .>  .5;
    i_r_l = r[][:,1] .<= .5;

    scatter!(ax,@lift($r[i_r_r,1]),(@lift $r[i_r_r,2]),markersize=2)
    scatter!(ax,(@lift $r[i_r_l,1]),(@lift $r[i_r_l,2]),markersize=2)
    hist!(ax2,@lift( .√( sum(x -> x^2,$vo, dims=2) )[:] ); bins=50)
    
    
    record(f, "gas.mp4", 1:ts, framerate=60) do frame

        Δx = r[][i_pairs[:,1],1] - r[][i_pairs[:,2],1];
        Δy = r[][i_pairs[:,1],2] - r[][i_pairs[:,2],2];
        d = .√(Δx.^2 + Δy.^2)
        i_collide = i_pairs[d .< d_cutoff,:]; 

        v[i_collide[:,1],:],v[i_collide[:,2],:] = copy.(update_v(v,r[],i_collide)); 

        inv_v = box[1] .< r[][:,1] .< box[2];
        v[.!inv_v,1] .*= -1;
    
        inv_v = box[3] .< r[][:,2] .< box[4];
        v[.!inv_v,2] .*= -1;

        r[] += v .* dt;
        vo[] = v;
        
        notify(r);
        notify(vo);
    end
    
end

gas_motion(Observable(r),v,i_pairs,1000,8e-6,2radius)