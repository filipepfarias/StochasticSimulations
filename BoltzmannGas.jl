## Boltzmann Distribution by Gas Simulation
## Source: https://github.com/lukepolson/youtube_channel/blob/main/Python%20GPU/multibody_boltzmann.ipynb

using GLMakie
using LinearAlgebra

n_particles = 2048;

box = (0,1,0,1);

r = rand(n_particles,2);
# r = [.2 .5; .8-.006 .5];

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
radius = 0.005;

# Calculating new velocities
function collision(v,r,i)
    i1,i2 = i[:,1],i[:,2];
    vnew1 = sum((v[i1,:] - v[i2,:]) .* (r[i1,:] - r[i2,:]), dims=2); 
    vnew1 ./= sum(x -> x^2,r[i1,:] - r[i2,:],dims=2); 
    vnew1 = v[i1,:] - vnew1 .* (r[i1,:] - r[i2,:])

    vnew2 = sum((v[i1,:] - v[i2,:]) .* (r[i1,:] - r[i2,:]), dims=2) 
    vnew2 ./= sum(x -> x^2,r[i2,:] - r[i1,:],dims=2);
    vnew2 = v[i2,:] - vnew2 .* (r[i2,:] - r[i1,:])

    return vnew1, vnew2
end

function update_particles(r,v,i_pairs,dt,d_cutoff)
    Δx = r[i_pairs[:,1],1] - r[i_pairs[:,2],1];
    Δy = r[i_pairs[:,1],2] - r[i_pairs[:,2],2];
    d = .√(Δx.^2 + Δy.^2)
    i_c = i_pairs[d .< d_cutoff,:]; 

    Δr = r[i_c[:,1],:] - r[i_c[:,2],:];
    Δv = v[i_c[:,1],:] - v[i_c[:,2],:];
    breakin = eachrow(Δr) .⋅ eachrow(Δv) .<= 0;
    i_c = i_c[breakin,:];
    
    v[i_c[:,1],:],v[i_c[:,2],:] = collision(v,r,i_c); 
    
    inv_v = box[1] .< r[:,1] .< box[2];
    v[.!inv_v,1] .*= -1;
    
    inv_v = box[3] .< r[:,2] .< box[4];
    v[.!inv_v,2] .*= -1;
    
    r += v .* dt;
    return r,v
end


function gas_motion(r, v, i_pairs, ts, dt, d_cutoff)

    f = Figure(size=(1000,500));
    ax = Axis(f[1,1], limits = box, aspect = 1);
    ax2 = Axis(f[1,2], limits = (0,1500,0,0.01),  aspect = 1);

    ro = Observable(r);
    vo = Observable(v);

    i_r_r = ro[][:,1] .>  .5;
    i_r_l = ro[][:,1] .<= .5;

    scatter!(ax,@lift($ro[i_r_r,1]),(@lift $ro[i_r_r,2]),
                    marker=Circle,
                    markerspace=:data,
                    markersize=2radius)
    scatter!(ax,(@lift $ro[i_r_l,1]),(@lift $ro[i_r_l,2]),
                    marker=Circle,
                    markerspace=:data,
                    markersize=2radius)
    hist!(ax2,@lift( .√( sum(x -> x^2,$vo, dims=2) )[:] ); bins=500, normalization = :pdf)
    lines!(ax2,0:1000, v -> 2/500^2 * v * exp( - 2/500^2 * v^2 / 2))
    
    
    record(f, "gas.mp4", 1:ts, framerate=60) do frame

        r,v = update_particles(r,v,i_pairs,dt,d_cutoff)
        ro[] = r;
        vo[] = v;
        
        notify(ro);
        notify(vo);
    end
    
end

gas_motion(r,v,i_pairs,1000,8e-6,2radius)