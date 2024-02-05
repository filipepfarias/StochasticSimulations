using GLMakie
using LinearAlgebra

function verify_collision(v,r,r1,r2,i,n)
    # i: i_pairs
    # n: n_particles
    ic = Vector{Bool}(undef,size(i,1));
    
    r_ = r[i[:,1],:] - r[i[:,2],:];
    d = norm.(eachrow(r_));

    ic[begin:n-1] = d[begin:n-1] .<= (r1+r2-1e-4r1);    
    ic[n:end] = d[n:end] .<= (2r2-1e-4r1);
    
    r_ = r[i[ic,1],:] - r[i[ic,2],:]
    v_ = v[i[ic,1],:] - v[i[ic,2],:]
    
    ic[ic] .*= eachrow(r_) .⋅ eachrow(v_) .< 0;
    
    return ic
end

function update_velocity!(v,r,i,ic,m)
    # i: i_pairs
    r1 = eachrow(@view r[i[ic,1],:]);
    r2 = eachrow(@view r[i[ic,2],:]);
    v1 = eachrow(@view v[i[ic,1],:]);
    v2 = eachrow(@view v[i[ic,2],:]);
    
    mi(i) = i == 1 ? m[1] : m[2];

    m1 = @. (2/(mi(i[ic,1]) + mi(i[ic,2]))) * mi(i[ic,2]);
    m2 = @. (2/(mi(i[ic,1]) + mi(i[ic,2]))) * mi(i[ic,1]);
    v_ = @. v1 - v2;
    r_ = @. r1 - r2;
    
    @. v1 = v1 - m1 * (v_) ⋅ r_ / norm(r_)^2 * r_
    @. v2 = v2 - m2 * (v_) ⋅ (r_) / norm(r_)^2 * (-r_)
end

# Settings
begin
    n_particles = 1024
    dt = 1e-6;
    box = (0.0,1.0,0.0,1.0);
    v₀ = 500; # m/s
    scale =  1e7 ;
    r₂ = .005#2.8e-10 * scale # water 
    r₁ = .05#27e-6 / 2 * 1e-3 *( 1e7 )# scalled polen
    m₂ = 1#(4π/3 * r₂^3) * (0.997); #(volume m³) * (density kg/m³)
    m₁ = 8#(4π/3 * r₁^3) * (1.435); #(volume m³) * (density kg/m³) http://www.jstor.org/stable/40586847
end

# Set initial conditions
begin
    r = zeros(n_particles,2);

    r[1,:] = [.5 .5]; # Center colloid
    
    v = [v₀*ones(n_particles) zeros(n_particles)];

    v[begin:(n_particles ÷ 2)] .*= -1;
    v[1,:] .*= 0.0; # Center colloid
    
    # Positioning smaller particles outside the bigger one
    ir = 2
    while ir <= n_particles
        ri = rand(1,2);
        drp = norm(ri - r[1:1,:]);
        if drp > (r₁ + r₂)
            r[ir,:] = ri;
            ir += 1;
        end
    end
    
    i_pairs = vcat([[x y] for x in 1:n_particles for y in x+1:n_particles]...,);
end;

# For plotting
begin
    ro = Observable(r);
    vo = Observable(v);
    vom = Observable([0.0]);
    rom = Observable([0.0]);
    gm = Observable([0.0]);
    
    f = Figure(size=(1000,500));
    ax = Axis(f[1:2,1], limits = box, aspect = 1);
    ax2 = Axis(f[1,2], limits = (0,5000,0,0.01),  aspect = 1);
    ax3 = Axis(f[2,2], limits = (-.025,.025,0,100),  aspect = 1);
    
    scatter!(ax,@lift($ro[1,1]),@lift($ro[1,2]),
    marker=Circle,
    markerspace=:data,
    markersize=2r₁)
    scatter!(ax,@lift($ro[2:end,1]),@lift($ro[2:end,2]),
    marker=Circle,
    markerspace=:data,
    markersize=2r₂)
    
    
    hist!(ax2,@lift( .√( sum(x -> x^2,$vo, dims=2) )[:] ); bins=50, normalization = :pdf)
    hist!(ax3, @lift($gm[:]) ; bins=10, normalization = :pdf)
    lines!(ax2,0:5000, v -> 2/v₀^2 * v * exp( - 2/v₀^2 * v^2 / 2), color=:black)
    vlines!(ax2,@lift(n_particles\sum(norm.(eachrow($vo)))), color=:red)
    vlines!(ax2,v₀, color=:black,)
    f
end


for fr in 1:10000
    ic = verify_collision(v,r,r₁,r₂,i_pairs,n_particles);
    
    update_velocity!(v,r,i_pairs,ic,[m₁,m₂]);
    
    inv_v = box[1]+r₁ .< r[1,1] .< box[2]-r₁;
    inv_v = vcat(inv_v, box[1]+r₂ .< r[2:end,1] .< box[2]-r₂)
    v[.!inv_v,1] .*= -1;
    
    inv_v = box[3]+r₁ .< r[1,2] .< box[4]-r₁;
    inv_v = vcat(inv_v, box[3]+r₂ .< r[2:end,2] .< box[4]-r₂)
    v[.!inv_v,2] .*= -1;
    
    r += v .* dt;
    
    
    ro[] = r;
    vo[] = v; 
    if (vom[] .!== v[1,1])[1]
        push!(gm[],r[1,1] - rom[][1]);  
        vom[] .= [v[1,1]];
        rom[] .= [r[1,1]];
        # display(gm[])
        # display(r[1,1] - rom[][end])
    end
    
    notify(ro);
    notify(vo);
    notify(gm);
    notify(rom)
    notify(vom);
    if fr%4 == 0
    sleep(1e-100)    
    end
end