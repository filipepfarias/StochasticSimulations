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

# Set initial conditions
function init_brownian(n_particles::Int64,v₀::Float64)
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

    return r,v,i_pairs
end;

function evolve_brownian!(r,v,dt,box,part_params)
    m₁,r₁,m₂,r₂ = part_params;

    update_velocity!(v,r,i_pairs,ic,[m₁,m₂]);
    
    inv_v = box[1]+r₁ .< r[1,1] .< box[2]-r₁;
    inv_v = vcat(inv_v, box[1]+r₂ .< r[2:end,1] .< box[2]-r₂)
    v[.!inv_v,1] .*= -1;
    
    inv_v = box[3]+r₁ .< r[1,2] .< box[4]-r₁;
    inv_v = vcat(inv_v, box[3]+r₂ .< r[2:end,2] .< box[4]-r₂)
    v[.!inv_v,2] .*= -1;
    
    r += v .* dt;
end