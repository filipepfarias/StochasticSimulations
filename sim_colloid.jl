using CairoMakie

include("BrownianMotion.jl")

# Settings
begin
    n_particles = 1024
    dt = 1e-6;
    box = (0.0,1.0,0.0,1.0);
    v₀ = 500; # m/s
    scale =  1e7 ;
    r₂ = .005#2.8e-10 * scale # water 
    r₁ = .025#27e-6 / 2 * 1e-3 *( 1e7 )# scalled polen
    m₂ = 1#(4π/3 * r₂^3) * (0.997); #(volume m³) * (density kg/m³)
    m₁ = 5#(4π/3 * r₁^3) * (1.435); #(volume m³) * (density kg/m³) http://www.jstor.org/stable/40586847
end

# For plotting
begin
    ro = Observable(r);
    vo = Observable(v);
    # vom = Observable([0.0]);
    rp = Observable(Point2f[]);
    gm = Observable([0.0]);
    
    set_theme!(theme_latexfonts())
    f = Figure(size=(500,500));
    ax = Axis(f[1:2,1], limits = box, aspect = 1);
    # ax2 = Axis(f[1,2], limits = (0,5000,0,0.01),  aspect = 1);
    # ax3 = Axis(f[2,2], limits = (-.025,.025,0,100),  aspect = 1);
    
    scatter!(ax,@lift($ro[1,1]),@lift($ro[1,2]),
    marker=Circle,
    markerspace=:data,
    color=:gray,
    strokewidth=1,
    markersize=2r₁)
    scatter!(ax,@lift($ro[2:end,1]),@lift($ro[2:end,2]),
    marker=Circle,
    markerspace=:data,
    color=:black,
    markersize=2r₂)
    lines!(ax,rp, color=:black)

    
    
    # hist!(ax2,@lift( .√( sum(x -> x^2,$vo, dims=2) )[:] ); bins=50, normalization = :pdf)
    # hist!(ax3, @lift($gm[:]) ; bins=1000, normalization = :pdf)
    # lines!(ax2,0:5000, v -> 2/v₀^2 * v * exp( - 2/v₀^2 * v^2 / 2), color=:black)
    # vlines!(ax2,@lift(n_particles\sum(norm.(eachrow($vo)))), color=:red)
    # vlines!(ax2,v₀, color=:black,)
    f
end


rom = 0.5;
# file_id = 0;
# save("outputs/brownian_f"*string(file_id)*".pdf",f);
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
    push!(rp[],Point2f(r[1,:]...))
    # if (vom[] .!== v[1,1])[1]
    #     push!(gm[],r[1,1] - rom[][1]);  
    #     vom[] .= [v[1,1]];
    #     rom[] .= [r[1,1]];
    #     # display(gm[])
    #     # display(r[1,1] - rom[][end])
    # end
    
    notify(ro);
    notify(vo);
    notify(rp)
    # notify(vom);
    if fr%200 == 0
    push!(gm[],gm[][1]-rom)
    # display(-rom + gm[][1])
    rom = r[1,1];
    notify(gm);
    # sleep(1e-100)  
    file_id += 1;  
    # save("outputs/brownian_f"*string(file_id)*".pdf",f);
    end
end