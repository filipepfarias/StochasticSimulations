using CUDA

N = 2^20
x = fill(1.0f0, N)  # a vector filled with 1.0 (Float32)
y = fill(2.0f0, N)  # a vector filled with 2.0

y .+= x             # increment each element of y with the corresponding element of x

all(y .== 3.0f0)

function sequential_add!(y, x)
    for i in eachindex(y, x)
        @inbounds y[i] += x[i]
    end
    return nothing
end

fill!(y, 2)
sequential_add!(y, x)
all(y .== 3.0f0)

function parallel_add!(y, x)
    Threads.@threads for i in eachindex(y, x)
        @inbounds y[i] += x[i]
    end
    return nothing
end

fill!(y, 2)
parallel_add!(y, x)
all(y .== 3.0f0)

@time sequential_add!(y, x)
@time parallel_add!(y, x)

x_d = CUDA.fill(1.0f0, N)  # a vector stored on the GPU filled with 1.0 (Float32)
y_d = CUDA.fill(2.0f0, N)  # a vector stored on the GPU filled with 2.0

y_d .+= x_d
all(Array(y_d) .== 3.0f0)

function add_broadcast!(y, x)
    CUDA.@sync y .+= x
    return
end

@time add_broadcast!(y_d, x_d)

function gpu_add1!(y, x)
    for i = 1:length(y)
        @inbounds y[i] += x[i]
    end
    return nothing
end

fill!(y_d, 2)
@cuda gpu_add1!(y_d, x_d)
all(Array(y_d) .== 3.0f0)

function bench_gpu1!(y, x)
    CUDA.@sync begin
        @cuda gpu_add1!(y, x)
    end
end

@time bench_gpu1!(y_d, x_d)

function gpu_add2!(y, x)
    index = threadIdx().x    # this example only requires linear indexing, so just use `x`
    stride = blockDim().x
    for i = index:stride:length(y)
        @inbounds y[i] += x[i]
    end
    return nothing
end

fill!(y_d, 2)
@cuda threads=256 gpu_add2!(y_d, x_d)
all(Array(y_d) .== 3.0f0)

function bench_gpu2!(y, x)
    CUDA.@sync begin
        @cuda threads=256 gpu_add2!(y, x)
    end
end

@time bench_gpu2!(y_d, x_d)

function gpu_add3!(y, x)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride = gridDim().x * blockDim().x
    for i = index:stride:length(y)
        @inbounds y[i] += x[i]
    end
    return
end

numblocks = ceil(Int, N/256)

fill!(y_d, 2)
@cuda threads=256 blocks=numblocks gpu_add3!(y_d, x_d)
all(Array(y_d) .== 3.0f0)

function bench_gpu3!(y, x)
    numblocks = ceil(Int, length(y)/256)
    CUDA.@sync begin
        @cuda threads=256 blocks=numblocks gpu_add3!(y, x)
    end
end

@time bench_gpu3!(y_d, x_d)

