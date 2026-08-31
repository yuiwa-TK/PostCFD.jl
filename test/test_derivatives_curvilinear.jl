@testset "Compact derivative on general curvilinear coordinates" begin
    # Curved mapping from computational coordinates:
    #   x = ξ + aη² + 10
    #   y = η + bζ² - 4
    #   z = ζ + 2
    # Its inverse derivatives are stored as metrics[physical, computational].
    for T in (Float64, Float32)
        dims = (9, 8, 7)
        a, b = T(1 // 20), T(1 // 30)
        f = Array{T}(undef, dims)
        metrics = zeros(T, dims..., 3, 3)

        gradient = T[5 // 4, -2, 3 // 4]
        for l in axes(f, 3), k in axes(f, 2), j in axes(f, 1)
            ξ, η, ζ = T(j - 1), T(k - 1), T(l - 1)
            metrics[j, k, l, 1, 1] = one(T)       # ξ_x
            metrics[j, k, l, 2, 1] = -2a * η     # ξ_y
            metrics[j, k, l, 2, 2] = one(T)       # η_y
            metrics[j, k, l, 3, 1] = 4a * b * η * ζ # ξ_z
            metrics[j, k, l, 3, 2] = -2b * ζ     # η_z
            metrics[j, k, l, 3, 3] = one(T)       # ζ_z

            x = ξ + a * η^2 + T(10)
            y = η + b * ζ^2 - T(4)
            z = ζ + T(2)
            f[j, k, l] = gradient[1] * x + gradient[2] * y +
                         gradient[3] * z + T(7)
        end

        original_f = copy(f)
        df = similar(f)
        for direction in 1:3
            result = derivative_curvilinear_inplace_compact!(
                df, f, metrics, direction
            )

            @test result === df
            @test df ≈ fill(gradient[direction], dims) rtol=100eps(T) atol=100eps(T)
            @test f == original_f
        end
    end
end
