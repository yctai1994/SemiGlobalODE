### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ a5881ab9-053c-47ac-b23c-da051ecc171d
import LinearAlgebra as la

# ╔═╡ 109c3c7c-ceb2-11ef-3483-89b230552d1f
const A0 = ComplexF64[
	complex(1.0, 2.0) complex(-5.0, 1.0) complex(4.0, 0.0);
	complex(6.0, 0.0) complex(16.0, 7.0) complex(-6.0, 8.0);
	complex(-4.0, 0.0) complex(2.0, 4.0) complex(-4.0, 1.0);
]

# ╔═╡ 0395bbcd-915a-48e9-8a63-64fb490bead7
let R = copy(A0), Q = similar(R)
	N = size(R, 2)
	n = Vector{ComplexF64}(undef, N)

	@inbounds for j in 1:N
		for i in 1:N
			Q[i,j] = complex(0.0, 0.0)
		end
		@inbounds Q[j,j] = complex(1.0, 0.0)
	end

	for k in 1:N-1
		xk = @view(R[k:N, k])
		nk = @view(n[k:N])

		sigma = la.norm(xk)
		phase = cis(angle(R[k,k]))
		let tmp = sigma * phase
			@inbounds for i in eachindex(nk)
				nk[i] = xk[i] + tmp * ifelse(i ≡ 1, 1.0, 0.0)
			end
		end

		nk_norm = la.norm(nk)
		@inbounds for i in eachindex(nk)
			nk[i] /= nk_norm
		end

		for j in k:N
			tmp = la.dot(nk, @view(R[k:N,j]))
			tmp *= 2.0
			for i in k:N
				R[i,j] -= tmp * n[i]
			end
		end

		for j in 1:N
			tmp = la.dot(nk, @view(Q[k:N,j]))
			tmp *= 2.0
			for i in k:N
				Q[i,j] -= tmp * n[i]
			end
		end
	end

	R, Q
end

# ╔═╡ a01dc808-324d-4c00-bcf5-153ad644db15
my_R, my_Q = let R = copy(A0), Q = similar(R)
	N = size(R, 2)

	@inbounds for j in 1:N
		for i in 1:N
			Q[i,j] = complex(0.0, 0.0)
		end
		@inbounds Q[j,j] = complex(1.0, 0.0)
	end

	for k in 1:N-1
		scale = maximum(abs, @view(R[k:N,k]))
		@inbounds for i in k:N
			R[i,k] /= scale
		end

		sigma = la.norm(@view(R[k:N,k]))
		phase = cis(angle(R[k,k]))
		beta = inv(sigma * (sigma + abs(R[k,k])))
		R[k,k] += sigma * phase

		for j in k+1:N
			tmp = la.dot(@view(R[k:N,k]), @view(R[k:N,j]))
			tmp *= beta
			for i in k:N
				R[i,j] -= tmp * R[i,k]
			end
		end

		for j in 1:N
			tmp = la.dot(@view(R[k:N,k]), @view(Q[k:N,j]))
			tmp *= beta
			for i in k:N
				Q[i,j] -= tmp * R[i,k]
			end
		end

		R[k,k] = -scale * sigma * phase
	end

	R, Q
end

# ╔═╡ 91a2d114-710c-446d-8e17-01916b936ba1
my_Q * my_Q'

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.2"
manifest_format = "2.0"
project_hash = "ac1187e548c6ab173ac57d4e72da1620216bce54"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"
"""

# ╔═╡ Cell order:
# ╠═a5881ab9-053c-47ac-b23c-da051ecc171d
# ╠═109c3c7c-ceb2-11ef-3483-89b230552d1f
# ╠═0395bbcd-915a-48e9-8a63-64fb490bead7
# ╠═a01dc808-324d-4c00-bcf5-153ad644db15
# ╠═91a2d114-710c-446d-8e17-01916b936ba1
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
