### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ a5881ab9-053c-47ac-b23c-da051ecc171d
import LinearAlgebra as la

# ╔═╡ 109c3c7c-ceb2-11ef-3483-89b230552d1f
const A0 = [ 12.0 -51.0 4.0; 6.0 167.0 -68.0; -4.0 24.0 -41.0; ]

# ╔═╡ 3515d200-b8ef-4149-8b59-d2dfffeb6ac9
let A = copy(A0)
	N = size(A, 2)
	n = Vector{Float64}(undef, N)

	for k in 1:N
		xk = @view(A[k:N, k])
		nk = @view(n[k:N])
	
		xk_norm = la.norm(xk)
		for i in eachindex(nk)
			@inbounds nk[i] = xk[i] - xk_norm * ifelse(i ≡ 1, 1.0, 0.0)
		end
	
		nk_norm = la.norm(nk)
		for i in eachindex(nk)
			nk[i] /= nk_norm
		end
	
		for j in k:N
			tmp = la.dot(nk, @view(A[k:N,j]))
			for i in k:N
				A[i,j] -= 2.0 * n[i] * tmp
			end
		end
	end

	A
end

# ╔═╡ a01dc808-324d-4c00-bcf5-153ad644db15
my_R, my_Q = let R = copy(A0), Q = similar(R)
	N = size(R, 2)

	@inbounds for j in 1:N
		for i in 1:N
			Q[i,j] = 0.0
		end
		@inbounds Q[j,j] = 1.0
	end

	for k in 1:N-1
		scale = maximum(abs, @view(R[k:N,k]))
		@inbounds for i in k:N
			R[i,k] /= scale
		end

		alpha = la.norm(@view(R[k:N,k]))
		R[k,k] += alpha
		alpha_beta = alpha * R[k,k]

		for j in k+1:N
			tmp = la.dot(@view(R[k:N,k]), @view(R[k:N,j]))
			tmp /= alpha_beta
			for i in k:N
				R[i,j] -= tmp * R[i,k]
			end
		end

		for j in 1:N
			tmp = la.dot(@view(R[k:N,k]), @view(Q[k:N,j]))
			tmp /= alpha_beta
			for i in k:N
				Q[i,j] -= tmp * R[i,k]
			end
		end

		R[k,k] = -scale * alpha
	end

	R, Q
end

# ╔═╡ 1429728d-aa2e-4a22-8cca-4529cd3d1a6f
ref_R, ref_Q = let QR = la.qr(copy(A0))
	QR.R, QR.Q
end

# ╔═╡ 763a2e06-d6ad-4a6f-8c3f-e4c2fbcfd273
my_Q * ref_Q

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
# ╠═3515d200-b8ef-4149-8b59-d2dfffeb6ac9
# ╠═a01dc808-324d-4c00-bcf5-153ad644db15
# ╠═1429728d-aa2e-4a22-8cca-4529cd3d1a6f
# ╠═763a2e06-d6ad-4a6f-8c3f-e4c2fbcfd273
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
