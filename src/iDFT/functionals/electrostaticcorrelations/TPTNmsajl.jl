
# struct TPT1msaFunctional{M, P} <: AbstractFunctional{TPT1{MSA_EOS}}

#     weight_fft_K :: Array{Float64, P}

#     Zeff_K :: Array{Float64, P}
#     kappa_K :: Array{Float64, M}

#     dZeff_K :: Array{Float64, P}
#     dkappa_K :: Array{Float64, P}

#     struct_EOS :: TPT1{MSA_EOS}

#     function TPT1msaFunctional(struct_iLST, geometry)

#         DP = struct_iLST.struct_molsys.properties.monomers["diameters"]
#         DR = geometry.bin_width

#         NP = geometry.NP
#         NB = length(DP)
        
#         dims_NB = (NP..., NB)

#         weight_fft_K = zeros(Float64, dims_NB)

#         weight_fft = zeros(Float64, NP)

#         Rsys = CartesianIndices(NP)


#         for j in eachindex(DP)
#             bval = DP[j] / 2.0  # This should actually be related to Gamma from MSA
#             delta_func_fft(weight_fft, DR, bval)
            
#             for K in Rsys
#                 weight_fft_K[K, j] = weight_fft[K]  / 4.0 / pi / bval^2 
#             end 
#         end

        
#         Zeff_K = zeros(Float64, dims_NB)
#         kappa_K = zeros(Float64, NP)

#         dims_derivs = (dims_NB..., NB)
#         dZeff_K = zeros(Float64, dims_derivs)
#         dkappa_K = zeros(Float64, dims_derivs)

#         struct_TPTel = filter(e -> isa(e, TPT1{MSA_EOS}), struct_iLST.struct_bulk.struct_EOS)[1]

#         new{ndims(weight_fft_K)}(weight_fft_K, Zeff_K, kappa_K, dZeff_K, dkappa_K, struct_TPTel)
#     end
# end




# function eval_fe_functional(fields, struct_iLST, struct_functional :: TPT1msaFunctional)

#     rho_hat = fields.ffts.rho_hat

#     planB = fields.ffts.plan_backward
#     f_hat = fields.ffts.f_hat

#     weight_fft_K = struct_functional.weight_fft_K

#     NP = size(f_hat)
#     NB = size(weight_fft_K)[end]
    
#     dims_NB = (NP..., NB)

#     wt_rho_K = zeros(Float64, dims_NB)


#     Rsys = CartesianIndices(NP)

#     # Calculate weighted density

#     for J = 1:NB
#         for K in Rsys
#             f_hat[K] = rho_hat[K,J] * weight_fft_K[K, J]
#         end
        
#         f_hat = planB * f_hat

#         for K in Rsys
#             wt_rho_K[K, J] = abs(real(f_hat[K]))
#         end
#     end


#     struct_molsys = struct_iLST.struct_molsys

#     DP = struct_molsys.properties.monomers["diameters"]
#     ZP = struct_molsys.properties.monomers["valences"]
#     BJ = struct_molsys.properties.system["bjerrum_length"]

#     Zeff_K = struct_functional.Zeff_K1
#     kappa_K = struct_functional.kappa_K

#     fe = 0.0
#     for (u, fixed) in enumerate(fields.fixed)
#         sequence = struct_molsys.configurations[u].sequence
#         MTS = length(sequence)
#         num_chains = length(fixed.configuration)

#         if any(x -> x == true, fields.fixed[u].segments)
#             norm = fields.fixed[u].density[1]/prod(DR)/num_chains
#         else
#             norm = struct_iLST.rho.species[u]/prod(DR)/num_chains
#         end

#         if struct_molsys.configurations[u].evaluator == SimulationEvaluation()
#             coordinates = fixed.coordinates

#             for (n, configuration) in enumerate(fixed.configuration)
#                 coord = coordinates[n]
            
#                 L = 0
#                 for config in eachslice(configuration; dims=2)
#                     L += 1

#                     j = sequence[L] + Int(config[4, L])

#                     K = CartesianIndex(coord[L])

#                     for L1 = L+1:MTS
#                         j1 = sequence[L1] + Int(config[4, L1])

#                         K1 = CartesianIndex(coord[L1])

#                         dij_seg = (DP[j] + DP[j1])/2.0

#                         rij_seg = 0.0
#                         for v in 1:3
#                             rij_seg += (config[v, L] - config[v, L1])^2
#                         end
#                         rij_seg = sqrt(rij_seg)


#                         # segment j
#                         lng_el = -BJ/dist*exp(-kappa_K[K]*(rij_seg-dij_seg))*Zeff_K[K, j]*Zeff_K[K, j1]

#                         fe -= norm*lng_el/2.0


#                         # segment j1
#                         lng_el = -BJ/dist*exp(-kappa_K[K1]*(rij_seg-dij_seg))*Zeff_K[K1, j]*Zeff_K[K1, j1]

#                         fe -= norm*lng_el/2.0
#                     end
#                 end
#             end
#         end
#     end

#     return fe
# end


# function eval_mu_functional(fields, struct_iLST, struct_functional :: TPT1msaFunctional)

#     rho_hat = fields.ffts.rho_hat

#     planB = fields.ffts.plan_backward
#     f_hat = fields.ffts.f_hat

#     weight_fft_K = struct_functional.weight_fft_K

#     NP = size(f_hat)
#     NB = size(weight_fft_K)[end]
    
#     dims_NB = (NP..., NB)
#     wt_rho_K = zeros(Float64, dims_NB)


#     Rsys = CartesianIndices(NP)

#     # Calculate weighted density

#     for J = 1:NB
#         for K in Rsys
#             f_hat[K] = rho_hat[K,J] * weight_fft_K[K, J]
#         end
        
#         f_hat = planB * f_hat

#         for K in Rsys
#             wt_rho_K[K, J] = abs(real(f_hat[K]))
#         end
#     end


#     # Calculate point-wise derivative
#     dPhi_K = zeros(Float64, dims_NB)

#     DP = struct_iLST.struct_molsys.properties.monomers["diameters"]
#     ZP = struct_iLST.struct_molsys.properties.monomers["valences"]
#     BJ = struct_iLST.struct_molsys.properties.system["bjerrum_length"]

#     struct_EOS = struct_functional.struct_EOS
    
#     wt_rho = struct_EOS.wt_rho
#     dZeff = struct_EOS.dZeff
#     dkappa = struct_EOS.dkappa
#     Zeff = struct_EOS.Zeff

#     for K in Rsys      
#         @views @. wt_rho = wt_rho_K[K, :]

#         kappa = derivFE_TPTNmsa(dZeff, dkappa, Zeff, wt_rho, DP, ZP, BJ)

#         @views @. Zeff_K[K, :] = Zeff
#         @views @. kappa_K[K] = kappa

#         @views @. dZeff_K[K, :, :] = dZeff
#         @views @. dkappa_K[K, :, :] = dkappa
#     end


#     # Calculate dPhi_K based off configurations

#     for (u, fixed) in enumerate(fields.fixed)
#         sequence = struct_molsys.configurations[u].sequence
#         MTS = length(sequence)
#         num_chains = length(fixed.configuration)

#         if any(x -> x == true, fields.fixed[u].segments)
#             norm = fields.fixed[u].density[1]/prod(DR)/num_chains
#         else
#             norm = struct_iLST.rho.species[u]/prod(DR)/num_chains
#         end

#         if struct_molsys.configurations[u].evaluator == SimulationEvaluation()
#             coordinates = fixed.coordinates

#             for (n, configuration) in enumerate(fixed.configuration)
#                 coord = coordinates[n]
            
#                 L = 0
#                 for config in eachslice(configuration; dims=2)
#                     L += 1

#                     j = sequence[L] + Int(config[4, L])

#                     K = CartesianIndex(coord[L])

#                     for L1 = L+1:MTS
#                         j1 = sequence[L1] + Int(config[4, L1])

#                         K1 = CartesianIndex(coord[L1])

#                         dij_seg = (DP[j] + DP[j1])/2.0

#                         rij_seg = 0.0
#                         for v in 1:3
#                             rij_seg += (config[v, L] - config[v, L1])^2
#                         end
#                         rij_seg = sqrt(rij_seg)

#                         for M in 1:NB
#                             dlny_el = -BJ/rij_seg*exp(-kappa_K[K]*(rij_seg-dij_seg)) *
#                                 (dZeff_K[K, j, M]*Zeff_K[K, j1]+Zeff_K[K, j]*dZeff_K[K, j1, M] 
#                                     + Zeff_K[K, j]*Zeff_K[K, j1]*dkappa_K[K, M]*(rij_seg-dij_seg))

#                             dPhi_K[K, M] -= norm*dlny_el/2.0

#                             dlny_el = -BJ/rij_seg*exp(-kappa_K[K1]*(rij_seg-dij_seg)) *
#                                 (dZeff_K[K1, j, M]*Zeff_K[K1, j1]+Zeff_K[K1, j]*dZeff_K[K1, j1, M] 
#                                     + Zeff_K[K1, j]*Zeff_K[K1, j1]*dkappa_K[K1, M]*(rij_seg-dij_seg))

#                             dPhi_K[K1, M] -= norm*dlny_el/2.0
#                         end
#                     end
#                 end
#             end
#         end
#     end



#     # Calculate excess chemical contribution by weighting point-wise derivative    
#     planF = fields.ffts.plan_forward
#     mu_ex_hat = fields.ffts.mu_ex_hat


#     # This is fine for one-body field
#     for J = 1:NB
#         for K in Rsys
#             f_hat[K] = dPhi_K[K, J]
#         end
#         f_hat = planF * f_hat
        
#         for K in Rsys
#             f_hat[K] *= weight_fft_K[K, J]
#         end

#         for K in Rsys
#             mu_ex_hat[K, J] += f_hat[K]
#         end
#     end
    





#     # This has to be redone
#     for u = 1:MB
#         j, j1 = struct_iLST.bonds[u]
        
#         for K in Rsys
#             f_hat[K] = lng_K[K, u]
#         end
#         f_hat = planF * f_hat
        
#         for K in Rsys
#             lng_hat[K, 1, u] += f_hat[K] * weight_fft_K[K, j]
#             lng_hat[K, 2, u] += f_hat[K] * weight_fft_K[K, j1]
#         end
#     end
# end




# function calc_2body_field(u, IB, conform, state, struct_scmf, struct_functional :: TPTNmsaFunctional)

#     Zeff_K = struct_functional.Zeff_K
#     kappa_K = struct_functional.Zeff_K
    
#     struct_iLST = struct_scmf.struct_iDFT.struct_iLST
#     DP = struct_iLST.struct_molsys.properties.monomers["diameters"]
#     BJ = struct_iLST.struct_molsys.properties.system["bjerrum_length"]

#     sequence = struct_iLST.struct_molsys.configurations[u].sequence
#     MTS = length(sequence)

#     idx = ()
#     for v in 1:ndim-1
#         idx = (floor(Int, conform[4-v, IB] / bin_width[ndim-v]) + 1, idx..., )
#     end

#     K = CartesianIndex(idx)

#     j = sequence[IB] + state[IB]

#     E_LR = 0.0
#     for L1 in 1:IB-1
#         j1 = sequence[L1] + state[L1]

#         idx1 = ()
#         for v in 1:ndim-1
#             idx1 = (floor(Int, conform[4-v, L1] / bin_width[ndim-v]) + 1, idx1..., )
#         end
    
#         K1 = CartesianIndex(idx1)

#         dist = 0.0
#         for v in 1:3
#             dist += (conform[v, IB] - conform[v, L1])^2
#         end
#         dist = sqrt(dist)

#         dij_seg = (DP[IB] + DP[L1])/2.0

#         lng = BJ/dist*exp(-kappa_K[K]*(dist-dij_seg))*Zeff_K[K, j]*Zeff_K[K, j1]
#         lng1 = BJ/dist*exp(-kappa_K[K1]*(dist-dij_seg))*Zeff_K[K1, j]*Zeff_K[K1, j1]

#         E_LR = (lng + lng1)/2.0
#     end

    
#     for L1 in IB+1:MTS
#         j1 = sequence[L1] + state[L1]

#         idx1 = ()
#         for v in 1:ndim-1
#             idx1 = (floor(Int, conform[4-v, L1] / bin_width[ndim-v]) + 1, idx1..., )
#         end
    
#         K1 = CartesianIndex(idx1)

#         dist = 0.0
#         for v in 1:3
#             dist += (conform[v, IB] - conform[v, L1])^2
#         end
#         dist = sqrt(dist)

#         dij_seg = (DP[IB] + DP[L1])/2.0

#         lng = BJ/dist*exp(-kappa_K[K]*(dist-dij_seg))*Zeff_K[K, j]*Zeff_K[K, j1]
#         lng1 = BJ/dist*exp(-kappa_K[K1]*(dist-dij_seg))*Zeff_K[K1, j]*Zeff_K[K1, j1]

#         E_LR = (lng + lng1)/2.0
#     end

#     return E_LR
# end