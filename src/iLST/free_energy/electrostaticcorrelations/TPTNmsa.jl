# struct TPTNmsa_Structure <: AbstractFreeEnergy
#     dPhi :: Vector{Float64}
#     Zeff :: Vector{Float64}
#     dZeff :: Matrix{Float64}
#     dkappa :: Vector{Float64}

#     wt_rho :: Vector{Float64}

#     function TPTNmsa_Structure(rho)
#         rho_beads = rho.beads

#         num = length(rho_beads)

#         dPhi = zeros(Float64, num)
#         Zeff = zeros(Float64, num)
#         dZeff = zeros(Float64, (num, num))
#         dkappa = zeros(Float64, num)

#         wt_rho = zeros(Float64, num)

#         new(dPhi, Zeff, dZeff, dkappa, wt_rho)
#     end
# end


# function free_energy(struct_iLST, struct_TPTNmsa :: TPTNmsa_Structure)


    
#     # Not meant for bulk at this time






    
#     # properties = struct_iLST.molsys.properties

#     # rho = struct_iLST.rho.beads

#     # DP = properties.monomers["diameters"]
#     # ZP = properties.monomers["valences"]

#     # BJ = properties.system["bjerrum_length"]

#     # Zeff = struct_EOS.Zeff

#     # f_ch = free_energy_TPTNmsa(Zeff, rho, DP, ZP, BJ)

#     # return f_ch
# end


# function chemical_potential(struct_iLST, struct_TPTNmsa :: TPTNmsa_Structure)

#     # Not meant for bulk at this time







#     # properties = struct_iLST.molsys.properties   

#     # rho = struct_iLST.rho.beads

#     # DP = properties.monomers["diameters"]
#     # ZP = properties.monomers["valences"]

#     # BJ = properties.system["bjerrum_length"]

#     # dPhi = struct_TPTNmsa.dPhi

#     # kappa = derivFE_TPTNmsa(dZeff, dkappa, Zeff, rho, DP, ZP, BJ)

#     # @views @. struct_iLST.bulk.mu_ex += dPhi
# end

# # Inhomogeneous and bulk term

# function free_energy_TPTNmsa(Zeff, rho, DP, ZP, BJ)

#     Gam, eta, _ = Gamma_MSA(rho, DP, ZP, BJ)

#     kappa = sqrt(4.0*pi*BJ*sum(@~ @. rho*ZP*ZP))

#     for J in eachindex(rho)
#         Zeff[J] = (ZP[J]-eta*DP[J]^2)/(1.0+Gam*DP[J])
#     end

#     return kappa
# end

# function derivFE_TPTNmsa(dZeff, dkappa, Zeff, rho, DP, ZP, BJ)

#     @. dZeff = 0.0
#     @. dkappa = 0.0
#     @. Zeff = ZP

#     Gam, eta, Hgam = Gamma_MSA(rho, DP, ZP, BJ)

#     derivFE_MSA(dPhi, rho, DP, ZP, BJ)

#     kappa = sqrt(4.0*pi*BJ*sum(@~ @. rho*ZP*ZP))
    
#     if kappa > 1e-10
#         @. dkappa = 2.0*pi*BJ*ZP^2/kappa
#     else
#         @. dkappa = 0.0
#     end

#     pdf_gam = -2.0*pi*BJ*
#         sum(@~ @. rho*(ZP-eta*DP*DP)/(1.0+Gam*DP)*(ZP-eta*DP*DP)/(1.0+Gam*DP)/(1.0+Gam*DP)*DP)

#     pdf_eta = -2.0*pi*BJ*
#         sum(@~ @. rho*DP*DP*(ZP-eta*DP*DP)/(1.0+Gam*DP)/(1.0+Gam*DP))

#     pdeta_gam = -1.0/Hgam*sum(@~ @. rho*ZP*DP*DP/(1.0+Gam*DP)/(1.0+Gam*DP))

#     pdeta_Hgam = -eta/Hgam

#     pdHgam_gam = -6.0/pi*sum(@~ @. pi/6.0*rho*DP*DP*DP*DP/(1.0+Gam*DP)/(1.0+Gam*DP))

#     deta_gam = pdeta_gam + pdeta_Hgam*pdHgam_gam

#     denom = 1.0 / (2.0*Gam-pdf_eta*deta_gam-pdf_gam)
#     for W in eachindex(rho)
#         inv_1_plus_GamDP = 1.0 / (1.0 + Gam*DP[W])
#         ZeffW = (ZP[W]-eta*DP[W]^2)*inv_1_plus_GamDP

#         pdf_rho = pi*BJ*ZeffW^2
#         pdeta_rho = DP[W]/Hgam*ZP[W]*inv_1_plus_GamDP
#         pdHgam_rho = DP[W]^3*(inv_1_plus_GamDP-1.0/3.0)

#         dGam_rho = (pdf_rho+pdf_eta*(pdeta_rho + pdf_eta*pdeta_Hgam*pdHgam_rho))*denom
#         deta_rho = deta_gam*dGam_rho

#         for J in eachindex(rho)
#             # J terms
#             inv_1_plus_GamDPJ = 1.0/(1.0+Gam*DP[J])

#             Zeff[J] = (ZP[J]-eta*DP[J]^2)*inv_1_plus_GamDPJ

#             pdZeffJ_eta = DP[J]^2 * inv_1_plus_GamDPJ
#             pdZeffJ_Gam = DP[J] * Zeff[J] * inv_1_plus_GamDPJ

#             dZeff[J, W] = -(pdZeffJ_eta * deta_rho + pdZeffJ_Gam * dGam_rho)
#         end
#     end
    
#     return kappa
# end