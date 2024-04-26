
mutable struct ScattTables
    zKuR::Array{Float64, 1}
    zKaR::Array{Float64, 1}
    dmr::Array{Float64, 1}
    rainRate::Array{Float64, 1}
    attKuR::Array{Float64, 1}
    attKaR::Array{Float64, 1}
    kextR::Array{Float64, 2}
    salbR::Array{Float64, 2}
    asymR::Array{Float64, 2}
    rwc::Array{Float64, 1}
    zKuS::Array{Float64, 1}
    zKaS::Array{Float64, 1}
    dms::Array{Float64, 1}
    snowRate::Array{Float64, 1}
    swc::Array{Float64, 1}
    attKuS::Array{Float64, 1}
    attKaS::Array{Float64, 1}
    kextS::Array{Float64, 2}
    salbS::Array{Float64, 2}
    asymS::Array{Float64, 2}
    zKuBB::Array{Float64, 1}
    zKaBB::Array{Float64, 1}
    dmBB::Array{Float64, 1}
    precRateBB::Array{Float64, 1}
    attKuBB::Array{Float64, 1}
    attKaBB::Array{Float64, 1}
    kextBB::Array{Float64, 2}
    salbBB::Array{Float64, 2}
    asymBB::Array{Float64, 2}
    dmg::Array{Float64, 1}
    zKuG::Array{Float64, 1}
    zKaG::Array{Float64, 1}
    graupRate::Array{Float64, 1}
    gwc::Array{Float64, 1}
    attKuG::Array{Float64, 1}
    attKaG::Array{Float64, 1}
    kextG::Array{Float64, 2}
    salbG::Array{Float64, 2}
    asymG::Array{Float64, 2}
end
function scattTables(pyTables_)
    zKuR=Float64.(pyTables_.zKuR)
    zKaR=Float64.(pyTables_.zKaR)
    dmr=Float64.(pyTables_.dmr)
    rainRate=Float64.(pyTables_.rainRate)
    rwc=Float64.(pyTables_.rwc)
    dms=Float64.(pyTables_.dms)
    snowRate=Float64.(pyTables_.snowRate)
    swc=Float64.(pyTables_.swc)
    dmBB=Float64.(pyTables_.dmBB)
    precRateBB=Float64.(pyTables_.precRateBB)
    gwc=Float64.(pyTables_.gwc)
    graupRate=Float64.(pyTables_.graupRate)
    zKuR=Float64.(pyTables_.zKuR)
    zKaR=Float64.(pyTables_.zKaR)
    zKuS=Float64.(pyTables_.zKuS)
    zKaS=Float64.(pyTables_.zKaS)
    zKuBB=Float64.(pyTables_.zKuBB)
    zKaBB=Float64.(pyTables_.zKaBB)
    zKuG=Float64.(pyTables_.zKuG)
    zKaG=Float64.(pyTables_.zKaG)
    attKuR=Float64.(pyTables_.attKuR)
    attKaR=Float64.(pyTables_.attKaR)
    attKuS=Float64.(pyTables_.attKuS)
    attKaS=Float64.(pyTables_.attKaS)
    attKuBB=Float64.(pyTables_.attKuBB)
    attKaBB=Float64.(pyTables_.attKaBB)
    attKuG=Float64.(pyTables_.attKuG)
    attKaG=Float64.(pyTables_.attKaG)
    kextR=Float64.(pyTables_.kextR)
    salbR=Float64.(pyTables_.salbR)
    asymR=Float64.(pyTables_.asymR)
    kextS=Float64.(pyTables_.kextS)
    salbS=Float64.(pyTables_.salbS)
    asymS=Float64.(pyTables_.asymS)
    kextBB=Float64.(pyTables_.kextBB)
    salbBB=Float64.(pyTables_.salbBB)
    asymBB=Float64.(pyTables_.asymBB)
    kextG=Float64.(pyTables_.kextG)
    salbG=Float64.(pyTables_.salbG)
    asymG=Float64.(pyTables_.asymG)
    dmg=Float64.(pyTables_.dmg)
    return ScattTables(zKuR, zKaR, dmr, rainRate, attKuR, attKaR, kextR, salbR, asymR, rwc, zKuS, zKaS, dms, snowRate, swc, attKuS, attKaS, kextS, salbS, asymS, zKuBB, zKaBB, dmBB, precRateBB, attKuBB, attKaBB, kextBB, salbBB, asymBB, dmg, zKuG, zKaG, graupRate, gwc, attKuG, attKaG, kextG, salbG, asymG)
    
end