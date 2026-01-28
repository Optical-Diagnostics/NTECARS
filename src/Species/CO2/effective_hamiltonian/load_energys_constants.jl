
######################################################################################################
#                                            Structs                                           
######################################################################################################
mutable struct VibRotConstants
    ωi   ::MVector{3, Float64}
    xij  ::MMatrix{3, 3, Float64}
    xll  ::Float64
    yijk ::MArray{Tuple{3, 3, 3}, Float64, 3, 27}
    yill ::MVector{3, Float64}
    zijkm::MArray{Tuple{3, 3, 3, 3}, Float64, 4, 81}
    zijll::MMatrix{3, 3, Float64}
    zllll::Float64
    Be   ::Float64
    αi   ::MVector{3, Float64}
    γij  ::MMatrix{3, 3, Float64}
    γll  ::Float64
    εijk ::MArray{Tuple{3, 3, 3}, Float64, 3, 27}
    εill ::MVector{3, Float64}
    De   ::Float64
    βi   ::MVector{3, Float64}
    ηij  ::MMatrix{3, 3, Float64}
    ηll  ::Float64
    He   ::Float64
    δi   ::MVector{3, Float64}
    Ge   ::Float64
    # Fermi resonance
    Fe   ::Float64
    Fi   ::MVector{3, Float64}
    FJ   ::Float64
    FJJ  ::Float64
    Fij  ::MMatrix{3, 3, Float64}
    FeIV ::Float64
    F1IV ::Float64
    F2IV ::Float64
    F3IV ::Float64
    FJIV ::Float64
    # l-type doubling
    Le  ::Float64
    LJ  ::Float64
    Li  ::MVector{3, Float64}
    Lij ::MMatrix{3, 3, Float64}
    # Coriolis interaction
    Ce  ::Float64
    CJ  ::Float64
    Ci  ::MVector{3, Float64}

    function VibRotConstants()
        return new(
        (@MArray zeros(3)),       #ωi   ::SVector{3, Float64}
        (@MArray zeros(3,3)),     #xij  ::SVector{3, 3, Float64}
        0.0,                      #xll  ::Float64
        (@MArray zeros(3,3,3)),   #yijk ::SVector{3, 3, 3, Float64}
        (@MArray zeros(3)),       #yill ::SVector{3, Float64}
        (@MArray zeros(3,3,3,3)), #zijkm::SVector{3, 3, 3, 3, Float64}
        (@MArray zeros(3,3)),     #zijll::SVector{3, 3, Float64}
        0.0,                      #zllll::Float64
        0.0,                      #Be   ::Float64
        (@MArray zeros(3)),       #αi   ::SVector{3, Float64}
        (@MArray zeros(3,3)),     #γij  ::SVector{3, 3, Float64}
        0.0,                      #γll  ::Float64
        (@MArray zeros(3, 3, 3)), #εijk ::SVector{3, 3, 3, Float64}
        (@MArray zeros(3)),       #εill ::SVector{3, Float64}
        0.0,                      #De   ::Float64
        (@MArray zeros(3)),       #βi   ::SVector{3, Float64}
        (@MArray zeros(3,3)),     #ηij  ::SVector{3, 3, Float64}
        0.0,                      #ηll  ::Float64
        0.0,                      #He   ::Float64
        (@MArray zeros(3)),       #δi   ::SVector{3, Float64}
        0.0,                      #Ge   ::Float64
        0.0,                      #Fe  ::Float64
        (@MArray zeros(3)),       #Fi  ::SVector{3, Float64}
        0.0,                      #FJ  ::Float64
        0.0,                      #FJJ  ::Float64
        (@MArray zeros(3,3)),     #Fij  ::SVector{3, 3, Float64}
        0.0,                      #FeIV::Float64
        0.0,                      #F1IV::Float64
        0.0,                      #F2IV::Float64
        0.0,                      #F3IV::Float64
        0.0,                      #FJIV::Float64
        0.0,                      #Le  ::Float64
        0.0,                      #LJ  ::Float64
        (@MArray zeros(3)),       #Li  ::MVector{3, Float64}
        (@MArray zeros(3,3)),     #Lij ::MMatrix{3, 3, Float64}
        0.0,
        0.0,
        (@MArray zeros(3))
        )
    end
end


######################################################################################################
#                                            functions                                          
######################################################################################################

function toml_string_to_indices(string_indices)
    indices = parse.(Int64, split(string_indices, "")) # "123" -> [1,2,3]
end

function set_values(quantity, constants_dict)
    for key in keys(constants_dict)
        quantity[toml_string_to_indices(key)...] = constants_dict[key]
    end
end

function get_CO2_energy_constants(;iso::Symbol = :O16C12O16)
    if iso != :O16C12O16
        error("Currently, only O16C12O16 is defined at not $(iso)")
    end

    if iso == :O16C12O16
        constants   = VibRotConstants()
        vibrot_data = TOML.parsefile(ROVIBRATIONAL_CONSTANTS_PATH)

        # diagonal terms
        constants.xll = vibrot_data["diagonal_terms"]["x_ll"]
        constants.Be  = vibrot_data["diagonal_terms"]["B_e"]
        constants.γll = vibrot_data["diagonal_terms"]["gamma_ll"]
        constants.De  = vibrot_data["diagonal_terms"]["D_e"]
        constants.ηll = vibrot_data["diagonal_terms"]["eta_ll"]
        constants.He  = vibrot_data["diagonal_terms"]["H_e"]
        set_values(constants.ωi,    vibrot_data["diagonal_terms"]["omega_i"])
        set_values(constants.xij,   vibrot_data["diagonal_terms"]["x_ij"])
        set_values(constants.yijk,  vibrot_data["diagonal_terms"]["y_ijk"])
        set_values(constants.yill,  vibrot_data["diagonal_terms"]["y_ill"])
        set_values(constants.zijkm, vibrot_data["diagonal_terms"]["z_ijkm"])
        set_values(constants.zijll, vibrot_data["diagonal_terms"]["z_ijll"])
        set_values(constants.αi,    vibrot_data["diagonal_terms"]["alpha_i"])
        set_values(constants.γij,   vibrot_data["diagonal_terms"]["gamma_ij"])
        set_values(constants.εijk,  vibrot_data["diagonal_terms"]["epsilon_ijk"])
        set_values(constants.εill,  vibrot_data["diagonal_terms"]["epsilon_ill"])
        set_values(constants.βi,    vibrot_data["diagonal_terms"]["beta_i"])
        set_values(constants.ηij,   vibrot_data["diagonal_terms"]["eta_ij"])

        # l-type doubling
        constants.Le   = vibrot_data["l_type_doubling"]["L_e"]
        constants.LJ   = vibrot_data["l_type_doubling"]["L_J"]
        set_values(constants.Li,  vibrot_data["l_type_doubling"]["L_i"])
        set_values(constants.Lij, vibrot_data["l_type_doubling"]["L_ij"])

        # fermi-interaction terms
        constants.Fe   = vibrot_data["Fermi_interaction_terms"]["F_e"]
        constants.FJ   = vibrot_data["Fermi_interaction_terms"]["F_J"]
        constants.FJJ  = vibrot_data["Fermi_interaction_terms"]["F_JJ"]
        constants.FeIV = vibrot_data["Fermi_interaction_terms"]["FeIV"]
        constants.F1IV = vibrot_data["Fermi_interaction_terms"]["F1IV"]
        constants.F2IV = vibrot_data["Fermi_interaction_terms"]["F2IV"]
        constants.F3IV = vibrot_data["Fermi_interaction_terms"]["F3IV"]
        constants.FJIV = vibrot_data["Fermi_interaction_terms"]["FJIV"]
        set_values(constants.Fi,  vibrot_data["Fermi_interaction_terms"]["F_i"])
        set_values(constants.Fij, vibrot_data["Fermi_interaction_terms"]["F_ij"])

        # Coriolis interaction
        constants.Ce   = vibrot_data["coriolis_interaction"]["C_e"]
        constants.CJ   = vibrot_data["coriolis_interaction"]["C_J"]
        set_values(constants.Ci,  vibrot_data["coriolis_interaction"]["C_i"])

        return constants
    end
end

function rotational_constants(qn::QuantumNumbers)
    # parse string into format of database
    qn_string = "$(qn.v₁)$(qn.v₂)$(qn.l₂)$(qn.v₃)$(qn.r)"
    if qn.l₂ > 0
        qn_string *= (qn.J+qn.l₂+qn.v₃) % 2 == 0 ? "e" : "f"
    end

    # load database
    constants = TOML.parsefile(VIBSTATE_SPECIFIC_ROT_CONSTANTS_PATH)

    # return result from database or default to ground state constants
    constants_of_state = qn_string in keys(constants) ? constants[qn_string] : constants["00001"]
    multipliers        = [1.0, 1e-7, 1e-13]
    return constants_of_state .* multipliers
end
