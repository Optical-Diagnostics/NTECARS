struct MoleculeMetadata
    iso_ID        ::Symbol  # ID of the isotopologue
    abundance     ::Real # mol abundance of the isotopologue
    mass          ::Real # mass of the isotopologue
    gᵢ            ::Real  # state independant statistical weight

    function MoleculeMetadata(iso::Symbol = :O16C12O16)
        if iso != :O16C12O16
            error("Currently, only O16C12O16 is defined at not $(iso)")
        end

        if iso == :O16C12O16
            isotoe_data = TOML.parsefile(ISOTOPOLOGUE_METADATA_PATH)
            abundance      = isotoe_data["abundance"]
            mass           = isotoe_data["mass"]*amu
            gᵢ             = isotoe_data["g_i"]
            return new(:O16C12O16, abundance, mass, gᵢ)
        end
    end
end


