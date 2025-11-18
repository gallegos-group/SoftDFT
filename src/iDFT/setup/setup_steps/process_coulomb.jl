function process_coulomb(bulk_system::IsingLST, geometry::CoordSystem)
    # Default: "coul" when nothing is specified
    model_name = "coul"

    # Get model type from registry
    ctype = get_coulomb_model(model_name)

    # Construct Coulomb object
    return construct_coulomb(
        ctype,
        bulk_system.molsys,
        bulk_system.bulk,
        geometry
    )
end
